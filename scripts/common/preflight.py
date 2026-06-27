# ============================================================
# Pipeline de Variantes Germinativas — Preflight Validation
# ============================================================
"""
Validações pre-flight abrangentes antes da execução do pipeline.

Verifica:
1. Referência genômica (FASTA + todos os índices)
2. Cache do VEP (diretório, assembly, arquivos)
3. FASTQs de entrada (existência, integridade gzip, formato)
4. Espaço em disco (estimativa baseada nos FASTQs)
5. Permissões de escrita (diretórios de saída)
6. Ferramentas (disponibilidade no PATH)
7. Consistência do config (assembly VEP = assembly referência)

Cada verificação retorna um PreflightResult com mensagem
de erro e orientação para correção (remediation).
"""

import os
import shutil
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any

from .validation import (
    validate_fastq,
    check_disk_space,
    check_write_permission,
    check_tool_installed,
    estimate_required_space,
)
from .reference_manager import ReferenceManager


# ============================================================
# DATA CLASSES
# ============================================================

@dataclass
class PreflightResult:
    """Resultado de uma verificação pre-flight."""
    check_name: str
    status: str         # "PASS", "WARN", "FAIL"
    message: str
    remediation: str    # Orientação para corrigir

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


# ============================================================
# PREFLIGHT RUNNER
# ============================================================

class PreflightRunner:
    """Orquestra todas as verificações pre-flight.

    Uso:
        runner = PreflightRunner(config, samples_dict)
        results = runner.run_all()
        if runner.has_failures:
            for r in runner.failures:
                print(f"FAIL: {r.message}")
                print(f"  → {r.remediation}")
    """

    def __init__(
        self,
        config: dict[str, Any],
        samples_dict: dict[str, dict[str, str]],
    ) -> None:
        self.config = config
        self.samples = samples_dict
        self.results: list[PreflightResult] = []

    def run_all(self) -> list[PreflightResult]:
        """Executa todas as verificações configuradas.

        Returns:
            Lista de PreflightResult.
        """
        self.results = []
        preflight_config = self.config.get("preflight", {})

        # Verificações sempre executadas
        self._check_tools()
        self._check_config_consistency()

        # Verificações configuráveis
        if preflight_config.get("validate_reference", True):
            self._check_reference()

        if preflight_config.get("validate_fastq", True):
            self._check_fastqs()

        if preflight_config.get("validate_vep_cache", True):
            self._check_vep_cache()

        # Disco e permissões
        self._check_disk_space()
        self._check_permissions()

        return self.results

    @property
    def has_failures(self) -> bool:
        """Verifica se há falhas críticas."""
        return any(r.status == "FAIL" for r in self.results)

    @property
    def failures(self) -> list[PreflightResult]:
        """Retorna apenas os resultados com falha."""
        return [r for r in self.results if r.status == "FAIL"]

    @property
    def warnings(self) -> list[PreflightResult]:
        """Retorna apenas os avisos."""
        return [r for r in self.results if r.status == "WARN"]

    @property
    def passes(self) -> list[PreflightResult]:
        """Retorna verificações que passaram."""
        return [r for r in self.results if r.status == "PASS"]

    def summary(self) -> str:
        """Retorna resumo formatado dos resultados."""
        total = len(self.results)
        passed = len(self.passes)
        warned = len(self.warnings)
        failed = len(self.failures)

        lines = [
            f"Preflight: {total} verificações — "
            f"✅ {passed} OK, ⚠️ {warned} avisos, ❌ {failed} falhas",
        ]

        for r in self.results:
            icon = {"PASS": "✅", "WARN": "⚠️", "FAIL": "❌"}[r.status]
            lines.append(f"  {icon} [{r.check_name}] {r.message}")
            if r.status != "PASS" and r.remediation:
                lines.append(f"     → {r.remediation}")

        return "\n".join(lines)

    # ---- Verificações individuais ----

    def _add(self, check_name: str, status: str, message: str, remediation: str = "") -> None:
        """Adiciona um resultado."""
        self.results.append(PreflightResult(
            check_name=check_name,
            status=status,
            message=message,
            remediation=remediation,
        ))

    def _check_reference(self) -> None:
        """Verifica referência genômica e índices."""
        genome = self.config.get("reference", {}).get("genome", "")

        if not genome:
            self._add(
                "REFERENCE", "FAIL",
                "Genoma de referência não configurado",
                "Defina 'reference.genome' no config.yaml",
            )
            return

        manager = ReferenceManager(genome, self.config)

        if not Path(genome).is_file():
            self._add(
                "REFERENCE_FASTA", "FAIL",
                f"FASTA não encontrado: {genome}",
                f"Baixe o genoma e coloque em: {genome}",
            )
            return

        self._add(
            "REFERENCE_FASTA", "PASS",
            f"FASTA encontrado: {genome}",
        )

        # Verifica índices
        info = manager.collect_info()
        missing = info.missing_indices
        if missing:
            self._add(
                "REFERENCE_INDICES", "WARN",
                f"Índices ausentes: {', '.join(missing)}",
                "Os índices serão criados automaticamente pelo pipeline. "
                "Ou execute: snakemake --until bwa_index samtools_faidx gatk_dict",
            )
        else:
            self._add(
                "REFERENCE_INDICES", "PASS",
                "Todos os índices presentes",
            )

    def _check_fastqs(self) -> None:
        """Verifica integridade dos FASTQs."""
        if not self.samples:
            self._add(
                "FASTQ", "FAIL",
                "Nenhuma amostra detectada",
                f"Coloque os FASTQs em: {self.config.get('dirs', {}).get('samples', 'samples/')}/",
            )
            return

        all_valid = True
        for sample_name, sample_info in self.samples.items():
            for key in ["r1", "r2"]:
                fq_path = sample_info.get(key, "")
                if not fq_path:
                    continue

                valid, msg = validate_fastq(fq_path)
                if not valid:
                    self._add(
                        f"FASTQ_{sample_name}_{key.upper()}", "FAIL",
                        f"[{sample_name}] {msg}",
                        "Verifique se o arquivo está completo e não corrompido",
                    )
                    all_valid = False

        if all_valid:
            total_files = sum(
                1 for s in self.samples.values()
                for k in ["r1", "r2"] if s.get(k)
            )
            self._add(
                "FASTQ", "PASS",
                f"Todos os {total_files} FASTQs validados com sucesso",
            )

    def _check_vep_cache(self) -> None:
        """Verifica cache do VEP."""
        annotation_enabled = self.config.get("pipeline", {}).get("annotation", True)
        vep_path = self.config.get("tools", {}).get("vep", "")

        if not annotation_enabled:
            self._add(
                "VEP_CACHE", "PASS",
                "Anotação VEP desabilitada — verificação ignorada",
            )
            return

        # Verifica se o VEP está instalado
        if not shutil.which(vep_path or "vep"):
            self._add(
                "VEP_INSTALL", "WARN",
                "VEP não encontrado no PATH",
                "Instale o VEP ou desabilite a anotação em config.yaml "
                "(pipeline.annotation: false)",
            )
            return

        cache_dir = self.config.get("vep", {}).get("cache_dir", "")
        assembly = self.config.get("vep", {}).get("assembly", "GRCh38")

        if not cache_dir:
            self._add(
                "VEP_CACHE", "WARN",
                "Diretório do cache VEP não configurado",
                "Defina 'vep.cache_dir' no config.yaml",
            )
            return

        if not Path(cache_dir).is_dir():
            self._add(
                "VEP_CACHE", "FAIL",
                f"Diretório do cache VEP não encontrado: {cache_dir}",
                f"Crie o diretório e instale o cache: "
                f"vep_install -a cf -s homo_sapiens -y {assembly} -c {cache_dir}",
            )
            return

        # Verifica se há cache para o assembly correto
        has_assembly = False
        for species_dir in ["homo_sapiens", "homo_sapiens_merged", "homo_sapiens_refseq"]:
            species_path = Path(cache_dir) / species_dir
            if species_path.is_dir():
                for subdir in species_path.iterdir():
                    if subdir.is_dir() and assembly in subdir.name:
                        has_assembly = True
                        break

        if has_assembly:
            self._add(
                "VEP_CACHE", "PASS",
                f"Cache VEP encontrado para assembly {assembly}",
            )
        else:
            self._add(
                "VEP_CACHE", "WARN",
                f"Cache VEP para {assembly} não encontrado em {cache_dir}",
                f"Instale o cache: vep_install -a cf -s homo_sapiens -y {assembly} -c {cache_dir}",
            )

    def _check_disk_space(self) -> None:
        """Verifica espaço em disco."""
        results_dir = self.config.get("dirs", {}).get("results", "results")
        min_gb = self.config.get("preflight", {}).get("min_disk_gb", 10)

        # Estima espaço necessário
        all_fastqs = [
            s.get(k, "")
            for s in self.samples.values()
            for k in ["r1", "r2"]
            if s.get(k)
        ]
        estimated_gb = estimate_required_space(all_fastqs)
        required_gb = max(estimated_gb, min_gb)

        # Verifica
        ok, msg = check_disk_space(
            results_dir if os.path.isdir(results_dir) else ".",
            required_gb,
        )

        if ok:
            self._add("DISK_SPACE", "PASS", msg)
        else:
            self._add(
                "DISK_SPACE", "WARN",
                msg,
                "Libere espaço em disco ou altere 'dirs.results' para outro volume",
            )

    def _check_permissions(self) -> None:
        """Verifica permissões de escrita."""
        dirs_to_check = ["results", "logs", "database"]
        all_ok = True

        for dir_key in dirs_to_check:
            dir_path = self.config.get("dirs", {}).get(dir_key, "")
            if dir_path:
                ok, msg = check_write_permission(dir_path)
                if not ok:
                    self._add(
                        f"PERMISSIONS_{dir_key.upper()}", "FAIL",
                        msg,
                        f"Execute: chmod 755 {dir_path} ou mude o proprietário",
                    )
                    all_ok = False

        if all_ok:
            self._add(
                "PERMISSIONS", "PASS",
                "Permissões de escrita OK em todos os diretórios",
            )

    def _check_tools(self) -> None:
        """Verifica disponibilidade de ferramentas obrigatórias."""
        required = ["bwa", "gatk", "samtools", "fastqc"]
        optional = ["bcftools", "multiqc"]

        all_required_ok = True
        for tool in required:
            tool_path = self.config.get("tools", {}).get(tool, "")
            ok, msg = check_tool_installed(tool, tool_path)
            if not ok:
                self._add(
                    f"TOOL_{tool.upper()}", "FAIL",
                    f"Ferramenta obrigatória não encontrada: {tool}",
                    f"Instale: conda install -c bioconda {tool}",
                )
                all_required_ok = False

        for tool in optional:
            tool_path = self.config.get("tools", {}).get(tool, "")
            ok, msg = check_tool_installed(tool, tool_path)
            if not ok:
                self._add(
                    f"TOOL_{tool.upper()}", "WARN",
                    f"Ferramenta opcional não encontrada: {tool}",
                    f"Instale: conda install -c bioconda {tool}",
                )

        if all_required_ok:
            self._add(
                "TOOLS_REQUIRED", "PASS",
                "Todas as ferramentas obrigatórias disponíveis",
            )

    def _check_config_consistency(self) -> None:
        """Verifica consistência interna da configuração."""
        ref_name = self.config.get("reference", {}).get("name", "")
        vep_assembly = self.config.get("vep", {}).get("assembly", "")

        # Verifica se assembly no VEP é consistente com referência
        if ref_name and vep_assembly:
            if ref_name.upper() not in vep_assembly.upper() and vep_assembly.upper() not in ref_name.upper():
                self._add(
                    "CONFIG_CONSISTENCY", "WARN",
                    f"Assembly da referência ({ref_name}) pode não corresponder "
                    f"ao assembly do VEP ({vep_assembly})",
                    "Verifique se 'reference.name' e 'vep.assembly' são compatíveis",
                )
                return

        self._add(
            "CONFIG_CONSISTENCY", "PASS",
            "Configuração internamente consistente",
        )
