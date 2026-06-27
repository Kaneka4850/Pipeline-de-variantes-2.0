# ============================================================
# Pipeline de Variantes Germinativas — Environment Capture
# ============================================================
"""
Captura completa do ambiente de execução para reprodutibilidade.

Registra:
- Sistema operacional, kernel, hostname, usuário
- Hardware: CPUs, memória RAM
- Versões de todas as ferramentas do pipeline
- Informações do VEP (versão, cache, assembly, plugins)
- Configuração utilizada
- Timestamp da captura

O resultado é um EnvironmentSnapshot serializável para JSON,
armazenável no SQLite e incluído no relatório de execução.
"""

import json
import os
import platform
import re
import shutil
import subprocess
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Any

from .utils import get_tool_version, get_hostname, get_username


# ============================================================
# DATA CLASSES
# ============================================================

@dataclass
class ToolInfo:
    """Informação de uma ferramenta."""
    name: str
    version: str
    path: str
    available: bool

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class VepInfo:
    """Informações detalhadas do Ensembl VEP."""
    version: str = ""
    cache_version: str = ""
    cache_dir: str = ""
    assembly: str = ""
    fasta_version: str = ""
    plugins_dir: str = ""
    plugins_installed: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class SystemInfo:
    """Informações do sistema operacional e hardware."""
    os_name: str = ""
    os_version: str = ""
    os_release: str = ""
    kernel: str = ""
    architecture: str = ""
    hostname: str = ""
    username: str = ""
    cpu_count: int = 0
    cpu_model: str = ""
    memory_total_mb: int = 0
    memory_available_mb: int = 0

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class EnvironmentSnapshot:
    """Snapshot completo do ambiente de execução."""
    timestamp: str = ""
    pipeline_version: str = ""
    system: SystemInfo = field(default_factory=SystemInfo)
    tools: dict[str, ToolInfo] = field(default_factory=dict)
    vep: VepInfo = field(default_factory=VepInfo)
    config_dump: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        result = {
            "timestamp": self.timestamp,
            "pipeline_version": self.pipeline_version,
            "system": self.system.to_dict(),
            "tools": {k: v.to_dict() for k, v in self.tools.items()},
            "vep": self.vep.to_dict(),
            "config": self.config_dump,
        }
        return result

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent, ensure_ascii=False)


# ============================================================
# CAPTURA DO SISTEMA
# ============================================================

def capture_system_info() -> SystemInfo:
    """Captura informações do sistema operacional e hardware."""
    info = SystemInfo(
        os_name=platform.system(),
        os_version=platform.version(),
        os_release=platform.release(),
        kernel=platform.platform(),
        architecture=platform.machine(),
        hostname=get_hostname(),
        username=get_username(),
        cpu_count=os.cpu_count() or 0,
        cpu_model=_get_cpu_model(),
        memory_total_mb=_get_memory_total_mb(),
        memory_available_mb=_get_memory_available_mb(),
    )
    return info


def _get_cpu_model() -> str:
    """Detecta o modelo da CPU."""
    try:
        if platform.system() == "Linux":
            with open("/proc/cpuinfo") as f:
                for line in f:
                    if line.startswith("model name"):
                        return line.split(":", 1)[1].strip()
        elif platform.system() == "Darwin":
            result = subprocess.run(
                ["sysctl", "-n", "machdep.cpu.brand_string"],
                capture_output=True, text=True, timeout=5,
            )
            return result.stdout.strip()
    except (OSError, subprocess.SubprocessError):
        pass
    return platform.processor() or "desconhecido"


def _get_memory_total_mb() -> int:
    """Detecta memória RAM total em MB."""
    try:
        if platform.system() == "Linux":
            with open("/proc/meminfo") as f:
                for line in f:
                    if line.startswith("MemTotal:"):
                        kb = int(line.split()[1])
                        return kb // 1024
        elif platform.system() == "Darwin":
            result = subprocess.run(
                ["sysctl", "-n", "hw.memsize"],
                capture_output=True, text=True, timeout=5,
            )
            return int(result.stdout.strip()) // (1024 * 1024)
    except (OSError, subprocess.SubprocessError, ValueError):
        pass
    return 0


def _get_memory_available_mb() -> int:
    """Detecta memória RAM disponível em MB."""
    try:
        if platform.system() == "Linux":
            with open("/proc/meminfo") as f:
                for line in f:
                    if line.startswith("MemAvailable:"):
                        kb = int(line.split()[1])
                        return kb // 1024
    except (OSError, ValueError):
        pass
    return 0


# ============================================================
# CAPTURA DE FERRAMENTAS
# ============================================================

# Ferramentas do pipeline e suas estratégias de detecção de versão
PIPELINE_TOOLS = [
    "python3",
    "bash",
    "bwa",
    "gatk",
    "samtools",
    "bcftools",
    "fastqc",
    "multiqc",
    "vep",
    "snakemake",
    "docker",
    "singularity",
    "conda",
]


def capture_tool_versions(
    config: dict[str, Any] | None = None,
) -> dict[str, ToolInfo]:
    """Captura versões de todas as ferramentas do pipeline.

    Args:
        config: Configuração do pipeline (para caminhos customizados).

    Returns:
        Dicionário {nome: ToolInfo}.
    """
    tools_config = (config or {}).get("tools", {})
    result: dict[str, ToolInfo] = {}

    for tool_name in PIPELINE_TOOLS:
        # Usa caminho do config se disponível
        configured_path = tools_config.get(tool_name, "")
        tool_path = configured_path or tool_name

        # Verifica disponibilidade
        found = shutil.which(tool_path)
        if found:
            version = get_tool_version(found)
            result[tool_name] = ToolInfo(
                name=tool_name,
                version=version,
                path=found,
                available=True,
            )
        else:
            result[tool_name] = ToolInfo(
                name=tool_name,
                version="não instalado",
                path="",
                available=False,
            )

    return result


# ============================================================
# CAPTURA DO VEP
# ============================================================

def capture_vep_info(config: dict[str, Any] | None = None) -> VepInfo:
    """Captura informações detalhadas do Ensembl VEP.

    Args:
        config: Configuração do pipeline (seção 'vep').

    Returns:
        VepInfo com dados do VEP.
    """
    vep_config = (config or {}).get("vep", {})
    vep_path = (config or {}).get("tools", {}).get("vep", "") or "vep"

    info = VepInfo(
        cache_dir=vep_config.get("cache_dir", ""),
        assembly=vep_config.get("assembly", ""),
        plugins_dir=vep_config.get("plugins_dir", ""),
    )

    # Versão do VEP
    found = shutil.which(vep_path)
    if found:
        version_str = get_tool_version(found)
        # Extrai número de versão
        match = re.search(r"(\d+)", version_str)
        if match:
            info.version = match.group(0)
        else:
            info.version = version_str

    # Versão do cache (inspeciona diretórios)
    if info.cache_dir and Path(info.cache_dir).is_dir():
        info.cache_version = _detect_vep_cache_version(
            info.cache_dir, info.assembly
        )

    # FASTA de referência do VEP
    fasta_ref = vep_config.get("fasta", "")
    if not fasta_ref:
        fasta_ref = (config or {}).get("reference", {}).get("genome", "")
    info.fasta_version = Path(fasta_ref).name if fasta_ref else ""

    # Plugins instalados
    if info.plugins_dir and Path(info.plugins_dir).is_dir():
        info.plugins_installed = _detect_vep_plugins(info.plugins_dir)

    return info


def _detect_vep_cache_version(cache_dir: str, assembly: str) -> str:
    """Detecta a versão do cache do VEP inspecionando diretórios.

    O cache do VEP segue a estrutura:
        cache_dir/homo_sapiens/112_GRCh38/

    Args:
        cache_dir: Diretório do cache.
        assembly: Assembly esperado (GRCh38, GRCh37).

    Returns:
        Versão detectada (ex: "112") ou string vazia.
    """
    cache_path = Path(cache_dir)

    # Procura em homo_sapiens/ e homo_sapiens_merged/ e homo_sapiens_refseq/
    for species_dir_name in [
        "homo_sapiens", "homo_sapiens_merged", "homo_sapiens_refseq"
    ]:
        species_dir = cache_path / species_dir_name
        if species_dir.is_dir():
            for subdir in sorted(species_dir.iterdir(), reverse=True):
                if subdir.is_dir() and assembly in subdir.name:
                    # Extrai versão do nome (ex: "112_GRCh38" → "112")
                    match = re.match(r"(\d+)_", subdir.name)
                    if match:
                        return match.group(1)

    return ""


def _detect_vep_plugins(plugins_dir: str) -> list[str]:
    """Lista plugins VEP instalados no diretório.

    Args:
        plugins_dir: Diretório dos plugins.

    Returns:
        Lista de nomes de plugins (ex: ["ClinVar", "REVEL"]).
    """
    plugins_path = Path(plugins_dir)
    plugins: list[str] = []

    for f in plugins_path.iterdir():
        if f.is_file() and f.suffix == ".pm":
            plugins.append(f.stem)

    return sorted(plugins)


# ============================================================
# CAPTURA COMPLETA
# ============================================================

def capture_environment(config: dict[str, Any] | None = None) -> EnvironmentSnapshot:
    """Captura snapshot completo do ambiente de execução.

    Args:
        config: Configuração completa do pipeline.

    Returns:
        EnvironmentSnapshot com todos os dados.
    """
    pipeline_config = config or {}
    pipeline_version = pipeline_config.get("pipeline", {}).get("version", "")

    snapshot = EnvironmentSnapshot(
        timestamp=datetime.now().isoformat(),
        pipeline_version=pipeline_version,
        system=capture_system_info(),
        tools=capture_tool_versions(config),
        vep=capture_vep_info(config),
        config_dump=_sanitize_config(pipeline_config),
    )

    return snapshot


def _sanitize_config(config: dict[str, Any]) -> dict[str, Any]:
    """Sanitiza a configuração para armazenamento.

    Remove informações sensíveis e garante serialização JSON.
    """
    # Por enquanto, apenas retorna uma cópia (sem dados sensíveis neste pipeline)
    try:
        return json.loads(json.dumps(config, default=str))
    except (TypeError, ValueError):
        return {"error": "Config não serializável"}
