# ============================================================
# Pipeline de Variantes Germinativas — Reference Manager
# ============================================================
"""
Sistema de gerenciamento e versionamento de referências genômicas.

Coleta, valida e registra metadados completos da referência:
- Nome, versão, organismo, assembly, fonte
- Checksums SHA256 (com cache para evitar recálculo)
- Status de todos os índices (BWA, samtools, GATK)
- Integridade geral

O checksum SHA256 é cacheado em references/.reference_meta.json
para evitar recalcular ~3 GB a cada execução.
"""

import hashlib
import json
import re
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Any


# ============================================================
# DATA CLASSES
# ============================================================

@dataclass
class IndexStatus:
    """Status de um índice da referência."""
    name: str
    extension: str
    path: str
    exists: bool
    size_bytes: int = 0

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class ReferenceInfo:
    """Metadados completos de uma referência genômica."""
    name: str                                    # "GRCh38", "hg19"
    version: str                                 # "p14", "release-110"
    organism: str                                # "Homo sapiens"
    assembly: str                                # "GRCh38" / "GRCh37"
    source: str                                  # "Ensembl", "NCBI", "UCSC"
    fasta_path: str                              # Caminho completo
    fasta_size_bytes: int = 0                    # Tamanho do FASTA
    checksum_sha256: str = ""                    # SHA256 do FASTA
    install_date: str = ""                       # Data de instalação
    indices: dict[str, IndexStatus] = field(default_factory=dict)
    integrity_ok: bool = True                    # Todos os requisitos OK

    def to_dict(self) -> dict[str, Any]:
        result = asdict(self)
        result["indices"] = {
            k: v.to_dict() if isinstance(v, IndexStatus) else v
            for k, v in self.indices.items()
        }
        return result

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), indent=2, ensure_ascii=False)

    @property
    def all_indices_present(self) -> bool:
        """Verifica se todos os índices obrigatórios existem."""
        return all(idx.exists for idx in self.indices.values())

    @property
    def missing_indices(self) -> list[str]:
        """Retorna lista de índices ausentes."""
        return [
            idx.name for idx in self.indices.values() if not idx.exists
        ]


# ============================================================
# REFERENCE MANAGER
# ============================================================

# Índices obrigatórios e suas extensões
REQUIRED_INDICES = {
    "bwa_amb": ".amb",
    "bwa_ann": ".ann",
    "bwa_bwt": ".bwt",
    "bwa_pac": ".pac",
    "bwa_sa": ".sa",
    "samtools_fai": ".fai",
}

# O dict do GATK tem lógica especial (troca extensão em vez de anexar)
GATK_DICT_EXT = ".dict"

# Cache de metadados
CACHE_FILENAME = ".reference_meta.json"


class ReferenceManager:
    """Gerencia metadados e integridade de referências genômicas.

    Uso:
        manager = ReferenceManager("references/hg38.fa", config)
        info = manager.collect_info()
        if not info.integrity_ok:
            print(f"Índices faltando: {info.missing_indices}")
    """

    def __init__(
        self,
        fasta_path: str,
        config: dict[str, Any] | None = None,
    ) -> None:
        """Inicializa o manager.

        Args:
            fasta_path: Caminho para o FASTA da referência.
            config: Configuração do pipeline (seção 'reference').
        """
        self.fasta_path = Path(fasta_path).resolve()
        self.config = config or {}
        self.ref_config = self.config.get("reference", {})
        self._cache_path = self.fasta_path.parent / CACHE_FILENAME

    def collect_info(self) -> ReferenceInfo:
        """Coleta todos os metadados da referência.

        Returns:
            ReferenceInfo com dados completos e status de integridade.
        """
        # Metadados do config ou auto-detectados
        name = self.ref_config.get("name", "") or self._detect_name()
        version = self.ref_config.get("version", "") or self._detect_version()
        organism = self.ref_config.get("organism", "Homo sapiens")
        assembly = self.ref_config.get("assembly", name) or name
        source = self.ref_config.get("source", "") or self._detect_source()

        # Tamanho e checksum
        fasta_size = 0
        checksum = ""
        if self.fasta_path.is_file():
            fasta_size = self.fasta_path.stat().st_size
            checksum = self._get_checksum()

        # Data de instalação
        install_date = ""
        if self.fasta_path.is_file():
            mtime = self.fasta_path.stat().st_mtime
            install_date = datetime.fromtimestamp(mtime).isoformat()

        # Índices
        indices = self._check_indices()

        # Integridade geral
        integrity_ok = (
            self.fasta_path.is_file()
            and fasta_size > 0
            and all(idx.exists for idx in indices.values())
        )

        return ReferenceInfo(
            name=name,
            version=version,
            organism=organism,
            assembly=assembly,
            source=source,
            fasta_path=str(self.fasta_path),
            fasta_size_bytes=fasta_size,
            checksum_sha256=checksum,
            install_date=install_date,
            indices=indices,
            integrity_ok=integrity_ok,
        )

    def validate(self) -> tuple[bool, list[str]]:
        """Valida a referência completamente.

        Returns:
            Tupla (válida, lista de erros/avisos).
        """
        errors: list[str] = []
        info = self.collect_info()

        if not self.fasta_path.is_file():
            errors.append(
                f"FAIL: Arquivo FASTA não encontrado: {self.fasta_path}\n"
                f"  → Ação: Baixe o genoma de referência e coloque em: {self.fasta_path}"
            )
            return False, errors

        if info.fasta_size_bytes == 0:
            errors.append(
                f"FAIL: Arquivo FASTA vazio: {self.fasta_path}\n"
                f"  → Ação: Verifique se o download foi concluído corretamente"
            )

        for idx_name, idx_status in info.indices.items():
            if not idx_status.exists:
                tool = self._index_to_tool(idx_name)
                errors.append(
                    f"FAIL: Índice ausente: {idx_status.path} ({idx_name})\n"
                    f"  → Ação: Execute a indexação com {tool} ou use "
                    f"'snakemake --until {self._index_to_rule(idx_name)}'"
                )

        return len(errors) == 0, errors

    # ---- Detecção automática ----

    def _detect_name(self) -> str:
        """Detecta o nome da referência pelo nome do arquivo."""
        stem = self.fasta_path.stem
        if stem.endswith(".fa"):
            stem = stem[:-3]

        # Padrões comuns
        patterns = {
            r"(?i)grch38|hg38": "GRCh38",
            r"(?i)grch37|hg19": "GRCh37",
            r"(?i)t2t|chm13": "T2T-CHM13",
            r"(?i)mm10|grcm38": "GRCm38",
            r"(?i)mm39|grcm39": "GRCm39",
        }
        for pattern, name in patterns.items():
            if re.search(pattern, stem):
                return name

        return stem

    def _detect_version(self) -> str:
        """Detecta a versão da referência pelo nome do arquivo."""
        stem = self.fasta_path.stem
        # Procura padrões como p14, release-110, v2.0
        match = re.search(r"[._-](p\d+|release[._-]?\d+|v\d+\.\d+)", stem, re.I)
        return match.group(1) if match else ""

    def _detect_source(self) -> str:
        """Detecta a fonte da referência pelo caminho/nome."""
        full = str(self.fasta_path).lower()
        if "ensembl" in full:
            return "Ensembl"
        if "ncbi" in full or "gencode" in full:
            return "NCBI/GENCODE"
        if "ucsc" in full:
            return "UCSC"
        return ""

    # ---- Checksums ----

    def _get_checksum(self) -> str:
        """Obtém SHA256, usando cache quando possível.

        O checksum é cacheado em .reference_meta.json. Se o tamanho
        ou timestamp do arquivo mudou, recalcula.
        """
        cached = self._load_cache()
        if cached:
            cached_size = cached.get("fasta_size_bytes", 0)
            cached_mtime = cached.get("fasta_mtime", 0)
            current_stat = self.fasta_path.stat()

            if (
                cached_size == current_stat.st_size
                and abs(cached_mtime - current_stat.st_mtime) < 1.0
                and cached.get("checksum_sha256")
            ):
                return cached["checksum_sha256"]

        # Precisa recalcular
        checksum = self._compute_sha256()
        self._save_cache(checksum)
        return checksum

    def _compute_sha256(self, chunk_size: int = 65536) -> str:
        """Calcula SHA256 do arquivo FASTA."""
        sha256 = hashlib.sha256()
        with open(self.fasta_path, "rb") as f:
            while chunk := f.read(chunk_size):
                sha256.update(chunk)
        return sha256.hexdigest()

    def _load_cache(self) -> dict[str, Any] | None:
        """Carrega cache de metadados."""
        if self._cache_path.is_file():
            try:
                with open(self._cache_path) as f:
                    return json.load(f)
            except (json.JSONDecodeError, OSError):
                return None
        return None

    def _save_cache(self, checksum: str) -> None:
        """Salva cache de metadados."""
        stat = self.fasta_path.stat()
        cache = {
            "fasta_path": str(self.fasta_path),
            "fasta_size_bytes": stat.st_size,
            "fasta_mtime": stat.st_mtime,
            "checksum_sha256": checksum,
            "cached_at": datetime.now().isoformat(),
        }
        try:
            with open(self._cache_path, "w") as f:
                json.dump(cache, f, indent=2)
        except OSError:
            pass  # Cache é opcional — não falha se não puder salvar

    # ---- Índices ----

    def _check_indices(self) -> dict[str, IndexStatus]:
        """Verifica existência de todos os índices."""
        indices: dict[str, IndexStatus] = {}

        # Índices que adicionam extensão ao FASTA
        for idx_name, ext in REQUIRED_INDICES.items():
            idx_path = Path(f"{self.fasta_path}{ext}")
            exists = idx_path.is_file()
            size = idx_path.stat().st_size if exists else 0
            indices[idx_name] = IndexStatus(
                name=idx_name,
                extension=ext,
                path=str(idx_path),
                exists=exists,
                size_bytes=size,
            )

        # GATK dict (troca extensão em vez de anexar)
        dict_path = self.fasta_path.with_suffix(GATK_DICT_EXT)
        exists = dict_path.is_file()
        size = dict_path.stat().st_size if exists else 0
        indices["gatk_dict"] = IndexStatus(
            name="gatk_dict",
            extension=GATK_DICT_EXT,
            path=str(dict_path),
            exists=exists,
            size_bytes=size,
        )

        return indices

    @staticmethod
    def _index_to_tool(idx_name: str) -> str:
        """Mapeia nome do índice para a ferramenta que o cria."""
        mapping = {
            "bwa_amb": "bwa index",
            "bwa_ann": "bwa index",
            "bwa_bwt": "bwa index",
            "bwa_pac": "bwa index",
            "bwa_sa": "bwa index",
            "samtools_fai": "samtools faidx",
            "gatk_dict": "gatk CreateSequenceDictionary",
        }
        return mapping.get(idx_name, "ferramenta desconhecida")

    @staticmethod
    def _index_to_rule(idx_name: str) -> str:
        """Mapeia nome do índice para a regra Snakemake."""
        mapping = {
            "bwa_amb": "bwa_index",
            "bwa_ann": "bwa_index",
            "bwa_bwt": "bwa_index",
            "bwa_pac": "bwa_index",
            "bwa_sa": "bwa_index",
            "samtools_fai": "samtools_faidx",
            "gatk_dict": "gatk_dict",
        }
        return mapping.get(idx_name, idx_name)
