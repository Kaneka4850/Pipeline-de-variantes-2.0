# ============================================================
# Pipeline de Variantes Germinativas — Utilitários
# ============================================================
"""
Funções utilitárias compartilhadas pelo Snakefile e scripts.

Responsabilidades:
- Resolução de caminhos de ferramentas
- Detecção automática de amostras (SE/PE)
- Detecção de versões de ferramentas
- Funções auxiliares de tempo e formatação
"""

import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Any


# ============================================================
# RESOLUÇÃO DE FERRAMENTAS
# ============================================================

def get_tool_path(config: dict[str, Any], tool_name: str) -> str:
    """Retorna o caminho do executável de uma ferramenta.

    Se o config define um caminho explícito, usa esse caminho.
    Caso contrário, tenta encontrar no PATH do sistema.

    Args:
        config: Dicionário de configuração do pipeline.
        tool_name: Nome da ferramenta (bwa, gatk, samtools, etc.).

    Returns:
        Caminho do executável ou o nome da ferramenta (fallback para PATH).
    """
    configured_path = config.get("tools", {}).get(tool_name, "")
    if configured_path:
        return configured_path

    # Tenta localizar no PATH
    found = shutil.which(tool_name)
    if found:
        return found

    # Fallback: retorna o nome da ferramenta (Snakemake falhará com msg clara)
    return tool_name


def get_tool_version(tool_path: str) -> str:
    """Detecta a versão de uma ferramenta de bioinformática.

    Tenta múltiplas estratégias, pois cada ferramenta reporta versão
    de forma diferente (--version, version, sem argumentos, etc.).

    Args:
        tool_path: Caminho ou nome do executável.

    Returns:
        String com a versão, ou "desconhecida" se não conseguir detectar.
    """
    strategies = [
        [tool_path, "--version"],
        [tool_path, "version"],
        [tool_path, "-v"],
    ]

    for cmd in strategies:
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=15,
            )
            output = result.stdout.strip() or result.stderr.strip()
            if output:
                # Extrai primeira linha relevante
                for line in output.splitlines():
                    line = line.strip()
                    if line and not line.startswith("Usage"):
                        return line
        except (subprocess.SubprocessError, FileNotFoundError, OSError):
            continue

    return "desconhecida"


# ============================================================
# DETECÇÃO AUTOMÁTICA DE AMOSTRAS
# ============================================================

def detect_samples(
    samples_dir: str,
    r1_patterns: list[str],
    r2_patterns: list[str],
) -> dict[str, dict[str, str]]:
    """Detecta amostras automaticamente no diretório de FASTQs.

    Escaneia o diretório procurando por padrões R1/R2 e agrupa
    arquivos em paired-end ou single-end.

    Args:
        samples_dir: Diretório onde estão os FASTQs.
        r1_patterns: Lista de sufixos R1 (ex: ["_R1.fastq.gz"]).
        r2_patterns: Lista de sufixos R2 (ex: ["_R2.fastq.gz"]).

    Returns:
        Dicionário {sample_name: {"r1": path, "r2": path|"", "type": "PE"|"SE"}}.

    Raises:
        FileNotFoundError: Se o diretório de amostras não existir.
    """
    samples_path = Path(samples_dir)
    if not samples_path.is_dir():
        return {}

    # Coletar todos os arquivos FASTQ
    all_files = []
    for f in samples_path.iterdir():
        if f.is_file():
            # Verifica se tem extensão de FASTQ (considerando .gz duplo)
            name = f.name
            if any(name.endswith(ext) for ext in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]):
                all_files.append(f)

    if not all_files:
        return {}

    samples: dict[str, dict[str, str]] = {}
    matched_files: set[str] = set()

    # Primeiro, detecta paired-end (procura pares R1/R2)
    for r1_pattern in r1_patterns:
        for f in all_files:
            if f.name in matched_files:
                continue
            if f.name.endswith(r1_pattern):
                sample_name = f.name[: -len(r1_pattern)]

                # Procura o R2 correspondente
                for r2_pattern in r2_patterns:
                    r2_file = samples_path / f"{sample_name}{r2_pattern}"
                    if r2_file.is_file():
                        samples[sample_name] = {
                            "r1": str(f),
                            "r2": str(r2_file),
                            "type": "PE",
                        }
                        matched_files.add(f.name)
                        matched_files.add(r2_file.name)
                        break

    # Depois, detecta single-end (arquivos sem par)
    for f in all_files:
        if f.name in matched_files:
            continue

        # Remove extensão para obter nome da amostra
        sample_name = f.name
        for ext in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
            if sample_name.endswith(ext):
                sample_name = sample_name[: -len(ext)]
                break

        # Verifica se não é um R2 órfão
        is_r2 = any(sample_name.endswith(p.replace(".fastq.gz", "").replace(".fastq", "").replace(".fq.gz", "").replace(".fq", "")) for p in r2_patterns)
        if not is_r2:
            samples[sample_name] = {
                "r1": str(f),
                "r2": "",
                "type": "SE",
            }
            matched_files.add(f.name)

    return dict(sorted(samples.items()))


# ============================================================
# FUNÇÕES DE TEMPO
# ============================================================

def format_duration(seconds: float) -> str:
    """Formata uma duração em segundos para formato legível.

    Args:
        seconds: Duração em segundos.

    Returns:
        String formatada (ex: "1h 23m 45s", "5m 12s", "3s").
    """
    hours, remainder = divmod(int(seconds), 3600)
    minutes, secs = divmod(remainder, 60)

    if hours > 0:
        return f"{hours}h {minutes:02d}m {secs:02d}s"
    elif minutes > 0:
        return f"{minutes}m {secs:02d}s"
    else:
        return f"{secs}s"


def get_timestamp() -> str:
    """Retorna timestamp atual no formato ISO 8601."""
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


# ============================================================
# FUNÇÕES DE SISTEMA
# ============================================================

def get_hostname() -> str:
    """Retorna o hostname da máquina."""
    import socket
    return socket.gethostname()


def get_username() -> str:
    """Retorna o nome do usuário atual."""
    import getpass
    return getpass.getuser()


def ensure_dir(path: str | Path) -> Path:
    """Cria diretório se não existir. Retorna o Path.

    Args:
        path: Caminho do diretório.

    Returns:
        Path do diretório criado/existente.
    """
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p
