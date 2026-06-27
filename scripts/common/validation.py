# ============================================================
# Pipeline de Variantes Germinativas — Validações
# ============================================================
"""
Validações robustas para o pipeline.

Verifica:
- Existência e integridade de arquivos FASTQ
- Integridade da referência genômica (FASTA + índices)
- Espaço em disco disponível
- Permissões de escrita
- Ferramentas instaladas no PATH
- Parâmetros obrigatórios
"""

import gzip
import os
import shutil
from pathlib import Path
from typing import Any


# ============================================================
# VALIDAÇÃO DE FASTQ
# ============================================================

def validate_fastq(file_path: str) -> tuple[bool, str]:
    """Valida integridade básica de um arquivo FASTQ.

    Verifica:
    - Arquivo existe
    - Arquivo não está vazio
    - Se .gz, verifica magic number do gzip
    - Primeira linha começa com '@' (cabeçalho FASTQ)

    Args:
        file_path: Caminho para o arquivo FASTQ.

    Returns:
        Tupla (válido, mensagem).
    """
    path = Path(file_path)

    if not path.exists():
        return False, f"Arquivo não encontrado: {file_path}"

    if not path.is_file():
        return False, f"Não é um arquivo regular: {file_path}"

    if path.stat().st_size == 0:
        return False, f"Arquivo vazio: {file_path}"

    # Verifica magic number do gzip
    if path.name.endswith(".gz"):
        try:
            with open(file_path, "rb") as f:
                magic = f.read(2)
                if magic != b"\x1f\x8b":
                    return False, f"Arquivo .gz corrompido (magic number inválido): {file_path}"
        except OSError as e:
            return False, f"Erro ao ler arquivo: {file_path}: {e}"

        # Verifica primeira linha do conteúdo descompactado
        try:
            with gzip.open(file_path, "rt") as f:
                first_line = f.readline().strip()
                if not first_line.startswith("@"):
                    return False, (
                        f"Formato FASTQ inválido (primeira linha não começa com '@'): "
                        f"{file_path}"
                    )
        except (gzip.BadGzipFile, OSError) as e:
            return False, f"Arquivo gzip corrompido: {file_path}: {e}"
    else:
        # FASTQ não compactado
        try:
            with open(file_path, "r") as f:
                first_line = f.readline().strip()
                if not first_line.startswith("@"):
                    return False, (
                        f"Formato FASTQ inválido (primeira linha não começa com '@'): "
                        f"{file_path}"
                    )
        except OSError as e:
            return False, f"Erro ao ler arquivo: {file_path}: {e}"

    return True, "OK"


# ============================================================
# VALIDAÇÃO DE REFERÊNCIA
# ============================================================

def validate_reference(genome_path: str) -> tuple[bool, list[str]]:
    """Valida a existência do genoma de referência.

    Verifica apenas se o FASTA principal existe.
    Os índices são gerados automaticamente pelo pipeline.

    Args:
        genome_path: Caminho para o FASTA do genoma.

    Returns:
        Tupla (válido, lista de mensagens de erro).
    """
    errors: list[str] = []
    path = Path(genome_path)

    if not path.exists():
        errors.append(f"Genoma de referência não encontrado: {genome_path}")
    elif path.stat().st_size == 0:
        errors.append(f"Genoma de referência vazio: {genome_path}")

    return len(errors) == 0, errors


def check_reference_indices(genome_path: str) -> dict[str, bool]:
    """Verifica quais índices da referência existem.

    Args:
        genome_path: Caminho para o FASTA do genoma.

    Returns:
        Dicionário com status de cada índice.
    """
    extensions = {
        "bwa_amb": ".amb",
        "bwa_ann": ".ann",
        "bwa_bwt": ".bwt",
        "bwa_pac": ".pac",
        "bwa_sa": ".sa",
        "samtools_fai": ".fai",
        "gatk_dict": ".dict",
    }

    # Para .dict, o nome é diferente (troca extensão do FASTA)
    base = Path(genome_path)
    results: dict[str, bool] = {}
    for name, ext in extensions.items():
        if name == "gatk_dict":
            dict_path = base.with_suffix(".dict")
            results[name] = dict_path.exists()
        else:
            results[name] = Path(f"{genome_path}{ext}").exists()

    return results


# ============================================================
# VALIDAÇÃO DE DISCO
# ============================================================

def check_disk_space(
    directory: str,
    required_gb: float = 10.0,
) -> tuple[bool, str]:
    """Verifica se há espaço em disco suficiente.

    Args:
        directory: Diretório para verificar.
        required_gb: Espaço mínimo requerido em GB.

    Returns:
        Tupla (suficiente, mensagem).
    """
    try:
        stat = shutil.disk_usage(directory)
        free_gb = stat.free / (1024 ** 3)
        if free_gb < required_gb:
            return False, (
                f"Espaço em disco insuficiente em '{directory}': "
                f"{free_gb:.1f} GB livre, {required_gb:.1f} GB necessário"
            )
        return True, f"Espaço disponível: {free_gb:.1f} GB"
    except OSError as e:
        return False, f"Erro ao verificar disco: {e}"


def estimate_required_space(fastq_paths: list[str]) -> float:
    """Estima o espaço necessário baseado nos FASTQs.

    Heurística: pipeline gera ~10x o tamanho dos FASTQs
    (BAMs intermediários, VCFs, etc.)

    Args:
        fastq_paths: Lista de caminhos para os FASTQs.

    Returns:
        Espaço estimado em GB.
    """
    total_bytes = sum(
        os.path.getsize(p)
        for p in fastq_paths
        if os.path.isfile(p)
    )
    # Fator de expansão ~10x para todo o pipeline
    estimated = (total_bytes * 10) / (1024 ** 3)
    # Mínimo de 5 GB
    return max(estimated, 5.0)


# ============================================================
# VALIDAÇÃO DE PERMISSÕES
# ============================================================

def check_write_permission(directory: str) -> tuple[bool, str]:
    """Verifica permissão de escrita em um diretório.

    Args:
        directory: Caminho do diretório.

    Returns:
        Tupla (pode_escrever, mensagem).
    """
    path = Path(directory)

    if not path.exists():
        # Tenta criar
        try:
            path.mkdir(parents=True, exist_ok=True)
            return True, f"Diretório criado: {directory}"
        except OSError as e:
            return False, f"Não foi possível criar diretório '{directory}': {e}"

    if not os.access(directory, os.W_OK):
        return False, f"Sem permissão de escrita em: {directory}"

    return True, "OK"


# ============================================================
# VALIDAÇÃO DE FERRAMENTAS
# ============================================================

def check_tool_installed(tool_name: str, tool_path: str = "") -> tuple[bool, str]:
    """Verifica se uma ferramenta está instalada e acessível.

    Args:
        tool_name: Nome da ferramenta.
        tool_path: Caminho explícito (se configurado).

    Returns:
        Tupla (instalada, caminho_ou_erro).
    """
    path_to_check = tool_path if tool_path else tool_name

    found = shutil.which(path_to_check)
    if found:
        return True, found

    if tool_path and os.path.isfile(tool_path) and os.access(tool_path, os.X_OK):
        return True, tool_path

    return False, f"Ferramenta '{tool_name}' não encontrada no PATH"


# ============================================================
# VALIDAÇÃO COMPLETA DO PIPELINE
# ============================================================

def validate_pipeline_config(config: dict[str, Any]) -> tuple[bool, list[str]]:
    """Executa todas as validações do pipeline.

    Args:
        config: Dicionário de configuração do pipeline.

    Returns:
        Tupla (válido, lista de mensagens de erro).
    """
    errors: list[str] = []
    warnings: list[str] = []

    # 1. Referência
    genome = config.get("reference", {}).get("genome", "")
    if not genome:
        errors.append("Genoma de referência não configurado (reference.genome)")
    else:
        valid, ref_errors = validate_reference(genome)
        if not valid:
            errors.extend(ref_errors)

    # 2. Diretório de amostras
    samples_dir = config.get("dirs", {}).get("samples", "")
    if not samples_dir:
        errors.append("Diretório de amostras não configurado (dirs.samples)")
    elif not Path(samples_dir).is_dir():
        errors.append(f"Diretório de amostras não existe: {samples_dir}")

    # 3. Permissões de escrita nos diretórios de saída
    for dir_key in ["results", "logs", "database"]:
        dir_path = config.get("dirs", {}).get(dir_key, "")
        if dir_path:
            can_write, msg = check_write_permission(dir_path)
            if not can_write:
                errors.append(msg)

    # 4. Ferramentas obrigatórias
    required_tools = ["bwa", "gatk", "samtools", "fastqc"]
    for tool in required_tools:
        tool_path = config.get("tools", {}).get(tool, "")
        installed, msg = check_tool_installed(tool, tool_path)
        if not installed:
            errors.append(msg)

    # 5. VEP (se anotação habilitada)
    if config.get("pipeline", {}).get("annotation", True):
        vep_path = config.get("tools", {}).get("vep", "")
        installed, msg = check_tool_installed("vep", vep_path)
        if not installed:
            warnings.append(f"VEP não encontrado — anotação será ignorada: {msg}")

    # 6. Recursos
    threads = config.get("resources", {}).get("threads", 0)
    if threads < 1:
        errors.append("Número de threads deve ser >= 1 (resources.threads)")

    memory = config.get("resources", {}).get("memory_mb", 0)
    if memory < 1000:
        errors.append("Memória deve ser >= 1000 MB (resources.memory_mb)")

    return len(errors) == 0, errors + warnings
