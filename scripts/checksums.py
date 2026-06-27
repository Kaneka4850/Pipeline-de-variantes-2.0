# ============================================================
# Pipeline de Variantes Germinativas — Checksums
# ============================================================
"""
Geração de checksums MD5 para arquivos críticos do pipeline.

Uso:
    python scripts/checksums.py results/Sample01/
"""

import argparse
import hashlib
import sys
from pathlib import Path


CHECKSUM_EXTENSIONS = {
    ".bam", ".bam.bai", ".bai",
    ".vcf.gz", ".vcf.gz.tbi",
    ".g.vcf.gz", ".g.vcf.gz.tbi",
    ".fastq.gz", ".fq.gz",
    ".tsv",
}


def compute_md5(file_path: Path, chunk_size: int = 8192) -> str:
    """Calcula o hash MD5 de um arquivo.

    Args:
        file_path: Caminho do arquivo.
        chunk_size: Tamanho do bloco de leitura.

    Returns:
        String hexadecimal do hash MD5.
    """
    md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        while chunk := f.read(chunk_size):
            md5.update(chunk)
    return md5.hexdigest()


def generate_checksums(directory: Path) -> dict[str, str]:
    """Gera checksums para todos os arquivos críticos em um diretório.

    Args:
        directory: Diretório para escanear recursivamente.

    Returns:
        Dicionário {caminho_relativo: checksum_md5}.
    """
    checksums: dict[str, str] = {}

    for file_path in sorted(directory.rglob("*")):
        if not file_path.is_file():
            continue

        # Verifica se a extensão é relevante
        name = file_path.name
        is_critical = any(name.endswith(ext) for ext in CHECKSUM_EXTENSIONS)
        if is_critical:
            rel_path = str(file_path.relative_to(directory))
            checksums[rel_path] = compute_md5(file_path)

    return checksums


def write_checksums_file(
    checksums: dict[str, str],
    output_path: Path,
) -> None:
    """Escreve arquivo de checksums no formato MD5SUM.

    Args:
        checksums: Dicionário {caminho: hash}.
        output_path: Caminho do arquivo de saída.
    """
    with open(output_path, "w") as f:
        for path, md5 in checksums.items():
            f.write(f"{md5}  {path}\n")


def main() -> None:
    """Ponto de entrada para execução via CLI."""
    parser = argparse.ArgumentParser(
        description="Gera checksums MD5 para arquivos críticos do pipeline"
    )
    parser.add_argument(
        "directory",
        help="Diretório para escanear",
    )
    parser.add_argument(
        "--output",
        default="",
        help="Arquivo de saída (default: {directory}/checksums.md5)",
    )
    args = parser.parse_args()

    directory = Path(args.directory)
    if not directory.is_dir():
        print(f"[ERRO] Diretório não encontrado: {directory}")
        sys.exit(1)

    output = Path(args.output) if args.output else directory / "checksums.md5"

    print(f"[INFO] Calculando checksums em: {directory}")
    checksums = generate_checksums(directory)

    if not checksums:
        print("[AVISO] Nenhum arquivo crítico encontrado.")
        return

    write_checksums_file(checksums, output)
    print(f"[INFO] {len(checksums)} checksum(s) gerado(s): {output}")
    for path, md5 in checksums.items():
        print(f"       {md5}  {path}")


if __name__ == "__main__":
    main()
