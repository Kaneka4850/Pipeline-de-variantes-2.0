# ============================================================
# Pipeline de Variantes Germinativas — Detecção de Amostras
# ============================================================
"""
Script para detecção automática de amostras e geração do
arquivo samples.tsv.

Pode ser executado diretamente:
    python scripts/parse_samples.py --samples-dir samples/ --output config/samples.tsv
"""

import argparse
import csv
import sys
from pathlib import Path

# Adiciona raiz do projeto ao path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from common.utils import detect_samples


def write_samples_tsv(
    samples: dict[str, dict[str, str]],
    output_path: str,
) -> None:
    """Escreve o arquivo samples.tsv.

    Args:
        samples: Dicionário de amostras detectadas.
        output_path: Caminho para o arquivo de saída.
    """
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sample", "r1", "r2", "type"])
        for name, info in samples.items():
            writer.writerow([
                name,
                info["r1"],
                info.get("r2", ""),
                info["type"],
            ])


def main() -> None:
    """Ponto de entrada para execução via linha de comando."""
    parser = argparse.ArgumentParser(
        description="Detecta amostras FASTQ automaticamente e gera samples.tsv"
    )
    parser.add_argument(
        "--samples-dir",
        default="samples",
        help="Diretório com arquivos FASTQ (default: samples/)",
    )
    parser.add_argument(
        "--output",
        default="config/samples.tsv",
        help="Caminho para o arquivo de saída (default: config/samples.tsv)",
    )
    args = parser.parse_args()

    # Padrões padrão
    r1_patterns = [
        "_R1_001.fastq.gz", "_R1_001.fastq",
        "_R1.fastq.gz", "_R1.fastq",
        "_1.fastq.gz", "_1.fastq",
    ]
    r2_patterns = [
        "_R2_001.fastq.gz", "_R2_001.fastq",
        "_R2.fastq.gz", "_R2.fastq",
        "_2.fastq.gz", "_2.fastq",
    ]

    samples = detect_samples(args.samples_dir, r1_patterns, r2_patterns)

    if not samples:
        print(f"[ERRO] Nenhuma amostra encontrada em '{args.samples_dir}/'")
        sys.exit(1)

    write_samples_tsv(samples, args.output)

    print(f"[INFO] {len(samples)} amostra(s) detectada(s):")
    for name, info in samples.items():
        print(f"       → {name} ({info['type']})")
    print(f"[INFO] Arquivo gerado: {args.output}")


if __name__ == "__main__":
    main()
