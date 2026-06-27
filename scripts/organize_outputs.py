# ============================================================
# Pipeline de Variantes Germinativas — Organização de Outputs
# ============================================================
"""
Script para organizar os outputs finais por amostra na
estrutura padronizada.

Estrutura gerada:
    results/{sample}/
        ├── fastq/        # Links simbólicos para FASTQs originais
        ├── qc/           # Relatórios FastQC
        ├── alignment/    # BAMs, índices, métricas
        ├── variants/     # VCFs, GVCFs
        ├── annotation/   # VCF anotado, TSV, HTML
        ├── reports/      # Relatórios finais
        └── logs/         # Logs da amostra
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from common.utils import ensure_dir


SAMPLE_SUBDIRS = [
    "fastq",
    "qc",
    "alignment",
    "variants",
    "annotation",
    "reports",
    "logs",
]


def create_sample_structure(results_dir: str, sample_name: str) -> dict[str, Path]:
    """Cria a estrutura de diretórios para uma amostra.

    Args:
        results_dir: Diretório raiz de resultados.
        sample_name: Nome da amostra.

    Returns:
        Dicionário com caminhos criados.
    """
    paths: dict[str, Path] = {}
    for subdir in SAMPLE_SUBDIRS:
        path = ensure_dir(Path(results_dir) / sample_name / subdir)
        paths[subdir] = path
    return paths


def main() -> None:
    """Cria estrutura de diretórios para amostras listadas."""
    parser = argparse.ArgumentParser(
        description="Cria estrutura de diretórios para outputs do pipeline"
    )
    parser.add_argument(
        "--results-dir",
        default="results",
        help="Diretório raiz de resultados (default: results/)",
    )
    parser.add_argument(
        "samples",
        nargs="+",
        help="Nomes das amostras",
    )
    args = parser.parse_args()

    for sample in args.samples:
        paths = create_sample_structure(args.results_dir, sample)
        print(f"[INFO] Estrutura criada para: {sample}")
        for name, path in paths.items():
            print(f"       → {path}")


if __name__ == "__main__":
    main()
