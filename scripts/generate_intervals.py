# ============================================================
# Pipeline de Variantes Germinativas — Geração de Intervalos
# ============================================================
"""
Script para gerar arquivos BED com intervalos genômicos 
otimizados para paralelização (WGS), a partir do índice .fai.
"""

import argparse
from pathlib import Path


def generate_windows(fai_path: str, bed_path: str, window_size: int = 50000000) -> None:
    """Gera janelas BED a partir do FAI.
    
    Apenas considera os cromossomos principais (1-22, X, Y, MT) para
    evitar sobrecarga com contigs alternativos e decoys pequenos.
    """
    valid_chroms = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM", "chrMT"}
    valid_chroms.update({str(i) for i in range(1, 23)} | {"X", "Y", "M", "MT"})

    Path(bed_path).parent.mkdir(parents=True, exist_ok=True)

    with open(fai_path, 'r') as f_in, open(bed_path, 'w') as f_out:
        for line in f_in:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            
            chrom = parts[0]
            if chrom not in valid_chroms:
                continue
                
            length = int(parts[1])
            start = 0
            while start < length:
                end = min(start + window_size, length)
                f_out.write(f"{chrom}\t{start}\t{end}\n")
                start = end


def main() -> None:
    parser = argparse.ArgumentParser(description="Gera intervalos BED a partir de um arquivo FAI")
    parser.add_argument("fai", help="Caminho para o arquivo .fai")
    parser.add_argument("bed", help="Caminho para o arquivo de saída .bed")
    parser.add_argument("--window", type=int, default=50000000, help="Tamanho da janela em pares de bases (default: 50MB)")
    
    args = parser.parse_args()
    generate_windows(args.fai, args.bed, args.window)


if __name__ == "__main__":
    main()
