# ============================================================
# Testes — Detecção de Amostras
# ============================================================
"""Testes para o módulo de detecção automática de amostras."""

import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

from common.utils import detect_samples


# Padrões padrão
R1_PATTERNS = [
    "_R1_001.fastq.gz", "_R1_001.fastq",
    "_R1.fastq.gz", "_R1.fastq",
    "_1.fastq.gz", "_1.fastq",
]
R2_PATTERNS = [
    "_R2_001.fastq.gz", "_R2_001.fastq",
    "_R2.fastq.gz", "_R2.fastq",
    "_2.fastq.gz", "_2.fastq",
]


class TestDetectSamples:
    def test_paired_end_r1_r2(self, tmp_path: Path) -> None:
        """Detecta par R1/R2 padrão."""
        (tmp_path / "Sample01_R1.fastq.gz").write_bytes(b"")
        (tmp_path / "Sample01_R2.fastq.gz").write_bytes(b"")

        result = detect_samples(str(tmp_path), R1_PATTERNS, R2_PATTERNS)

        assert "Sample01" in result
        assert result["Sample01"]["type"] == "PE"
        assert "R1" in result["Sample01"]["r1"]
        assert "R2" in result["Sample01"]["r2"]

    def test_paired_end_illumina(self, tmp_path: Path) -> None:
        """Detecta par com padrão Illumina (_R1_001)."""
        (tmp_path / "SampleA_S1_L001_R1_001.fastq.gz").write_bytes(b"")
        (tmp_path / "SampleA_S1_L001_R2_001.fastq.gz").write_bytes(b"")

        result = detect_samples(str(tmp_path), R1_PATTERNS, R2_PATTERNS)

        assert "SampleA_S1_L001" in result
        assert result["SampleA_S1_L001"]["type"] == "PE"

    def test_single_end(self, tmp_path: Path) -> None:
        """Detecta amostra single-end."""
        (tmp_path / "Sample02.fastq.gz").write_bytes(b"")

        result = detect_samples(str(tmp_path), R1_PATTERNS, R2_PATTERNS)

        assert "Sample02" in result
        assert result["Sample02"]["type"] == "SE"
        assert result["Sample02"]["r2"] == ""

    def test_mixed_se_pe(self, tmp_path: Path) -> None:
        """Detecta mix de SE e PE."""
        (tmp_path / "PE_Sample_R1.fastq.gz").write_bytes(b"")
        (tmp_path / "PE_Sample_R2.fastq.gz").write_bytes(b"")
        (tmp_path / "SE_Sample.fastq.gz").write_bytes(b"")

        result = detect_samples(str(tmp_path), R1_PATTERNS, R2_PATTERNS)

        assert len(result) == 2
        assert result["PE_Sample"]["type"] == "PE"
        assert result["SE_Sample"]["type"] == "SE"

    def test_uncompressed_fastq(self, tmp_path: Path) -> None:
        """Detecta FASTQs não compactados."""
        (tmp_path / "Sample03_R1.fastq").write_bytes(b"")
        (tmp_path / "Sample03_R2.fastq").write_bytes(b"")

        result = detect_samples(str(tmp_path), R1_PATTERNS, R2_PATTERNS)

        assert "Sample03" in result
        assert result["Sample03"]["type"] == "PE"

    def test_empty_directory(self, tmp_path: Path) -> None:
        """Retorna vazio para diretório sem FASTQs."""
        result = detect_samples(str(tmp_path), R1_PATTERNS, R2_PATTERNS)
        assert len(result) == 0

    def test_nonexistent_directory(self) -> None:
        """Retorna vazio para diretório inexistente."""
        result = detect_samples("/nonexistent/dir", R1_PATTERNS, R2_PATTERNS)
        assert len(result) == 0

    def test_multiple_samples(self, tmp_path: Path) -> None:
        """Detecta múltiplas amostras corretamente."""
        for i in range(5):
            (tmp_path / f"Sample{i:02d}_R1.fastq.gz").write_bytes(b"")
            (tmp_path / f"Sample{i:02d}_R2.fastq.gz").write_bytes(b"")

        result = detect_samples(str(tmp_path), R1_PATTERNS, R2_PATTERNS)

        assert len(result) == 5
        for i in range(5):
            assert f"Sample{i:02d}" in result

    def test_ignores_non_fastq_files(self, tmp_path: Path) -> None:
        """Ignora arquivos que não são FASTQ."""
        (tmp_path / "Sample_R1.fastq.gz").write_bytes(b"")
        (tmp_path / "Sample_R2.fastq.gz").write_bytes(b"")
        (tmp_path / "README.md").write_text("test")
        (tmp_path / "data.csv").write_text("a,b,c")
        (tmp_path / "script.py").write_text("print('hello')")

        result = detect_samples(str(tmp_path), R1_PATTERNS, R2_PATTERNS)

        assert len(result) == 1
        assert "Sample" in result

    def test_underscore_numbered_pattern(self, tmp_path: Path) -> None:
        """Detecta padrão _1/_2."""
        (tmp_path / "SRR12345_1.fastq.gz").write_bytes(b"")
        (tmp_path / "SRR12345_2.fastq.gz").write_bytes(b"")

        result = detect_samples(str(tmp_path), R1_PATTERNS, R2_PATTERNS)

        assert "SRR12345" in result
        assert result["SRR12345"]["type"] == "PE"
