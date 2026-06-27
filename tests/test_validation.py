# ============================================================
# Testes — Validações
# ============================================================
"""Testes para o módulo de validação do pipeline."""

import sys
import gzip
from pathlib import Path

import pytest

# Adiciona scripts ao path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

from common.validation import (
    validate_fastq,
    validate_reference,
    check_disk_space,
    check_write_permission,
    check_tool_installed,
    estimate_required_space,
)


# ============================================================
# FIXTURES
# ============================================================

@pytest.fixture
def tmp_dir(tmp_path: Path) -> Path:
    """Diretório temporário para testes."""
    return tmp_path


@pytest.fixture
def valid_fastq(tmp_dir: Path) -> Path:
    """Cria um FASTQ válido não compactado."""
    fq = tmp_dir / "test.fastq"
    fq.write_text(
        "@SEQ_ID_1\n"
        "ACGTACGTACGT\n"
        "+\n"
        "IIIIIIIIIIII\n"
    )
    return fq


@pytest.fixture
def valid_fastq_gz(tmp_dir: Path) -> Path:
    """Cria um FASTQ válido compactado."""
    fq = tmp_dir / "test.fastq.gz"
    content = b"@SEQ_ID_1\nACGTACGT\n+\nIIIIIIII\n"
    with gzip.open(fq, "wb") as f:
        f.write(content)
    return fq


@pytest.fixture
def valid_reference(tmp_dir: Path) -> Path:
    """Cria um FASTA de referência válido."""
    ref = tmp_dir / "ref.fa"
    ref.write_text(">chr1\nACGTACGT\n")
    return ref


# ============================================================
# TESTES: validate_fastq
# ============================================================

class TestValidateFastq:
    def test_valid_fastq(self, valid_fastq: Path) -> None:
        valid, msg = validate_fastq(str(valid_fastq))
        assert valid is True
        assert msg == "OK"

    def test_valid_fastq_gz(self, valid_fastq_gz: Path) -> None:
        valid, msg = validate_fastq(str(valid_fastq_gz))
        assert valid is True
        assert msg == "OK"

    def test_nonexistent_file(self) -> None:
        valid, msg = validate_fastq("/nonexistent/file.fastq")
        assert valid is False
        assert "não encontrado" in msg

    def test_empty_file(self, tmp_dir: Path) -> None:
        empty = tmp_dir / "empty.fastq"
        empty.write_text("")
        valid, msg = validate_fastq(str(empty))
        assert valid is False
        assert "vazio" in msg.lower()

    def test_invalid_fastq_format(self, tmp_dir: Path) -> None:
        bad = tmp_dir / "bad.fastq"
        bad.write_text("NOT_A_FASTQ\nACGT\n+\nIIII\n")
        valid, msg = validate_fastq(str(bad))
        assert valid is False
        assert "não começa com '@'" in msg

    def test_corrupted_gz(self, tmp_dir: Path) -> None:
        bad_gz = tmp_dir / "bad.fastq.gz"
        bad_gz.write_bytes(b"NOT_GZIP_DATA")
        valid, msg = validate_fastq(str(bad_gz))
        assert valid is False
        assert "magic number" in msg.lower() or "corrompido" in msg.lower()


# ============================================================
# TESTES: validate_reference
# ============================================================

class TestValidateReference:
    def test_valid_reference(self, valid_reference: Path) -> None:
        valid, errors = validate_reference(str(valid_reference))
        assert valid is True
        assert len(errors) == 0

    def test_nonexistent_reference(self) -> None:
        valid, errors = validate_reference("/nonexistent/ref.fa")
        assert valid is False
        assert len(errors) > 0

    def test_empty_reference(self, tmp_dir: Path) -> None:
        empty = tmp_dir / "empty.fa"
        empty.write_text("")
        valid, errors = validate_reference(str(empty))
        assert valid is False


# ============================================================
# TESTES: check_disk_space
# ============================================================

class TestCheckDiskSpace:
    def test_sufficient_space(self, tmp_dir: Path) -> None:
        ok, msg = check_disk_space(str(tmp_dir), required_gb=0.001)
        assert ok is True

    def test_insufficient_space(self, tmp_dir: Path) -> None:
        ok, msg = check_disk_space(str(tmp_dir), required_gb=999999)
        assert ok is False
        assert "insuficiente" in msg.lower()


# ============================================================
# TESTES: check_write_permission
# ============================================================

class TestCheckWritePermission:
    def test_writable_dir(self, tmp_dir: Path) -> None:
        ok, msg = check_write_permission(str(tmp_dir))
        assert ok is True

    def test_create_new_dir(self, tmp_dir: Path) -> None:
        new_dir = tmp_dir / "new_subdir"
        ok, msg = check_write_permission(str(new_dir))
        assert ok is True
        assert new_dir.is_dir()


# ============================================================
# TESTES: check_tool_installed
# ============================================================

class TestCheckToolInstalled:
    def test_python_installed(self) -> None:
        ok, path = check_tool_installed("python3")
        # Python3 pode ou não estar disponível no Windows
        # Tenta python como fallback
        if not ok:
            ok, path = check_tool_installed("python")
        assert ok is True

    def test_nonexistent_tool(self) -> None:
        ok, msg = check_tool_installed("ferramenta_inexistente_xyz")
        assert ok is False


# ============================================================
# TESTES: estimate_required_space
# ============================================================

class TestEstimateRequiredSpace:
    def test_estimate_with_files(self, valid_fastq: Path) -> None:
        estimate = estimate_required_space([str(valid_fastq)])
        assert estimate >= 5.0  # Mínimo de 5 GB

    def test_estimate_empty_list(self) -> None:
        estimate = estimate_required_space([])
        assert estimate >= 5.0

    def test_estimate_nonexistent_files(self) -> None:
        estimate = estimate_required_space(["/nonexistent/file.fastq.gz"])
        assert estimate >= 5.0
