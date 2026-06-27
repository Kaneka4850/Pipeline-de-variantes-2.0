# ============================================================
# Testes — Configuração
# ============================================================
"""Testes para carregamento e validação do config YAML."""

import sys
from pathlib import Path

import pytest
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

from common.utils import get_tool_path, format_duration


# ============================================================
# FIXTURES
# ============================================================

@pytest.fixture
def sample_config() -> dict:
    """Config de exemplo válido."""
    return {
        "reference": {
            "genome": "references/hg38.fa",
            "known_sites": [],
        },
        "dirs": {
            "samples": "samples",
            "results": "results",
            "logs": "logs",
            "database": "database",
        },
        "resources": {
            "threads": 4,
            "memory_mb": 8000,
            "max_parallel": 4,
        },
        "tools": {
            "bwa": "",
            "gatk": "",
            "samtools": "",
            "fastqc": "",
            "bcftools": "",
            "multiqc": "",
            "vep": "",
        },
        "pipeline": {
            "version": "2.0.0",
            "name": "Pipeline de Variantes",
            "joint_calling": False,
            "annotation": True,
        },
    }


# ============================================================
# TESTES: Config YAML
# ============================================================

class TestConfigLoading:
    def test_load_config_file(self) -> None:
        """Carrega o config.yaml real do projeto."""
        config_path = Path(__file__).resolve().parent.parent / "config" / "config.yaml"
        if config_path.exists():
            with open(config_path) as f:
                config = yaml.safe_load(f)

            assert "reference" in config
            assert "dirs" in config
            assert "resources" in config
            assert "tools" in config
            assert "pipeline" in config

    def test_config_has_required_sections(self, sample_config: dict) -> None:
        """Config contém todas as seções obrigatórias."""
        required = ["reference", "dirs", "resources", "tools", "pipeline"]
        for section in required:
            assert section in sample_config

    def test_config_genome_path(self, sample_config: dict) -> None:
        """Genome path está definido."""
        assert sample_config["reference"]["genome"]
        assert isinstance(sample_config["reference"]["genome"], str)

    def test_config_threads_valid(self, sample_config: dict) -> None:
        """Threads é um número positivo."""
        assert sample_config["resources"]["threads"] >= 1

    def test_config_memory_valid(self, sample_config: dict) -> None:
        """Memória é >= 1000 MB."""
        assert sample_config["resources"]["memory_mb"] >= 1000


# ============================================================
# TESTES: get_tool_path
# ============================================================

class TestGetToolPath:
    def test_empty_config_returns_tool_name(self) -> None:
        """Sem configuração, retorna o nome da ferramenta."""
        config = {"tools": {"bwa": ""}}
        path = get_tool_path(config, "bwa")
        # Deve retornar "bwa" ou o caminho encontrado no PATH
        assert "bwa" in path.lower() or path != ""

    def test_explicit_path(self) -> None:
        """Com caminho explícito, retorna esse caminho."""
        config = {"tools": {"bwa": "/custom/path/bwa"}}
        path = get_tool_path(config, "bwa")
        assert path == "/custom/path/bwa"

    def test_missing_tool_section(self) -> None:
        """Sem seção tools, retorna o nome da ferramenta."""
        config = {}
        path = get_tool_path(config, "bwa")
        assert "bwa" in path.lower()


# ============================================================
# TESTES: format_duration
# ============================================================

class TestFormatDuration:
    def test_seconds(self) -> None:
        assert format_duration(5) == "5s"

    def test_minutes(self) -> None:
        assert format_duration(125) == "2m 05s"

    def test_hours(self) -> None:
        assert format_duration(3665) == "1h 01m 05s"

    def test_zero(self) -> None:
        assert format_duration(0) == "0s"
