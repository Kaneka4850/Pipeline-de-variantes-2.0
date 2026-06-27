# ============================================================
# Testes — Integração do Pipeline (Snakemake)
# ============================================================
"""Testes de integração para verificar se o DAG do Snakemake compila corretamente."""

import subprocess
from pathlib import Path
import yaml
import shutil

import pytest


@pytest.fixture
def integration_env(tmp_path: Path) -> dict:
    """Prepara um ambiente temporário com dados mock para rodar o Snakemake."""
    # Cria diretório de projeto temporário
    project_dir = tmp_path / "project"
    project_dir.mkdir()
    
    # Copia os mocks
    samples_dir = project_dir / "samples"
    samples_dir.mkdir()
    (samples_dir / "sample_test_R1.fastq.gz").touch()
    (samples_dir / "sample_test_R2.fastq.gz").touch()
    
    ref_dir = project_dir / "reference"
    ref_dir.mkdir()
    (ref_dir / "mock_ref.fa").touch()
    
    # Cria uma config.yaml para o teste
    config = {
        "pipeline": {
            "name": "Integration Test",
            "version": "1.0",
            "annotation": False
        },
        "dirs": {
            "samples": str(samples_dir),
            "results": str(project_dir / "results"),
            "logs": str(project_dir / "logs"),
            "database": str(project_dir / "database")
        },
        "reference": {
            "genome": str(ref_dir / "mock_ref.fa"),
            "name": "mock"
        },
        "sample_detection": {
            "r1_patterns": ["_R1.fastq.gz"],
            "r2_patterns": ["_R2.fastq.gz"]
        },
        "preflight": {
            "enabled": False # Desabilita o preflight para não falhar por faltar VEP ou espaço real
        },
        "resources": {
            "threads": 1,
            "memory_mb": 1024
        },
        "read_groups": {
            "platform": "illumina",
            "library": "lib1",
            "platform_unit": "pu1"
        }
    }
    
    config_path = project_dir / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config, f)
        
    return {
        "project_dir": project_dir,
        "config_path": config_path
    }


class TestPipelineIntegration:
    @pytest.mark.skipif(shutil.which("snakemake") is None, reason="Snakemake no est instalado no PATH")
    def test_snakemake_dry_run(self, integration_env: dict) -> None:
        """Executa um dry-run do Snakemake para verificar se todas as regras são compiladas sem erros de sintaxe ou I/O."""
        
        # O Snakefile principal fica na raiz do repositório
        repo_root = Path(__file__).resolve().parent.parent
        snakefile_path = repo_root / "Snakefile"
        
        config_path = integration_env["config_path"]
        
        # Executa snakemake no modo dry-run (-n)
        cmd = [
            "snakemake",
            "-s", str(snakefile_path),
            "--configfile", str(config_path),
            "-n",  # Dry run
            "--quiet" # Reduz a saída verbosa
        ]
        
        # Subprocess
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(repo_root))
        
        # Verifica se o Snakemake executou sem erros
        assert result.returncode == 0, f"Snakemake dry-run falhou: {result.stderr}"
        
        # Verifica se o log menciona a criação dos arquivos esperados
        assert "sample_test" in result.stdout or "sample_test" in result.stderr
        
