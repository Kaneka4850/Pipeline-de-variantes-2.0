import pytest
from scripts.common.database import PipelineDatabase

@pytest.fixture
def temp_db(tmp_path):
    """Fixture que fornece um banco de dados SQLite temporário."""
    db_path = tmp_path / "pipeline.db"
    db = PipelineDatabase(str(db_path))
    yield db
    db.close()

@pytest.fixture
def sample_config(tmp_path):
    """Fixture com um dicionário de configuração simulado."""
    return {
        "pipeline": {"version": "2.1.0", "annotation": True},
        "reference": {
            "genome": str(tmp_path / "ref.fa"),
            "name": "GRCh38",
            "version": "110",
        },
        "vep": {
            "cache_dir": str(tmp_path / "vep_cache"),
            "assembly": "GRCh38",
        },
        "dirs": {
            "results": str(tmp_path / "results"),
            "logs": str(tmp_path / "logs"),
            "database": str(tmp_path / "database"),
        },
        "tools": {
            "bwa": "bwa",
            "gatk": "gatk",
            "samtools": "samtools",
        }
    }
