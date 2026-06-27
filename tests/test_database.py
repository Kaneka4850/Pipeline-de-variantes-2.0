# ============================================================
# Testes — Banco de Dados
# ============================================================
"""Testes para o módulo de banco de dados SQLite."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

from common.database import PipelineDatabase


@pytest.fixture
def db(tmp_path: Path) -> PipelineDatabase:
    """Cria banco de dados temporário para testes."""
    db_path = str(tmp_path / "test_pipeline.db")
    database = PipelineDatabase(db_path)
    yield database
    database.close()


class TestPipelineDatabase:
    def test_create_database(self, db: PipelineDatabase) -> None:
        """Banco é criado com todas as tabelas."""
        tables = db.conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"
        ).fetchall()
        table_names = {row["name"] for row in tables}

        assert "pipeline_runs" in table_names
        assert "samples" in table_names
        assert "steps" in table_names
        assert "files" in table_names
        assert "tools" in table_names

    def test_start_and_finish_run(self, db: PipelineDatabase) -> None:
        """Registra e finaliza uma execução."""
        run_id = db.start_run("2.0.0", "config.yaml", total_samples=3)
        assert run_id > 0

        db.finish_run(run_id, "SUCCESS")

        row = db.conn.execute(
            "SELECT * FROM pipeline_runs WHERE id = ?", (run_id,)
        ).fetchone()
        assert row["status"] == "SUCCESS"
        assert row["end_time"] is not None

    def test_start_and_finish_sample(self, db: PipelineDatabase) -> None:
        """Registra e finaliza uma amostra."""
        run_id = db.start_run("2.0.0")
        sample_id = db.start_sample(
            run_id, "Sample01", "PE", "r1.fq.gz", "r2.fq.gz"
        )
        assert sample_id > 0

        db.finish_sample(sample_id, "SUCCESS")

        row = db.conn.execute(
            "SELECT * FROM samples WHERE id = ?", (sample_id,)
        ).fetchone()
        assert row["status"] == "SUCCESS"
        assert row["total_duration"] is not None

    def test_steps_lifecycle(self, db: PipelineDatabase) -> None:
        """Registra etapas com sucesso e erro."""
        run_id = db.start_run("2.0.0")
        sample_id = db.start_sample(run_id, "Sample01", "PE", "r1.fq.gz")

        # Etapa com sucesso
        step_id = db.start_step(sample_id, "ALIGNMENT", "bwa mem ...")
        db.finish_step(step_id, return_code=0)

        # Etapa com erro
        step_id2 = db.start_step(sample_id, "CALLER", "gatk HaplotypeCaller ...")
        db.finish_step(step_id2, return_code=1, error_message="Out of memory")

        steps = db.query_steps(sample_id)
        assert len(steps) == 2
        assert steps[0]["status"] == "SUCCESS"
        assert steps[1]["status"] == "ERROR"
        assert steps[1]["error_message"] == "Out of memory"

    def test_register_file(self, db: PipelineDatabase, tmp_path: Path) -> None:
        """Registra um arquivo com checksum."""
        # Cria arquivo de teste
        test_file = tmp_path / "test.bam"
        test_file.write_bytes(b"fake bam content")

        run_id = db.start_run("2.0.0")
        sample_id = db.start_sample(run_id, "Sample01", "SE", "r1.fq.gz")

        file_id = db.register_file(
            sample_id, str(test_file), "BAM", compute_checksum=True
        )
        assert file_id > 0

        row = db.conn.execute(
            "SELECT * FROM files WHERE id = ?", (file_id,)
        ).fetchone()
        assert row["file_type"] == "BAM"
        assert row["checksum_md5"] is not None
        assert row["file_size"] > 0

    def test_register_tool(self, db: PipelineDatabase) -> None:
        """Registra ferramenta e versão."""
        run_id = db.start_run("2.0.0")
        db.register_tool(run_id, "bwa", "0.7.17", "/usr/bin/bwa")

        rows = db.conn.execute(
            "SELECT * FROM tools WHERE run_id = ?", (run_id,)
        ).fetchall()
        assert len(rows) == 1
        assert rows[0]["tool_name"] == "bwa"
        assert rows[0]["tool_version"] == "0.7.17"

    def test_query_sample(self, db: PipelineDatabase) -> None:
        """Consulta histórico de uma amostra."""
        run_id = db.start_run("2.0.0")
        db.start_sample(run_id, "SampleX", "PE", "r1.fq.gz", "r2.fq.gz")

        results = db.query_sample("SampleX")
        assert len(results) == 1
        assert results[0]["sample_name"] == "SampleX"

    def test_query_all_samples(self, db: PipelineDatabase) -> None:
        """Lista todas as amostras."""
        run_id = db.start_run("2.0.0")
        db.start_sample(run_id, "A", "PE", "a_r1.fq.gz")
        db.start_sample(run_id, "B", "SE", "b.fq.gz")

        results = db.query_all_samples()
        assert len(results) == 2

    def test_get_summary(self, db: PipelineDatabase) -> None:
        """Resumo geral do banco."""
        run_id = db.start_run("2.0.0")
        s1 = db.start_sample(run_id, "S1", "PE", "r1.fq.gz")
        s2 = db.start_sample(run_id, "S2", "SE", "r1.fq.gz")
        db.finish_sample(s1, "SUCCESS")
        db.finish_sample(s2, "ERROR", "Failed")

        summary = db.get_summary()
        assert summary["total_runs"] == 1
        assert summary["total_samples"] == 2
        assert summary["successful"] == 1
        assert summary["errors"] == 1

    def test_sample_error_flow(self, db: PipelineDatabase) -> None:
        """Fluxo completo com erro."""
        run_id = db.start_run("2.0.0")
        sample_id = db.start_sample(run_id, "BadSample", "PE", "r1.fq.gz")
        step_id = db.start_step(sample_id, "ALIGNMENT")
        db.finish_step(step_id, return_code=137, error_message="Killed (OOM)")
        db.finish_sample(sample_id, "ERROR", "Pipeline falhou na etapa ALIGNMENT")
        db.finish_run(run_id, "PARTIAL")

        row = db.conn.execute(
            "SELECT * FROM pipeline_runs WHERE id = ?", (run_id,)
        ).fetchone()
        assert row["status"] == "PARTIAL"
