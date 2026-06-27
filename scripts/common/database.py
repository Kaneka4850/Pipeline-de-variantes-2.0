# ============================================================
# Pipeline de Variantes Germinativas — Banco de Dados SQLite
# ============================================================
"""
Interface para banco de dados SQLite de rastreabilidade.

Registra execuções do pipeline, etapas por amostra, arquivos
gerados, versões de ferramentas e checksums.
"""

import hashlib
import os
import sqlite3
from datetime import datetime
from pathlib import Path
from typing import Any

from .utils import ensure_dir, get_hostname, get_username


# ============================================================
# SCHEMA DO BANCO
# ============================================================

SCHEMA_SQL = """
-- Execuções do pipeline
CREATE TABLE IF NOT EXISTS pipeline_runs (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    version         TEXT    NOT NULL,
    hostname        TEXT    NOT NULL,
    username        TEXT    NOT NULL,
    start_time      TEXT    NOT NULL,
    end_time        TEXT,
    status          TEXT    DEFAULT 'RUNNING',
    total_samples   INTEGER DEFAULT 0,
    config_path     TEXT,
    os_info         TEXT,
    cpu_count       INTEGER,
    memory_total_mb INTEGER,
    config_json     TEXT
);

-- Amostras processadas
CREATE TABLE IF NOT EXISTS samples (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id          INTEGER NOT NULL,
    sample_name     TEXT    NOT NULL,
    seq_type        TEXT    NOT NULL CHECK(seq_type IN ('PE', 'SE')),
    r1_path         TEXT    NOT NULL,
    r2_path         TEXT,
    start_time      TEXT    NOT NULL,
    end_time        TEXT,
    status          TEXT    DEFAULT 'RUNNING',
    total_duration  REAL,
    error_message   TEXT,
    FOREIGN KEY (run_id) REFERENCES pipeline_runs(id)
);

-- Etapas executadas por amostra
CREATE TABLE IF NOT EXISTS steps (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id       INTEGER NOT NULL,
    step_name       TEXT    NOT NULL,
    command         TEXT,
    start_time      TEXT    NOT NULL,
    end_time        TEXT,
    duration        REAL,
    status          TEXT    DEFAULT 'RUNNING',
    return_code     INTEGER,
    error_message   TEXT,
    log_file        TEXT,
    FOREIGN KEY (sample_id) REFERENCES samples(id)
);

-- Arquivos gerados
CREATE TABLE IF NOT EXISTS files (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id       INTEGER NOT NULL,
    file_path       TEXT    NOT NULL,
    file_type       TEXT    NOT NULL,
    file_size       INTEGER,
    checksum_md5    TEXT,
    created_at      TEXT    NOT NULL,
    FOREIGN KEY (sample_id) REFERENCES samples(id)
);

-- Ferramentas e versões
CREATE TABLE IF NOT EXISTS tools (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id          INTEGER NOT NULL,
    tool_name       TEXT    NOT NULL,
    tool_version    TEXT    NOT NULL,
    tool_path       TEXT,
    FOREIGN KEY (run_id) REFERENCES pipeline_runs(id)
);

-- Índices para consultas frequentes
CREATE INDEX IF NOT EXISTS idx_samples_name ON samples(sample_name);
CREATE INDEX IF NOT EXISTS idx_samples_status ON samples(status);
CREATE INDEX IF NOT EXISTS idx_steps_sample ON steps(sample_id);
CREATE INDEX IF NOT EXISTS idx_files_sample ON files(sample_id);

-- ==========================================
-- NOVAS TABELAS DE REPRODUTIBILIDADE E VERSIONAMENTO
-- ==========================================

-- Referências genômicas utilizadas
CREATE TABLE IF NOT EXISTS reference_metadata (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id          INTEGER NOT NULL,
    name            TEXT    NOT NULL,
    version         TEXT,
    organism        TEXT    DEFAULT 'Homo sapiens',
    assembly        TEXT    NOT NULL,
    source          TEXT,
    fasta_path      TEXT    NOT NULL,
    checksum_sha256 TEXT,
    install_date    TEXT,
    indices_status  TEXT,
    integrity_ok    BOOLEAN DEFAULT 1,
    FOREIGN KEY (run_id) REFERENCES pipeline_runs(id)
);

-- Informações do VEP por execução
CREATE TABLE IF NOT EXISTS vep_info (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id          INTEGER NOT NULL,
    vep_version     TEXT,
    cache_version   TEXT,
    cache_dir       TEXT,
    assembly        TEXT,
    fasta_version   TEXT,
    plugins_dir     TEXT,
    plugins_list    TEXT,
    FOREIGN KEY (run_id) REFERENCES pipeline_runs(id)
);

-- Snapshot do ambiente de execução
CREATE TABLE IF NOT EXISTS environment_snapshots (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id          INTEGER NOT NULL,
    os_info         TEXT,
    hostname        TEXT,
    username        TEXT,
    cpu_count       INTEGER,
    memory_total_mb INTEGER,
    memory_free_mb  INTEGER,
    snapshot_json   TEXT,
    FOREIGN KEY (run_id) REFERENCES pipeline_runs(id)
);

-- Resultados do preflight
CREATE TABLE IF NOT EXISTS preflight_results (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id          INTEGER NOT NULL,
    check_name      TEXT    NOT NULL,
    status          TEXT    NOT NULL CHECK(status IN ('PASS', 'WARN', 'FAIL')),
    message         TEXT,
    remediation     TEXT,
    FOREIGN KEY (run_id) REFERENCES pipeline_runs(id)
);
"""


# ============================================================
# CLASSE PRINCIPAL
# ============================================================

class PipelineDatabase:
    """Interface para o banco de dados SQLite de rastreabilidade.

    Uso:
        db = PipelineDatabase("database/pipeline.db")
        run_id = db.start_run("2.0.0", "config/config.yaml")
        sample_id = db.start_sample(run_id, "Sample01", "PE", "r1.fq.gz", "r2.fq.gz")
        step_id = db.start_step(sample_id, "ALIGNMENT", "bwa mem ...")
        db.finish_step(step_id, 0)
        db.finish_sample(sample_id)
        db.finish_run(run_id)
    """

    def __init__(self, db_path: str) -> None:
        """Inicializa a conexão com o banco.

        Args:
            db_path: Caminho para o arquivo .db (criado automaticamente).
        """
        ensure_dir(Path(db_path).parent)
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)
        self.conn.row_factory = sqlite3.Row
        self.conn.execute("PRAGMA journal_mode=WAL")
        self.conn.execute("PRAGMA foreign_keys=ON")
        self._init_schema()

    def _init_schema(self) -> None:
        """Cria as tabelas se não existirem e aplica migrações se necessário."""
        self.conn.executescript(SCHEMA_SQL)
        # Adiciona colunas extras se o banco já existia antes da v2.1
        try:
            self.conn.execute("ALTER TABLE pipeline_runs ADD COLUMN os_info TEXT")
            self.conn.execute("ALTER TABLE pipeline_runs ADD COLUMN cpu_count INTEGER")
            self.conn.execute("ALTER TABLE pipeline_runs ADD COLUMN memory_total_mb INTEGER")
            self.conn.execute("ALTER TABLE pipeline_runs ADD COLUMN config_json TEXT")
        except sqlite3.OperationalError:
            pass # Colunas já existem
        self.conn.commit()

    def close(self) -> None:
        """Fecha a conexão com o banco."""
        self.conn.close()

    # ---- Pipeline Runs ----

    def start_run(
        self,
        version: str,
        config_path: str = "",
        total_samples: int = 0,
    ) -> int:
        """Registra o início de uma execução do pipeline.

        Returns:
            ID da execução criada.
        """
        cursor = self.conn.execute(
            """INSERT INTO pipeline_runs
               (version, hostname, username, start_time, total_samples, config_path)
               VALUES (?, ?, ?, ?, ?, ?)""",
            (
                version,
                get_hostname(),
                get_username(),
                datetime.now().isoformat(),
                total_samples,
                config_path,
            ),
        )
        self.conn.commit()
        return cursor.lastrowid  # type: ignore[return-value]

    def update_run_environment(
        self,
        run_id: int,
        os_info: str,
        cpu_count: int,
        memory_total_mb: int,
        config_json: str,
    ) -> None:
        """Atualiza a execução com dados do ambiente e configuração."""
        self.conn.execute(
            """UPDATE pipeline_runs
               SET os_info = ?, cpu_count = ?, memory_total_mb = ?, config_json = ?
               WHERE id = ?""",
            (os_info, cpu_count, memory_total_mb, config_json, run_id),
        )
        self.conn.commit()

    def finish_run(self, run_id: int, status: str = "SUCCESS") -> None:
        """Registra o fim de uma execução do pipeline."""
        self.conn.execute(
            """UPDATE pipeline_runs
               SET end_time = ?, status = ?
               WHERE id = ?""",
            (datetime.now().isoformat(), status, run_id),
        )
        self.conn.commit()

    # ---- Amostras ----

    def start_sample(
        self,
        run_id: int,
        sample_name: str,
        seq_type: str,
        r1_path: str,
        r2_path: str = "",
    ) -> int:
        """Registra o início do processamento de uma amostra.

        Returns:
            ID da amostra criada.
        """
        cursor = self.conn.execute(
            """INSERT INTO samples
               (run_id, sample_name, seq_type, r1_path, r2_path, start_time)
               VALUES (?, ?, ?, ?, ?, ?)""",
            (
                run_id,
                sample_name,
                seq_type,
                r1_path,
                r2_path,
                datetime.now().isoformat(),
            ),
        )
        self.conn.commit()
        return cursor.lastrowid  # type: ignore[return-value]

    def finish_sample(
        self,
        sample_id: int,
        status: str = "SUCCESS",
        error_message: str = "",
    ) -> None:
        """Registra o fim do processamento de uma amostra."""
        # Calcula duração total
        row = self.conn.execute(
            "SELECT start_time FROM samples WHERE id = ?", (sample_id,)
        ).fetchone()
        duration = None
        if row:
            start = datetime.fromisoformat(row["start_time"])
            duration = (datetime.now() - start).total_seconds()

        self.conn.execute(
            """UPDATE samples
               SET end_time = ?, status = ?, total_duration = ?, error_message = ?
               WHERE id = ?""",
            (
                datetime.now().isoformat(),
                status,
                duration,
                error_message,
                sample_id,
            ),
        )
        self.conn.commit()

    # ---- Etapas ----

    def start_step(
        self,
        sample_id: int,
        step_name: str,
        command: str = "",
        log_file: str = "",
    ) -> int:
        """Registra o início de uma etapa.

        Returns:
            ID da etapa criada.
        """
        cursor = self.conn.execute(
            """INSERT INTO steps
               (sample_id, step_name, command, start_time, log_file)
               VALUES (?, ?, ?, ?, ?)""",
            (
                sample_id,
                step_name,
                command,
                datetime.now().isoformat(),
                log_file,
            ),
        )
        self.conn.commit()
        return cursor.lastrowid  # type: ignore[return-value]

    def finish_step(
        self,
        step_id: int,
        return_code: int = 0,
        error_message: str = "",
    ) -> None:
        """Registra o fim de uma etapa."""
        row = self.conn.execute(
            "SELECT start_time FROM steps WHERE id = ?", (step_id,)
        ).fetchone()
        duration = None
        if row:
            start = datetime.fromisoformat(row["start_time"])
            duration = (datetime.now() - start).total_seconds()

        status = "SUCCESS" if return_code == 0 else "ERROR"
        self.conn.execute(
            """UPDATE steps
               SET end_time = ?, duration = ?, status = ?,
                   return_code = ?, error_message = ?
               WHERE id = ?""",
            (
                datetime.now().isoformat(),
                duration,
                status,
                return_code,
                error_message,
                step_id,
            ),
        )
        self.conn.commit()

    # ---- Arquivos ----

    def register_file(
        self,
        sample_id: int,
        file_path: str,
        file_type: str,
        compute_checksum: bool = True,
    ) -> int:
        """Registra um arquivo gerado pelo pipeline.

        Args:
            sample_id: ID da amostra.
            file_path: Caminho do arquivo.
            file_type: Tipo (FASTQ, BAM, VCF, etc.).
            compute_checksum: Se True, calcula MD5.

        Returns:
            ID do registro criado.
        """
        file_size = None
        checksum = None

        if os.path.isfile(file_path):
            file_size = os.path.getsize(file_path)
            if compute_checksum:
                checksum = _compute_md5(file_path)

        cursor = self.conn.execute(
            """INSERT INTO files
               (sample_id, file_path, file_type, file_size, checksum_md5, created_at)
               VALUES (?, ?, ?, ?, ?, ?)""",
            (
                sample_id,
                file_path,
                file_type,
                file_size,
                checksum,
                datetime.now().isoformat(),
            ),
        )
        self.conn.commit()
        return cursor.lastrowid  # type: ignore[return-value]

    # ---- Ferramentas ----

    def register_tool(
        self,
        run_id: int,
        tool_name: str,
        tool_version: str,
        tool_path: str = "",
    ) -> None:
        """Registra uma ferramenta e sua versão."""
        self.conn.execute(
            """INSERT INTO tools (run_id, tool_name, tool_version, tool_path)
               VALUES (?, ?, ?, ?)""",
            (run_id, tool_name, tool_version, tool_path),
        )
        self.conn.commit()

    # ---- Versionamento e Preflight ----

    def register_reference(self, run_id: int, ref_info: dict[str, Any]) -> int:
        """Registra metadados da referência genômica utilizada."""
        import json
        indices_json = json.dumps(ref_info.get("indices", {}))
        cursor = self.conn.execute(
            """INSERT INTO reference_metadata
               (run_id, name, version, organism, assembly, source, fasta_path,
                checksum_sha256, install_date, indices_status, integrity_ok)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (
                run_id,
                ref_info.get("name"),
                ref_info.get("version"),
                ref_info.get("organism"),
                ref_info.get("assembly"),
                ref_info.get("source"),
                ref_info.get("fasta_path"),
                ref_info.get("checksum_sha256"),
                ref_info.get("install_date"),
                indices_json,
                1 if ref_info.get("integrity_ok") else 0,
            ),
        )
        self.conn.commit()
        return cursor.lastrowid  # type: ignore[return-value]

    def register_vep_info(self, run_id: int, vep_info: dict[str, Any]) -> int:
        """Registra informações do VEP."""
        import json
        plugins_json = json.dumps(vep_info.get("plugins_installed", []))
        cursor = self.conn.execute(
            """INSERT INTO vep_info
               (run_id, vep_version, cache_version, cache_dir, assembly,
                fasta_version, plugins_dir, plugins_list)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
            (
                run_id,
                vep_info.get("version"),
                vep_info.get("cache_version"),
                vep_info.get("cache_dir"),
                vep_info.get("assembly"),
                vep_info.get("fasta_version"),
                vep_info.get("plugins_dir"),
                plugins_json,
            ),
        )
        self.conn.commit()
        return cursor.lastrowid  # type: ignore[return-value]

    def register_environment(self, run_id: int, snapshot: dict[str, Any]) -> int:
        """Registra snapshot do ambiente."""
        import json
        sys_info = snapshot.get("system", {})
        cursor = self.conn.execute(
            """INSERT INTO environment_snapshots
               (run_id, os_info, hostname, username, cpu_count,
                memory_total_mb, memory_free_mb, snapshot_json)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
            (
                run_id,
                f"{sys_info.get('os_name')} {sys_info.get('os_version')}",
                sys_info.get("hostname"),
                sys_info.get("username"),
                sys_info.get("cpu_count"),
                sys_info.get("memory_total_mb"),
                sys_info.get("memory_available_mb"),
                json.dumps(snapshot),
            ),
        )
        self.conn.commit()
        return cursor.lastrowid  # type: ignore[return-value]

    def register_preflight(self, run_id: int, results: list[dict[str, Any]]) -> None:
        """Registra resultados do preflight."""
        for r in results:
            self.conn.execute(
                """INSERT INTO preflight_results
                   (run_id, check_name, status, message, remediation)
                   VALUES (?, ?, ?, ?, ?)""",
                (
                    run_id,
                    r.get("check_name"),
                    r.get("status"),
                    r.get("message"),
                    r.get("remediation"),
                ),
            )
        self.conn.commit()

    # ---- Consultas ----

    def query_sample(self, sample_name: str) -> list[dict[str, Any]]:
        """Consulta o histórico de uma amostra.

        Returns:
            Lista de dicionários com dados das execuções da amostra.
        """
        rows = self.conn.execute(
            """SELECT s.*, pr.version as pipeline_version
               FROM samples s
               JOIN pipeline_runs pr ON s.run_id = pr.id
               WHERE s.sample_name = ?
               ORDER BY s.start_time DESC""",
            (sample_name,),
        ).fetchall()
        return [dict(row) for row in rows]

    def query_steps(self, sample_id: int) -> list[dict[str, Any]]:
        """Consulta as etapas executadas para uma amostra.

        Returns:
            Lista de dicionários com dados das etapas.
        """
        rows = self.conn.execute(
            """SELECT * FROM steps WHERE sample_id = ? ORDER BY id""",
            (sample_id,),
        ).fetchall()
        return [dict(row) for row in rows]

    def query_all_samples(self) -> list[dict[str, Any]]:
        """Consulta todas as amostras registradas.

        Returns:
            Lista de dicionários com dados resumidos.
        """
        rows = self.conn.execute(
            """SELECT s.sample_name, s.seq_type, s.status,
                      s.start_time, s.total_duration,
                      pr.version as pipeline_version
               FROM samples s
               JOIN pipeline_runs pr ON s.run_id = pr.id
               ORDER BY s.start_time DESC"""
        ).fetchall()
        return [dict(row) for row in rows]

    def get_summary(self) -> dict[str, Any]:
        """Retorna um resumo geral do banco.

        Returns:
            Dicionário com contagens e estatísticas.
        """
        total_runs = self.conn.execute(
            "SELECT COUNT(*) FROM pipeline_runs"
        ).fetchone()[0]
        total_samples = self.conn.execute(
            "SELECT COUNT(*) FROM samples"
        ).fetchone()[0]
        success = self.conn.execute(
            "SELECT COUNT(*) FROM samples WHERE status = 'SUCCESS'"
        ).fetchone()[0]
        errors = self.conn.execute(
            "SELECT COUNT(*) FROM samples WHERE status = 'ERROR'"
        ).fetchone()[0]

        return {
            "total_runs": total_runs,
            "total_samples": total_samples,
            "successful": success,
            "errors": errors,
        }


# ============================================================
# FUNÇÕES AUXILIARES
# ============================================================

def _compute_md5(file_path: str, chunk_size: int = 8192) -> str:
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
