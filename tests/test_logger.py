# ============================================================
# Testes — Sistema de Logging
# ============================================================
"""Testes para o módulo de logging do pipeline."""

import sys
import logging
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

from common.logger import (
    setup_pipeline_logger,
    setup_sample_logger,
    get_sample_adapter,
    StepTimer
)


class TestLogger:
    @pytest.fixture(autouse=True)
    def clean_loggers(self):
        """Limpa handlers dos loggers para garantir isolamento entre testes."""
        yield
        for name in logging.root.manager.loggerDict:
            if name == "pipeline" or name.startswith("sample."):
                logging.getLogger(name).handlers.clear()

    def test_setup_pipeline_logger(self, tmp_path: Path) -> None:
        """Testa a inicialização e gravação no logger do pipeline."""
        logs_dir = tmp_path / "logs"
        logger = setup_pipeline_logger(str(logs_dir), "2.0.0")
        
        logger.info("Test message")
        
        log_file = logs_dir / "pipeline.log"
        assert log_file.is_file()
        content = log_file.read_text()
        assert "Test message" in content
        assert "PIPELINE" in content
        
    def test_setup_sample_logger(self, tmp_path: Path) -> None:
        """Testa o logger dedicado para amostras."""
        sample_log_dir = tmp_path / "sample_logs"
        logger = setup_sample_logger("SampleX", str(sample_log_dir))
        
        logger.error("Sample error message", extra={"sample": "SampleX", "step": "CALLER"})
        
        log_file = sample_log_dir / "SampleX.log"
        assert log_file.is_file()
        content = log_file.read_text()
        assert "Sample error message" in content
        assert "SampleX" in content
        
    def test_get_sample_adapter(self, tmp_path: Path) -> None:
        """Testa se o adapter injeta corretamente o sample_name e step_name."""
        logs_dir = tmp_path / "logs"
        results_dir = tmp_path / "results"
        adapter = get_sample_adapter("SampleY", "ALIGNMENT", str(logs_dir), str(results_dir))
        
        adapter.info("Alignment finished")
        
        sample_log_file = results_dir / "SampleY" / "logs" / "SampleY.log"
        assert sample_log_file.is_file()
        content = sample_log_file.read_text()
        assert "SampleY" in content
        assert "ALIGNMENT" in content
        assert "Alignment finished" in content
        
    def test_step_timer(self, tmp_path: Path) -> None:
        """Testa o context manager StepTimer."""
        logs_dir = tmp_path / "logs"
        logger = setup_pipeline_logger(str(logs_dir), "2.0.0")
        
        with StepTimer("TEST_STEP", "SampleZ", logger) as timer:
            import time
            time.sleep(0.1)
            
        assert timer.duration_seconds >= 0.1
        
        log_file = logs_dir / "pipeline.log"
        content = log_file.read_text()
        assert "Iniciando etapa: TEST_STEP" in content
        assert "Etapa conclu" in content
        
    def test_step_timer_exception(self, tmp_path: Path) -> None:
        """Testa se o StepTimer captura corretamente os logs quando há exceção."""
        logs_dir = tmp_path / "logs"
        logger = setup_pipeline_logger(str(logs_dir), "2.0.0")
        
        try:
            with StepTimer("FAIL_STEP", "SampleW", logger):
                raise ValueError("Something went wrong")
        except ValueError:
            pass
            
        log_file = logs_dir / "pipeline.log"
        content = log_file.read_text()
        assert "Etapa falhou" in content
        assert "ValueError: Something went wrong" in content
