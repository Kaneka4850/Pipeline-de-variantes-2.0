# ============================================================
# Pipeline de Variantes Germinativas — Sistema de Logging
# ============================================================
"""
Sistema de logging profissional para o pipeline.

Registra eventos com:
- Timestamp, nível, amostra, etapa, duração
- Usuário, hostname, versão do pipeline
- Comando executado, código de retorno, mensagem de erro
- Log por amostra + log geral do pipeline
"""

import logging
import os
import time
from pathlib import Path
from typing import Any

from .utils import format_duration, get_hostname, get_username, ensure_dir


# ============================================================
# FORMATO DO LOG
# ============================================================

LOG_FORMAT = (
    "[%(asctime)s] [%(levelname)-8s] [%(sample)-20s] [%(step)-20s] %(message)s"
)
LOG_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


class SampleLogAdapter(logging.LoggerAdapter):
    """Adapter que adiciona campos de amostra e etapa ao log."""

    def process(
        self, msg: str, kwargs: dict[str, Any]
    ) -> tuple[str, dict[str, Any]]:
        kwargs.setdefault("extra", {})
        kwargs["extra"]["sample"] = self.extra.get("sample", "PIPELINE")
        kwargs["extra"]["step"] = self.extra.get("step", "GERAL")
        return msg, kwargs


class SampleStepFilter(logging.Filter):
    """Filtro que garante que campos sample e step existam no record."""

    def filter(self, record: logging.LogRecord) -> bool:
        if not hasattr(record, "sample"):
            record.sample = "PIPELINE"  # type: ignore[attr-defined]
        if not hasattr(record, "step"):
            record.step = "GERAL"  # type: ignore[attr-defined]
        return True


# ============================================================
# CONFIGURAÇÃO DOS LOGGERS
# ============================================================

def setup_pipeline_logger(
    logs_dir: str,
    pipeline_version: str,
    level: int = logging.INFO,
) -> logging.Logger:
    """Configura o logger geral do pipeline.

    Args:
        logs_dir: Diretório para o arquivo de log geral.
        pipeline_version: Versão do pipeline para registro.
        level: Nível de logging.

    Returns:
        Logger configurado.
    """
    ensure_dir(logs_dir)
    log_file = Path(logs_dir) / "pipeline.log"

    logger = logging.getLogger("pipeline")
    logger.setLevel(level)

    # Evita handlers duplicados em re-importações
    if logger.handlers:
        return logger

    # Handler para arquivo
    file_handler = logging.FileHandler(str(log_file), mode="a", encoding="utf-8")
    file_handler.setLevel(level)
    file_handler.setFormatter(logging.Formatter(LOG_FORMAT, LOG_DATE_FORMAT))
    file_handler.addFilter(SampleStepFilter())
    logger.addHandler(file_handler)

    # Handler para console (menos verboso)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(
        logging.Formatter(
            "[%(asctime)s] [%(levelname)-8s] %(message)s", LOG_DATE_FORMAT
        )
    )
    console_handler.addFilter(SampleStepFilter())
    logger.addHandler(console_handler)

    # Mensagem inicial
    logger.info(
        f"Pipeline v{pipeline_version} iniciado | "
        f"Host: {get_hostname()} | "
        f"Usuário: {get_username()}",
        extra={"sample": "PIPELINE", "step": "INIT"},
    )

    return logger


def setup_sample_logger(
    sample_name: str,
    sample_log_dir: str,
    level: int = logging.DEBUG,
) -> logging.Logger:
    """Configura logger dedicado para uma amostra.

    Args:
        sample_name: Nome identificador da amostra.
        sample_log_dir: Diretório para o log da amostra.
        level: Nível de logging.

    Returns:
        Logger configurado para a amostra.
    """
    ensure_dir(sample_log_dir)
    log_file = Path(sample_log_dir) / f"{sample_name}.log"

    logger = logging.getLogger(f"sample.{sample_name}")
    logger.setLevel(level)

    if logger.handlers:
        return logger

    handler = logging.FileHandler(str(log_file), mode="a", encoding="utf-8")
    handler.setLevel(level)
    handler.setFormatter(logging.Formatter(LOG_FORMAT, LOG_DATE_FORMAT))
    handler.addFilter(SampleStepFilter())
    logger.addHandler(handler)

    return logger


def get_sample_adapter(
    sample_name: str,
    step: str,
    logs_dir: str,
    results_dir: str,
) -> SampleLogAdapter:
    """Retorna um LogAdapter para uma amostra e etapa específica.

    Registra em ambos os loggers: geral do pipeline e da amostra.

    Args:
        sample_name: Nome da amostra.
        step: Nome da etapa (ex: "ALIGNMENT", "CALLER").
        logs_dir: Diretório de logs gerais.
        results_dir: Diretório de resultados.

    Returns:
        SampleLogAdapter configurado.
    """
    # Garante que o logger da amostra existe
    sample_log_dir = os.path.join(results_dir, sample_name, "logs")
    sample_logger = setup_sample_logger(sample_name, sample_log_dir)

    # Adiciona os handlers do logger geral ao logger da amostra para registrar em ambos
    pipeline_logger = logging.getLogger("pipeline")
    for handler in pipeline_logger.handlers:
        if handler not in sample_logger.handlers:
            sample_logger.addHandler(handler)

    return SampleLogAdapter(
        sample_logger, {"sample": sample_name, "step": step}
    )


# ============================================================
# CONTEXT MANAGER PARA CRONOMETRAR ETAPAS
# ============================================================

class StepTimer:
    """Context manager que cronometra e registra uma etapa.

    Uso:
        with StepTimer("ALIGNMENT", sample_name, logger):
            # ... código da etapa ...

    Registra automaticamente:
    - Início e fim da etapa
    - Duração total
    - Status (sucesso/erro)
    - Mensagem de erro (se houver)
    """

    def __init__(
        self,
        step_name: str,
        sample_name: str,
        logger: logging.Logger | SampleLogAdapter,
    ) -> None:
        self.step_name = step_name
        self.sample_name = sample_name
        self.logger = logger
        self.start_time: float = 0.0
        self.end_time: float = 0.0

    def __enter__(self) -> "StepTimer":
        self.start_time = time.time()
        self.logger.info(
            f"Iniciando etapa: {self.step_name}",
            extra={"sample": self.sample_name, "step": self.step_name},
        )
        return self

    def __exit__(
        self,
        exc_type: type | None,
        exc_val: BaseException | None,
        exc_tb: Any,
    ) -> bool:
        self.end_time = time.time()
        duration = format_duration(self.end_time - self.start_time)

        if exc_type is None:
            self.logger.info(
                f"Etapa concluída com sucesso ({duration})",
                extra={"sample": self.sample_name, "step": self.step_name},
            )
        else:
            self.logger.error(
                f"Etapa falhou após {duration}: {exc_type.__name__}: {exc_val}",
                extra={"sample": self.sample_name, "step": self.step_name},
            )
        return False  # Não suprime exceções

    @property
    def duration_seconds(self) -> float:
        """Retorna a duração em segundos."""
        if self.end_time > 0:
            return self.end_time - self.start_time
        return time.time() - self.start_time
