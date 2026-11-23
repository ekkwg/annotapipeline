"""configure logging for the command-line adapters"""

from __future__ import annotations

import logging
import sys
from typing import TextIO

LOG_FORMAT = "%(asctime)s | %(name)s | %(levelname)s | %(message)s"
DATE_FORMAT = "%d/%m/%Y %H:%M:%S"


def configure(
    filename: str, *, name: str | None = None, level: int = logging.INFO
) -> logging.Logger:
    """write log records to ``filename`` using the legacy layout"""
    logging.basicConfig(
        level=level,
        format=LOG_FORMAT,
        datefmt=DATE_FORMAT,
        filename=filename,
        filemode="a",
        encoding="utf-8",
    )
    return logging.getLogger(name)


class StreamToLogger:
    """file-like stream that forwards complete lines to a logger"""

    def __init__(self, logger: logging.Logger, level: int = logging.ERROR) -> None:
        self.logger = logger
        self.level = level
        self.buffer = ""

    def write(self, text: str) -> None:
        pending = self.buffer + text
        lines = pending.splitlines(keepends=True)
        self.buffer = ""
        if lines and not lines[-1].endswith(("\n", "\r")):
            self.buffer = lines.pop()
        for line in lines:
            self.logger.log(self.level, line.rstrip())

    def flush(self) -> None:
        if self.buffer:
            self.logger.log(self.level, self.buffer.rstrip())
        self.buffer = ""

    def isatty(self) -> bool:
        return False


def redirect_stderr(logger: logging.Logger) -> TextIO:
    """replace ``sys.stderr`` with a logger-backed stream"""
    original = sys.stderr
    sys.stderr = StreamToLogger(logger)
    return original
