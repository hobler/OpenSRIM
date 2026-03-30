from __future__ import annotations

from .bus import LogBus, LogHandler

_BUS = LogBus()


def subscribe(handler: LogHandler) -> None:
    _BUS.subscribe(handler)


def unsubscribe(handler: LogHandler) -> None:
    _BUS.unsubscribe(handler)


def log(message: str) -> None:
    """Emit a log message to all subscribers (e.g. the UI)."""

    _BUS.emit(message)


__all__ = [
    "LogBus",
    "LogHandler",
    "log",
    "subscribe",
    "unsubscribe",
]
