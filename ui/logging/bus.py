from __future__ import annotations

from collections.abc import Callable
import sys
from threading import RLock
import traceback
from typing import Any

LogHandler = Callable[[str], Any]


class LogBus:
    """Lightweight publish/subscribe log bus.

    - Importable from non-UI modules (no Qt dependency).
    - Thread-safe subscription management.
    - Emits to all subscribers; subscriber exceptions are isolated.
    """

    def __init__(self) -> None:
        self._lock = RLock()
        self._subscribers: list[LogHandler] = []

    def subscribe(self, handler: LogHandler) -> None:
        if not callable(handler):
            raise TypeError("log handler must be callable")
        with self._lock:
            if handler not in self._subscribers:
                self._subscribers.append(handler)

    def unsubscribe(self, handler: LogHandler) -> None:
        with self._lock:
            self._subscribers = [h for h in self._subscribers if h != handler]

    def emit(self, message: str) -> None:
        msg = str(message)
        with self._lock:
            subscribers = list(self._subscribers)
        for handler in subscribers:
            try:
                handler(msg)
            except Exception:
                traceback.print_exc(file=sys.stderr)
