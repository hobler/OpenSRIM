# UI Logging Bus

This folder provides a tiny, Qt-free logging API that simulation code can use to write messages into the UI log.

## For simulators (recommended)

```python
from ui.logging import log

log("OpenTRIM: simulation started")
log("OpenTRIM: finished")
```

Notes:
- `log(...)` accepts anything; it is converted to `str(...)`.
- The UI turns messages into timestamped entries via `AppState.add_log(...)`.
- You can call `log(...)` from background threads; the UI subscribes via a Qt signal bridge to stay thread-safe.

## Optional: subscribing to logs

If you want to consume log messages (e.g. tests, CLI runner), you can subscribe a handler:

```python
from ui.logging import subscribe, unsubscribe, log

def handler(message: str) -> None:
    print("LOG:", message)

subscribe(handler)
log("Hello")
unsubscribe(handler)
```

Important: subscriber handlers run in the thread that called `log(...)`. If your handler touches Qt widgets, you must
forward to the Qt main thread yourself (the UI main window already does this).

