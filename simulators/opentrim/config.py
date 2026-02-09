import sys
import os
# os.environ["NUMBA_DISABLE_JIT"] = "1"
if getattr(sys, 'frozen', False):   # Different caching location when using PyInstaller
    os.environ["NUMBA_CACHE_LOCATOR_CLASSES"] = "UserWideCacheLocator"

ENABLE_CACHING = getattr(sys, 'frozen', False)  # Caching disabled unless packed using PyInstaller