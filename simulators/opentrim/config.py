import sys
import os

#os.environ["NUMBA_DISABLE_JIT"] = "1"
#os.environ["NUMBA_FULL_TRACEBACK"] = "1"
#os.environ["NUMBA_DEBUG"] = "1"  # usually too much output
#os.environ["NUMBA_BOUNDSCHECK"] = "1"
#os.environ["NUMBA_DISABLE_ARRAY_EXPR_REWRITE"] = "1"
#os.environ["NUMBA_OPT"] = "0"

if getattr(sys, "frozen", False):
    # Different caching location when using PyInstaller
    os.environ["NUMBA_CACHE_LOCATOR_CLASSES"] = "UserWideCacheLocator"

# Caching disabled unless packed using PyInstaller
ENABLE_CACHING = getattr(sys, "frozen", False)  
# Parallel processing in Numba
PARALLEL = True
DEBUG = False
#PARALLEL = False
#DEBUG = True
# Enable profiling of jitted functions if called with a profiler
#os.environ["NUMBA_ENABLE_PROFILING"] = "1"
