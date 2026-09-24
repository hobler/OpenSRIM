@echo off
setlocal

python --version >nul 2>&1
if errorlevel 1 (
    echo Python interpreter missing, download it here: https://www.python.org/downloads/
    pause
    exit /b 1
)

python -c "import sys; raise SystemExit(sys.version_info < (3, 11))" >nul 2>&1
if errorlevel 1 (
    echo Outdated Python version, download a new one here: https://www.python.org/downloads/
    pause
    exit /b 1
)

if not exist "venv\" (
    echo Installing dependencies...
    python -m venv venv
    venv\Scripts\python.exe -m pip install -r requirements.txt
)

call venv\Scripts\activate.bat
python main_window.py
