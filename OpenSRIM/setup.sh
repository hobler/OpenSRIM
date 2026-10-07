#!/usr/bin/env bash
# setup.sh — Install all dependencies for OpenSRIM on Ubuntu (uses a venv)
set -e

echo "=== OpenSRIM Setup ==="

# ── System packages ─────────────────────────────────────────────────────────
# Qt/X11 libraries PyQt6 needs at runtime; these aren't pip-installable.
echo "[1/3] Installing system packages..."
sudo apt-get update -qq
sudo apt-get install -y \
    python3-pip \
    python3-dev \
    python3-venv \
    libgl1 \
    libxcb-cursor0 \
    libxkbcommon-x11-0 \
    libxcb-icccm4 \
    libxcb-image0 \
    libxcb-keysyms1 \
    libxcb-randr0 \
    libxcb-render-util0 \
    libxcb-shape0

# ── Virtual environment + Python packages ───────────────────────────────────
echo "[2/3] Creating virtual environment (./venv)..."
python3 -m venv venv

echo "[3/3] Installing Python packages into ./venv..."
./venv/bin/pip install --upgrade pip
./venv/bin/pip install -r requirements.txt

echo ""
echo "=== Setup complete. Run the app with: ==="
echo "  source venv/bin/activate"
echo "  python3 main_window.py"
