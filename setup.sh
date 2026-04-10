#!/usr/bin/env bash
# setup.sh — Install all dependencies for OpenSRIM on Ubuntu (no venv)
set -e

echo "=== OpenSRIM Setup ==="

# ── System packages ─────────────────────────────────────────────────────────
echo "[1/2] Installing system packages..."
sudo apt-get update -qq
sudo apt-get install -y \
    python3-pip \
    python3-dev \
    libgl1 \
    libxcb-cursor0 \
    libxkbcommon-x11-0 \
    libxcb-icccm4 \
    libxcb-image0 \
    libxcb-keysyms1 \
    libxcb-randr0 \
    libxcb-render-util0 \
    libxcb-shape0

# ── Python packages ──────────────────────────────────────────────────────────
# numpy, scipy, matplotlib are already present as Ubuntu system packages;
# pip will skip them if up-to-date.
echo "[2/2] Installing Python packages..."
pip3 install --break-system-packages \
    PyQt6 \
    numpy \
    scipy \
    matplotlib \
    numba \
    tomli

echo ""
echo "=== Setup complete. Run the app with: ==="
echo "  python3 main_window.py"
