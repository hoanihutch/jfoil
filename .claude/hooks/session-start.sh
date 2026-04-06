#!/bin/bash
set -euo pipefail

# Only run in remote (Claude Code on the web) environments
if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

JULIA_VERSION="1.11.3"
JULIA_MAJOR_MINOR="1.11"
JULIA_INSTALL_DIR="/opt/julia-${JULIA_VERSION}"
JULIA_TARBALL="julia-${JULIA_VERSION}-linux-x86_64.tar.gz"
JULIA_URL="https://julialang-s3.julialang.org/bin/linux/x64/${JULIA_MAJOR_MINOR}/${JULIA_TARBALL}"
PROJECT_DIR="${CLAUDE_PROJECT_DIR:-/home/user/jfoil}"

# ── Install Julia if not already available ──
if ! command -v julia &>/dev/null; then
  echo "Installing Julia ${JULIA_VERSION}..."

  # Try primary S3 URL, fall back to GitHub releases
  if wget -q --timeout=30 "${JULIA_URL}" -O "/tmp/${JULIA_TARBALL}" 2>/dev/null; then
    echo "Downloaded from julialang-s3."
  else
    JULIA_GH_URL="https://github.com/JuliaLang/julia/releases/download/v${JULIA_VERSION}/${JULIA_TARBALL}"
    if wget -q --timeout=30 "${JULIA_GH_URL}" -O "/tmp/${JULIA_TARBALL}" 2>/dev/null; then
      echo "Downloaded from GitHub releases."
    else
      echo "WARNING: Could not download Julia ${JULIA_VERSION}. Julia hosts may be blocked by proxy."
      echo "WARNING: Julia tests will not be available this session."
      # Continue without Julia — don't fail the hook
      JULIA_INSTALL_FAILED=true
    fi
  fi

  if [ "${JULIA_INSTALL_FAILED:-}" != "true" ]; then
    tar xzf "/tmp/${JULIA_TARBALL}" -C /opt/
    ln -sf "${JULIA_INSTALL_DIR}/bin/julia" /usr/local/bin/julia
    rm -f "/tmp/${JULIA_TARBALL}"
    echo "Julia $(julia --version) installed."
  fi
else
  echo "Julia already installed: $(julia --version)"
fi

# ── Instantiate Julia project dependencies ──
if command -v julia &>/dev/null; then
  echo "Instantiating Julia project..."
  julia --project="${PROJECT_DIR}" -e 'using Pkg; Pkg.instantiate()'
fi

# ── Install Python dependencies for reference testing ──
if ! python3 -c "import numpy, scipy, matplotlib" &>/dev/null; then
  echo "Installing Python dependencies (numpy, scipy, matplotlib)..."
  pip install --quiet numpy scipy matplotlib
else
  echo "Python dependencies already installed."
fi

# ── Export Julia to PATH for the session ──
if [ -d "${JULIA_INSTALL_DIR}/bin" ]; then
  echo "export PATH=\"${JULIA_INSTALL_DIR}/bin:\$PATH\"" >> "${CLAUDE_ENV_FILE:-/dev/null}"
fi

echo "Session start hook complete."
