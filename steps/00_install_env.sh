#!/usr/bin/env bash
# =============================================================================
# 00_install_env.sh — One-time environment setup
# Uses micromamba (a fast, drop-in conda alternative; conda/mamba also work).
# =============================================================================
# !! RUN INTERACTIVELY ON A LOGIN NODE — NOT via sbatch !!
#
#   cd /path/to/phagefactor/
#   bash steps/00_install_env.sh
#
# What this does:
#   1. Installs micromamba binary to ~/.mamba/bin/
#   2. Creates the 'phagefactor' env (python=3.11) with pharokka, phold, foldseek
#   3. Downloads the phold reference database (phold install)
#   4. Creates a separate 'phynteny' env (python=3.10)
#      ⚠️  phynteny requires Python < 3.11 — cannot share the phagefactor env
#   5. Verifies all key tools
#
# Rerunning is safe — each step checks if already done.
# Takes ~20–40 min depending on network speed.
# =============================================================================
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/config.sh"   # for MAMBA_ROOT_PREFIX, MAMBA_ENV

MAMBA_BIN="${MAMBA_ROOT_PREFIX}/bin/micromamba"
PHYNTENY_ENV="phynteny"

# Colours
GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; NC='\033[0m'
ok()   { echo -e "${GREEN}[OK]${NC}  $*"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $*"; }
info() { echo -e "       $*"; }

echo ""
echo "════════════════════════════════════════════════════════"
echo "  phageFACTor — environment setup"
echo "  $(date '+%Y-%m-%d %H:%M:%S')"
echo "════════════════════════════════════════════════════════"
echo ""

# ---------------------------------------------------------------------------
# Step 1 — Install micromamba
# ---------------------------------------------------------------------------
echo "── Step 1: micromamba ──────────────────────────────────"

if [[ -x "${MAMBA_BIN}" ]]; then
    ok "micromamba already installed: $(${MAMBA_BIN} --version)"
else
    info "Downloading micromamba to ${MAMBA_ROOT_PREFIX}/bin/ ..."
    mkdir -p "${MAMBA_ROOT_PREFIX}/bin"
    curl -fsSL "https://micro.mamba.pm/api/micromamba/linux-64/latest" \
        | tar -xvj -C "${MAMBA_ROOT_PREFIX}/bin" --strip-components=1 bin/micromamba
    chmod +x "${MAMBA_BIN}"
    ok "micromamba installed: $(${MAMBA_BIN} --version)"
    "${MAMBA_BIN}" shell init --shell bash --root-prefix "${MAMBA_ROOT_PREFIX}"
    warn "Shell hook added to ~/.bashrc (effective for new shells)."
fi

export MAMBA_ROOT_PREFIX
eval "$("${MAMBA_BIN}" shell hook --shell bash)"
ok "micromamba shell hook loaded"

# ---------------------------------------------------------------------------
# Step 2 — Create the phagefactor env (python=3.11)
# Includes: pharokka, phold, foldseek, and Python dependencies.
# NOTE: phynteny requires Python <3.11 and CANNOT be installed here.
#       It lives in its own 'phynteny' env (Step 4 below).
# ---------------------------------------------------------------------------
echo ""
echo "── Step 2: pharokka environment (python=3.11) ──────────"

if [[ -d "${MAMBA_ROOT_PREFIX}/envs/${MAMBA_ENV}" ]]; then
    ok "Environment '${MAMBA_ENV}' already exists"
    warn "To reinstall: micromamba env remove -n ${MAMBA_ENV} && bash 00_install_env.sh"
else
    info "Creating '${MAMBA_ENV}' (pharokka + phold + foldseek) ..."
    info "This takes 10–20 min."
    micromamba create -y -n "${MAMBA_ENV}" \
        -c conda-forge -c bioconda \
        python=3.11 \
        pharokka \
        foldseek \
        diamond \
        mmseqs2 \
        biopython \
        pandas \
        openpyxl \
        tqdm
    ok "Environment '${MAMBA_ENV}' created."
fi

micromamba activate "${MAMBA_ENV}"

# ---------------------------------------------------------------------------
# Step 3 — Install phold via pip (latest version, avoids bioconda lag)
# ---------------------------------------------------------------------------
echo ""
echo "── Step 3: phold (pip) ──────────────────────────────────"

if phold --version &>/dev/null; then
    ok "phold already installed: $(phold --version 2>&1 | head -1)"
else
    info "Installing phold via pip..."
    pip install phold
    ok "phold installed: $(phold --version 2>&1 | head -1)"
fi

# ---------------------------------------------------------------------------
# Step 3b — Install pharokka via pip (guards against pre-existing env
#            that was created before pharokka was added to the conda block)
# ---------------------------------------------------------------------------
echo ""
echo "── Step 3b: pharokka (pip) ──────────────────────────────"

if pharokka --version &>/dev/null; then
    ok "pharokka already installed: $(pharokka --version 2>&1 | head -1)"
else
    info "Installing pharokka via pip..."
    pip install pharokka
    ok "pharokka installed: $(pharokka --version 2>&1 | head -1)"
fi

# ---------------------------------------------------------------------------
# Step 4 — Download phold database
# ---------------------------------------------------------------------------
echo ""
echo "── Step 4: phold database ───────────────────────────────"

if phold install --check 2>/dev/null; then
    ok "phold database already installed"
else
    info "Downloading phold database (~2 GB) ..."
    phold install
    ok "phold database installed"
fi

micromamba deactivate

# ---------------------------------------------------------------------------
# Step 5 — Create phynteny env (Phynteny Transformer)
# 2026-06: switched from classic phynteny (LSTM, py<3.11, 120-gene cap) to
# phynteny_transformer (py>=3.9, no gene cap, more accurate — Grigson 2025).
# Kept in a SEPARATE env (torch/transformers are heavy) that step 06 activates.
# ---------------------------------------------------------------------------
echo ""
echo "── Step 5: phynteny_transformer environment ──────────"

if [[ -d "${MAMBA_ROOT_PREFIX}/envs/${PHYNTENY_ENV}" ]]; then
    ok "Environment '${PHYNTENY_ENV}' already exists"
    warn "To reinstall: micromamba env remove -n ${PHYNTENY_ENV} && bash 00_install_env.sh"
else
    info "Creating '${PHYNTENY_ENV}' with python=3.10 ..."
    micromamba create -y -n "${PHYNTENY_ENV}" \
        -c conda-forge -c bioconda \
        python=3.10 \
        biopython \
        pandas
    ok "Environment '${PHYNTENY_ENV}' created."
fi

micromamba activate "${PHYNTENY_ENV}"

# setuptools provides pkg_resources (needed by install_models). torch>=2.12
# pins setuptools<82, so install a compatible setuptools (still has
# pkg_resources). We also drop classic `phynteny` first so its install_models
# entry-point and its scikit-learn<=1.2.2 pin don't linger.
pip uninstall -y phynteny 2>/dev/null || true
pip install "setuptools<82"

if python -c "import phynteny_transformer" 2>/dev/null || command -v phynteny_transformer >/dev/null 2>&1; then
    ok "phynteny_transformer already installed"
else
    info "Installing phynteny_transformer via pip..."
    pip install phynteny_transformer
    ok "phynteny_transformer installed"
fi

# CRITICAL compat fix: phynteny_transformer pulls numpy 2.x (via torch), but the
# scikit-learn left at 1.2.2 by classic phynteny was compiled against numpy 1.x
# -> "ValueError: numpy.dtype size changed (Expected 96, got 88)" at import.
# Upgrade scikit-learn to a numpy-2-compatible build so its C extensions match.
info "Aligning scikit-learn with numpy 2.x (binary-compat fix)..."
pip install -U "scikit-learn>=1.4"
python -c "import sklearn, numpy; print('  sklearn', sklearn.__version__, '/ numpy', numpy.__version__)" \
    || warn "sklearn/numpy still incompatible — try: pip install -U 'scikit-learn>=1.4'"

# Transformer ships NO model weights by default -> download once via install_models.
# (Fallback: download the Zenodo tarball
#  https://zenodo.org/records/15584824/files/phynteny_transformer_models_2025-06-03.tar.gz
#  and untar, then pass it to phynteny_transformer with -m <dir>.)
info "Downloading phynteny_transformer pre-trained models (install_models)..."
install_models || warn "install_models failed — install setuptools, or download models from Zenodo (2025-06-03) and untar"

micromamba deactivate

# ---------------------------------------------------------------------------
# Step 6 — Verification
# ---------------------------------------------------------------------------
echo ""
echo "── Step 6: Verification ─────────────────────────────────"

micromamba activate "${MAMBA_ENV}"

TOOLS_OK=true
check_tool() {
    local tool="$1" cmd="${2:-$1}"
    if command -v "$cmd" &>/dev/null; then
        local ver; ver=$($cmd --version 2>&1 | head -1 | tr -d '\n')
        ok "${tool}: ${ver}"
    else
        echo -e "${RED}[MISSING]${NC} ${tool}"
        TOOLS_OK=false
    fi
}

check_tool "pharokka"
check_tool "phold"
check_tool "foldseek"
check_tool "mmseqs2"  "mmseqs"
check_tool "diamond"
check_tool "python"

python -c "import Bio" 2>/dev/null && ok "biopython imported" \
    || { warn "biopython not importable — run: pip install biopython"; TOOLS_OK=false; }
python -c "import pandas" 2>/dev/null && ok "pandas imported" \
    || { warn "pandas not importable — run: pip install pandas"; TOOLS_OK=false; }

micromamba deactivate

# Check phynteny env
micromamba activate "${PHYNTENY_ENV}"
python -c "import phynteny_utils" 2>/dev/null && ok "phynteny env OK" \
    || { warn "phynteny not importable in '${PHYNTENY_ENV}' — check step 5"; TOOLS_OK=false; }
micromamba deactivate

echo ""
if [[ "${TOOLS_OK}" == "true" ]]; then
    echo -e "${GREEN}════ Setup complete ════${NC}"
    echo ""
    echo "  phagefactor env → pharokka, phold, foldseek (python 3.11)"
    echo "  phynteny env  → phynteny (python 3.10, separate due to version constraint)"
else
    echo -e "${YELLOW}════ Setup finished with warnings — check above ════${NC}"
fi
echo ""
echo "Next step: bash steps/00b_download_foldseek_dbs.sh"
echo "           (pdb100, afdb-swissprot, and afdb50 — see docs/databases.md)"
