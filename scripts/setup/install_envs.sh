#!/usr/bin/env bash
# =============================================================================
# install_envs.sh — create the conda environments phageFACTor runs its tools in.
# Called by `phagefactor install-envs [--tools] [--phynteny]`.
#   --tools     envs/phagefactor.yaml (pharokka, phold, foldseek); not needed
#               when phageFACTor itself was installed with conda
#   --phynteny  envs/phynteny.yaml + the phynteny_transformer models
# No flag: each env that is missing. Uses micromamba, else mamba, else conda.
# =============================================================================
set -euo pipefail
CODE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

log() { echo "[$(date '+%H:%M:%S')] $*"; }
die() { echo "[ERROR] $*" >&2; exit 1; }

if   [[ -n "${MAMBA_EXE:-}" && -x "${MAMBA_EXE}" ]]; then M="${MAMBA_EXE}"
elif command -v micromamba &>/dev/null; then M="$(command -v micromamba)"
elif command -v mamba &>/dev/null; then M="$(command -v mamba)"
elif command -v conda &>/dev/null; then M="$(command -v conda)"
else die "no micromamba, mamba or conda found: install one first (https://mamba.readthedocs.io)"; fi

env_exists() { "$M" env list 2>/dev/null | awk '{print $1}' | grep -qx "$1"; }
run_in() { local e="$1"; shift
    if [[ "$(basename "$M")" == conda ]]; then "$M" run --no-capture-output -n "$e" "$@"; else "$M" run -n "$e" "$@"; fi; }
create() {  # create <name> <yaml>
    if env_exists "$1"; then log "env '$1' already exists (remove it first to rebuild)"; return 1; fi
    log "creating env '$1' from $2 (10-20 min)"
    if [[ "$(basename "$M")" == micromamba ]]; then "$M" create -y -n "$1" -f "$2"
    else "$M" env create -n "$1" -f "$2"; fi
}

DO_TOOLS=0; DO_PHY=0
for a in "$@"; do
    case "$a" in
        --tools) DO_TOOLS=1 ;;
        --phynteny) DO_PHY=1 ;;
        *) die "unknown flag: $a (valid: --tools --phynteny)" ;;
    esac
done
if (( DO_TOOLS + DO_PHY == 0 )); then
    command -v phold &>/dev/null && command -v pharokka.py &>/dev/null || DO_TOOLS=1
    command -v phynteny_transformer &>/dev/null || DO_PHY=1
fi

if (( DO_TOOLS )); then
    create phagefactor "${CODE}/envs/phagefactor.yaml" || true
    run_in phagefactor phold --version && run_in phagefactor pharokka.py --version | head -1 \
        && run_in phagefactor foldseek version
fi
if (( DO_PHY )); then
    create phynteny "${CODE}/envs/phynteny.yaml" || true
    log "downloading the phynteny_transformer models"
    run_in phynteny install_models \
        || die "install_models failed: download the models from Zenodo (phynteny_transformer_models_2025-06-03) instead"
fi
log "Done. Next: phagefactor databases <dir> --all"
