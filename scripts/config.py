# phageFACTor config shim
# ========================
# All analysis scripts in scripts/ do `from config import ...`. This shim loads
# the single source of truth -- config/config.yaml via config/config.py (relative
# paths, env-overridable) -- and re-exports its names, so the scripts need ZERO
# changes — the single source of truth is config/config.yaml (no hardcoded paths).
import importlib.util as _ilu
from pathlib import Path as _Path

_cfg = _Path(__file__).resolve().parents[1] / "config" / "config.py"
_spec = _ilu.spec_from_file_location("_phagefactor_config", _cfg)
_mod = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
globals().update({_k: _v for _k, _v in vars(_mod).items() if not _k.startswith("__")})
