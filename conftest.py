# conftest.py -- pytest session setup; loads integration config before tests.
# resolution order (setdefault: first write wins; shell exports beat all):
#   1. shell export: SPRINX_CANONICAL_CM=... pytest
#   2. .env at project root (committed, gitignored only if it contains secrets;
#      these are paths not secrets so committing is fine -- see .env.example)
# absolute paths required; relative paths fail silently when cwd differs.
import os

def _load_env(path):
    if not os.path.exists(path):
        return
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                k, _, v = line.partition("=")
                os.environ.setdefault(k.strip(), v.strip())

_load_env(os.path.join(os.path.dirname(__file__), ".env"))

for _var in ("SPRINX_CANONICAL_CM", "SPRINX_ARMLESS_CM_DIR"):
    _v = os.environ.get(_var)
    if _v and not os.path.isabs(_v):
        raise ValueError(f"{_var}={_v!r}: must be an absolute path")
