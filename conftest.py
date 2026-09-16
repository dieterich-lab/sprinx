# integration tests skip if SPRINX_CANONICAL_CM /SPRINX_ARMLESS_CM_DIR are unset

# do `cp env.example .env` and set both to absolute paths
# before running integratoin tests

# OR do shell export of the variables before running pytest
import os

def _load_env(path):
    if not os.path.exists(path):
        return
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                k, _, v = line.partition("=")
                # setdefault: a variable keeps its value if already exported
                os.environ.setdefault(k.strip(), v.strip())

_load_env(os.path.join(os.path.dirname(__file__), ".env"))

for _var in ("SPRINX_CANONICAL_CM", "SPRINX_ARMLESS_CM_DIR"):
    _v = os.environ.get(_var)
    if _v and not os.path.isabs(_v):
        raise ValueError(f"{_var}={_v!r}: must be an absolute path")
