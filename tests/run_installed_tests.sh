#!/usr/bin/env bash
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="${MONSDA_TEST_HOME:-${TMPDIR:-/tmp}/monsda_installed_test}"
ENVDIR="${MONSDA_TEST_ENV:-$WORK/env}"
ENVFILE="${MONSDA_TEST_ENVFILE:-$REPO/envs/monsda_min.yaml}"
CONDA="${CONDA_EXE:-conda}"
DIST="$WORK/dist"
STAGE="$WORK/stage"

if [ "${1:-}" = "--help" ]; then
    echo "usage: $0 [pytest args]"
    echo "builds a wheel from the checkout, installs it into a conda env made from envs/monsda_min.yaml"
    echo "and runs the tests against the installed package (not the working tree)"
    echo "env: MONSDA_TEST_HOME work dir, MONSDA_TEST_ENV existing env prefix, MONSDA_TEST_ENVFILE env yaml,"
    echo "     MONSDA_TEST_FRESH=1 recreate the env"
    exit 0
fi

mkdir -p "$WORK"

if [ "${MONSDA_TEST_FRESH:-0}" = "1" ] && [ -d "$ENVDIR" ]; then
    rm -rf "$ENVDIR"
fi
if [ ! -x "$ENVDIR/bin/python" ]; then
    "$CONDA" env create -y -p "$ENVDIR" -f "$ENVFILE"
fi
PY="$ENVDIR/bin/python"

rm -rf "$DIST"
"$PY" -m pip install -q --upgrade build
(cd "$REPO" && "$PY" -m build --wheel --outdir "$DIST" .)
WHEEL="$(ls "$DIST"/monsda-*.whl)"

"$PY" -m pip uninstall -y -q monsda 2>/dev/null || true
"$PY" -m pip install -q --no-deps --force-reinstall "$WHEEL"
"$PY" -m pip install -q pytest httpx

SHARE="$("$PY" -c 'import sys; print(sys.prefix)')/share/MONSDA"

rm -rf "$STAGE"
mkdir -p "$STAGE"
cp -r "$REPO/tests" "$STAGE/tests"
rm -rf "$STAGE/tests/LOGS" "$STAGE/tests/tool_audit/__pycache__"
for d in workflows envs scripts configs; do
    ln -s "$SHARE/$d" "$STAGE/$d"
done
for d in containers docker .github; do
    [ -e "$REPO/$d" ] && ln -s "$REPO/$d" "$STAGE/$d"
done

cd "$STAGE/tests"
env -u PYTHONPATH PYTHONNOUSERSITE=1 "$PY" - <<'EOF'
import MONSDA
assert "site-packages" in MONSDA.__file__, MONSDA.__file__
print("testing installed MONSDA", MONSDA.__version__, "from", MONSDA.__file__)
EOF
export PATH="$ENVDIR/bin:$PATH"
env -u PYTHONPATH PYTHONNOUSERSITE=1 "$PY" -m pytest -p no:cacheprovider "$@"
