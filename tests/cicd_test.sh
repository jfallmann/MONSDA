#!/usr/bin/env bash
set -euo pipefail

# Optional: activate pre-created conda env if available
if command -v conda >/dev/null 2>&1; then
	eval "$(conda shell.bash hook)"
	if conda env list | awk '{print $1}' | grep -qx monsda-test; then
		conda activate monsda-test
	fi
fi

# Get the directory of the current script
SCRIPT_DIR=$(dirname "$(realpath "$0")")
echo ${SCRIPT_DIR}
cd ${SCRIPT_DIR}

# Keep local execution aligned with CI
python -m pip install --upgrade pip
python -m pip install -e ..
python -m pip install pytest
pytest -q test_Utils.py

# Optional integration smoke test (opt-in)
if [[ "${RUN_INTEGRATION_TESTS:-0}" == "1" ]]; then
	mkdir -p CONDALIB
	monsda -j 6 -c config_Test.json --directory ${PWD} --use-conda --conda-prefix CONDALIB --save
fi