#!/usr/bin/env bash -el
conda activate monsda-test

# Get the directory of the current script
SCRIPT_DIR=$(dirname "$(realpath "$0")")
echo ${SCRIPT_DIR}
cd ${SCRIPT_DIR}
mkdir -p CONDALIB
# Run the MONSDA command with the updated -c parameter
monsda -j 6 -c config_Test.json --directory ${PWD} --use-conda --conda-prefix CONDALIB --save && echo "MONSDA test passed" || echo "MONSDA test failed"