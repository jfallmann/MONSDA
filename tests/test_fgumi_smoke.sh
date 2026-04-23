#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(dirname "$(realpath "$0")")
cd "${SCRIPT_DIR}"

# Mirror existing integration style: expose tests/data as working inputs
ln -fs data/* .

# Generate tiny gzipped UMI FASTQ inputs
python data/make_fgumi_umi_fixture.py

# Keep version in sync with installed MONSDA
VERSION=$(monsda --version 2>&1 | sed 's/MONSDA version //g')
sed -i "s/\"VERSION\": \"FIXME\"/\"VERSION\": \"${VERSION}\"/g" config_fgumi_test.json

mkdir -p CONDALIB

# --save keeps this as a lightweight workflow generation smoke test
monsda -j 2 -c config_fgumi_test.json --directory "${PWD}" --use-conda --conda-prefix CONDALIB --save

echo "FGUMI smoke test workflow generation completed."
