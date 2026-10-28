#!/usr/bin/env bash
# Snakemake post-deploy hook for envs/rustar.yaml.
#
# rustar-aligner is not available on bioconda or conda-forge, so the conda env
# only provides the Rust toolchain and the binary is installed from crates.io
# here. Snakemake runs this automatically after creating the environment.
set -euo pipefail

RUSTAR_VERSION="${RUSTAR_VERSION:-0.1.0}"
PREFIX="${CONDA_PREFIX:?CONDA_PREFIX is not set}"

if [ -x "$PREFIX/bin/rustar-aligner" ]; then
    echo "rustar-aligner already present in $PREFIX/bin"
else
    export CARGO_HOME="$PREFIX/share/cargo"
    cargo install rustar-aligner --version "$RUSTAR_VERSION" --locked \
        --root "$PREFIX"
fi

# MONSDA calls the mapper as 'rustar'
ln -sf "$PREFIX/bin/rustar-aligner" "$PREFIX/bin/rustar"
"$PREFIX/bin/rustar" --version
