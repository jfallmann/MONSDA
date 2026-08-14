#!/usr/bin/env bash
set -uo pipefail
shopt -s nullglob

repo=${HOME}/Gits/Monsda
version=$(cd "$repo" && python scripts/container_matrix.py --version | cut -d= -f2)
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
: > push-ghcr-results.tsv

for sif in "$script_dir"/*-${version}.sif; do
    printf 'Preparing %s\n' "$sif"
    base_name=${sif##*/}
    environment=${base_name%-${version}.sif}
    target="oras://ghcr.io/jfallmann/monsda:${environment}-${version}"

    if ! "$script_dir/smoke_test.sh" "$sif" push-ghcr.log; then
        printf '%s\tsmoke-failed\t%s\n' "$environment" "$target" >> push-ghcr-results.tsv
        continue
    fi
    printf 'Pushing %s\n' "$target"
    if conda run -n apptainer apptainer push --allow-unsigned "$sif" "$target"; then
        printf '%s\tpassed\t%s\n' "$environment" "$target" \
            >> push-ghcr-results.tsv
    else
        printf '%s\tfailed\t%s\n' "$environment" "$target" \
            >> push-ghcr-results.tsv
    fi
done 2>&1 | tee push-ghcr.log
