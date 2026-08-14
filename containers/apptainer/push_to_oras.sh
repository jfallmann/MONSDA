#!/usr/bin/env bash
set -uo pipefail
shopt -s nullglob

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=${HOME}/Gits/Monsda

if [[ ! -f "$repo_root/scripts/container_matrix.py" ]]; then
    printf 'Could not locate MONSDA repo root from %s\n' "$script_dir" >&2
    exit 1
fi

version=$(python "$repo_root/scripts/container_matrix.py" --version | cut -d= -f2)
search_dir=${1:-$PWD}
log_file="$search_dir/push-ghcr.log"
results_file="$search_dir/push-ghcr-results.tsv"
: > "$results_file"

matches=("$search_dir"/*-${version}.sif)
if (( ${#matches[@]} == 0 )); then
    printf 'No SIF files matching %s/*-%s.sif found in %s\n' "$search_dir" "$version" "$search_dir" | tee "$log_file" >&2
    exit 1
fi

for sif in "${matches[@]}"; do
    printf 'Preparing %s\n' "$sif"
    base_name=${sif##*/}
    environment=${base_name%-${version}.sif}
    target="oras://ghcr.io/jfallmann/monsda:${environment}-${version}"

    if ! "$script_dir/smoke_test.sh" "$sif" "$log_file"; then
        printf '%s\tsmoke-failed\t%s\n' "$environment" "$target" >> "$results_file"
        continue
    fi
    printf 'Pushing %s\n' "$target"
    if conda run -n apptainer apptainer push --allow-unsigned "$sif" "$target"; then
        printf '%s\tpassed\t%s\n' "$environment" "$target" >> "$results_file"
    else
        printf '%s\tfailed\t%s\n' "$environment" "$target" >> "$results_file"
    fi
done 2>&1 | tee "$log_file"
