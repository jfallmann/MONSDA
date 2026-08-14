#!/usr/bin/env bash
set -uo pipefail

repo=${HOME}/Gits/Monsda
out=${HOME}/Work/Container/MONSDA
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
version=$(cd "$repo" && python scripts/container_matrix.py --version | cut -d= -f2)
revision=$(cd "$repo" && git rev-parse HEAD)
mapfile -t environments < <(cd "$repo" && python scripts/container_matrix.py | python -c 'import json,sys; print("\n".join(json.load(sys.stdin)["environment"]))')
mkdir -p "$out/logs" "$out/tmp"
printf 'environment\tstatus\timage\n' > "$out/results.tsv"

for environment in "${environments[@]}"; do
    sif="$out/${environment}-${version}.sif"
    log="$out/logs/${environment}-${version}.log"
    image="local/monsda:${environment}-${version}"
    available=$(df --output=avail -B1 "$out" | tail -1)
    if (( available < 16106127360 )); then
        printf '%s\t%s\t%s\n' "$environment" "stopped-low-disk" "$sif" | tee -a "$out/results.tsv"
        break
    fi
    if [[ -s "$sif" ]] && "$script_dir/smoke_test.sh" "$sif" "$log"; then
        printf '%s\t%s\t%s\n' "$environment" "existing-smoke-passed" "$sif" | tee -a "$out/results.tsv"
        continue
    fi
    rm -f "$sif"
    printf '%s\t%s\n' "$environment" "building" | tee -a "$out/build.log"
    if ! docker build --pull --build-arg "ENVIRONMENT=$environment" --build-arg "MONSDA_VERSION=$version" --build-arg "VCS_REF=$revision" -t "$image" -f "$repo/docker/Dockerfile" "$repo" >>"$log" 2>&1; then
        printf '%s\t%s\t%s\n' "$environment" "docker-build-failed" "$sif" | tee -a "$out/results.tsv"
        continue
    fi
    if APPTAINER_TMPDIR="$out/tmp" conda run -n apptainer apptainer build --force --fakeroot "$sif" "docker-daemon://$image" >>"$log" 2>&1 && "$script_dir/smoke_test.sh" "$sif" "$log"; then
        printf '%s\t%s\t%s\n' "$environment" "smoke-passed" "$sif" | tee -a "$out/results.tsv"
    else
        rm -f "$sif"
        printf '%s\t%s\t%s\n' "$environment" "smoke-failed" "$sif" | tee -a "$out/results.tsv"
    fi
    docker image rm "$image" >>"$log" 2>&1 || true
done

rm -rf "$out/tmp"
printf 'finished\n' | tee -a "$out/build.log"
