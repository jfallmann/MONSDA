#!/usr/bin/env bash
set -uo pipefail

usage() {
    printf 'Usage: %s [-e environment[,environment...]] [repo] [out]\n' "$(basename "${BASH_SOURCE[0]}")" >&2
    printf '  -e  Build only the given environment(s) instead of the full matrix.\n' >&2
    printf '      May be repeated or given as a comma-separated list, e.g.\n' >&2
    printf '      -e deseq2_DE -e edger_DE  or  -e deseq2_DE,edger_DE\n' >&2
}

only_environments=()
while getopts ":e:h" opt; do
    case "$opt" in
        e) IFS=',' read -ra envs <<< "$OPTARG"; only_environments+=("${envs[@]}") ;;
        h) usage; exit 0 ;;
        \?) printf 'Unknown option: -%s\n' "$OPTARG" >&2; usage; exit 1 ;;
        :) printf 'Option -%s requires an argument\n' "$OPTARG" >&2; usage; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

repo=${1:-${HOME}/Gits/Monsda}
out=${2:-${HOME}/Work2/Container/MONSDA}
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
version=$(cd "$repo" && python scripts/container_matrix.py --version | cut -d= -f2)
revision=$(cd "$repo" && git rev-parse HEAD)
if (( ${#only_environments[@]} > 0 )); then
    mapfile -t all_environments < <(cd "$repo" && python scripts/container_matrix.py | python -c 'import json,sys; print("\n".join(json.load(sys.stdin)["environment"]))')
    environments=()
    for environment in "${only_environments[@]}"; do
        [[ -z "$environment" ]] && continue
        if [[ ! " ${all_environments[*]} " == *" $environment "* ]]; then
            printf 'Unknown environment: %s (available: %s)\n' "$environment" "${all_environments[*]}" >&2
            exit 1
        fi
        environments+=("$environment")
    done
else
    mapfile -t environments < <(cd "$repo" && python scripts/container_matrix.py | python -c 'import json,sys; print("\n".join(json.load(sys.stdin)["environment"]))')
fi
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
