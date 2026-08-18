#!/usr/bin/env bash
set -uo pipefail

sif=${1:?Usage: smoke_test.sh IMAGE.sif [LOG]}
log=${2:-/dev/stderr}

if command -v apptainer >/dev/null 2>&1; then
    runner=(apptainer)
else
    runner=(conda run -n apptainer apptainer)
fi

"${runner[@]}" exec --cleanenv "$sif" sh -c '
    report=$(mktemp)
    test -x /opt/MONSDA/scripts/Analysis/build_count_table_simple.py || exit 1
    found=0
    for bindir in /opt/conda/envs/*/bin; do
        test -d "$bindir" || continue
        found=1
        find "$bindir" -type f -perm /111 -exec sh -c '\''
            for executable do
                magic=$(od -An -tx1 -N4 "$executable" 2>/dev/null | tr -d " \n")
                test "$magic" = 7f454c46 || continue
                output=$(ldd "$executable" 2>&1 || true)
                if printf "%s\n" "$output" | grep -Eq "version .+ not found"; then
                    printf "%s\n%s\n" "$executable" "$output"
                fi
            done
        '\'' sh {} + >> "$report"
    done
    if test "$found" -eq 0 || test -s "$report"; then
        cat "$report" >&2
        rm -f "$report"
        exit 1
    fi
    rm -f "$report"
' >>"$log" 2>&1
