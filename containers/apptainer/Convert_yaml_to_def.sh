#!/usr/bin/env bash
set -uo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

for i in ${HOME}/MONSDA/envs/*.yaml; do
    a=${i%.yaml}
    b=${a##*/}
    sed -e "s,##environment.yaml##,${i} /opt/envs/,g" -e "s,##environment_short.yaml##,/opt/envs/${b}.yaml,g" -e "s/##ENV_NAME##/${b}/g" "$script_dir/Container_skeleton.def" > "$script_dir/${b}.def"
done
sed -i 's,=base,=monsda_base,g' "$script_dir/base.def"
mkdir -p "$script_dir/TMP_BUILDDIR"
rm -f "$script_dir/build.log"
for i in "$script_dir"/*.def; do
    a=${i%.def}
    sif=${a}.sif
    if APPTAINER_TMPDIR="$script_dir/TMP_BUILDDIR" conda run -n apptainer apptainer build --force --fakeroot "$sif" "$i" >> "$script_dir/build.log" 2>&1 && "$script_dir/smoke_test.sh" "$sif" "$script_dir/build.log"; then
        printf '%s\tsmoke-passed\n' "${a##*/}" | tee -a "$script_dir/build.log"
    else
        rm -f "$sif"
        printf '%s\tsmoke-failed\n' "${a##*/}" | tee -a "$script_dir/build.log"
    fi
done
rm -rf "$script_dir/TMP_BUILDDIR"
