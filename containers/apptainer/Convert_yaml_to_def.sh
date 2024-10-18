for i in ${HOME}/MONSDA/envs/*.yaml;do a=${i%.yaml};b=${a##*/};cat Container_skeleton.def |sed -e "s,##environment.yaml##,${i} /opt/envs/,g" -e "s,##environment_short.yaml##,/opt/envs/${b}.yaml,g" -e "s/##ENV_NAME##/${b}/g" > ${b}.def;done
sed -i 's,=base,=monsda_base,g' base.def
conda activate apptainer
mkdir -p TMP_BUILDDIR
for i in *.def;do a=${i%.def};APPTAINER_TMPDIR="TMP_BUILDDIR"; apptainer build --fakeroot ${a}.sif ${i};done
rm -rf TMP_BUILDDIR