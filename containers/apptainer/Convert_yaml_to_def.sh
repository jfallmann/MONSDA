# zsh -i Convert_yaml_to_def.sh 

for i in ${HOME}/MONSDA/envs/*.yaml;do a=${i%.yaml};b=${a##*/};cat Container_skeleton.def |sed -e "s,##environment.yaml##,${i} /opt/envs/,g" -e "s,##environment_short.yaml##,/opt/envs/${b}.yaml,g" -e "s/##ENV_NAME##/${b}/g" > ${b}.def;done
sed -i 's,=base,=monsda_base,g' base.def
conda activate apptainer
mkdir -p TMP_BUILDDIR
for i in *.def;do a=${i%.def};APPTAINER_TMPDIR="TMP_BUILDDIR"; apptainer build --force --fakeroot ${a}.sif ${i};done
rm -rf TMP_BUILDDIR

# for i in bbduk.sif bbmap.sif bwa2.sif bwameth.sif bwa.sif ciri2.sif cutadapt.sif deseq2_DE.sif dexseq_DEU.sif diego_DAS.sif dorado.sif drimseq_DTU.sif edger_DE.sif fastp.sif fastqc.sif guppy.sif hisat2.sif index.sif kallisto.sif macs.sif minimap.sif perl.sif picard.sif piranha.sif qc.sif salmon.sif samtools.sif scyphy.sif segemehl3.sif segemehl.sif sra.sif star.sif summary.sif test.sif trimgalore.sif trimmomatic.sif trnascan.sif ucsc.sif umitools.sif zip.sif ;do a=${i%.sif};apptainer push ${i} oras://docker.io/jfallmann/monsda:${a};done