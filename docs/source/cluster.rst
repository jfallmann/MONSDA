==============
RUN ON CLUSTER
==============

Snakemake
=========

SLURM
-----

You can either use the slurm profile adapted from
[Snakemake-Profiles](https://github.com/Snakemake-Profiles/slurm) that comes with this repository, or go
through the process of manually creating one, either using the cookiecutter example in the
``Snakemake-Profiles`` repository or on your own. To use the preconfigured example that comes with this
repository simply adapt the call below to your needs.

``
python nextsnakes/RunSnakemake.py -j ${cpus} --configfile ${config.json} --directory ${PWD} --profile nextsnakes/slurm --cluster-config nextsnakes/cluster/config_slurm.yamlx
``

Further adaptions like grouping of jobs and advanced configs for rule
based performance increase will be tackled in future releases.

SGE(outdated)
-------------

Define the cluster config file and for SGE support simply append

``
--cluster "qsub -q ${QUEUENAME} -e ${PWD}/sgeerror -o ${PWD}/sgeout -N ${JOBNAME}" --jobscript nextsnakes/cluster/sge.sh
``

Nextflow
========

Cluster config for Nextflow is currently not implemented and will follow in a future release
