RUN ON CLUSTER
==============

##SLURM

You can either use the slurm profile adapted from [Snakemake-Profiles](https://github.com/Snakemake-Profiles/slurm) that comes with this repository, or go through the process of manually creating one, either using the cookiecutter example in the ```Snakemake-Profiles``` repository or on your own. To use the preconfigured example that comes with this repository simply adapt the call below to your needs.

```python snakes/RunSnakemake.py -j ${cpus} --configfile ${config.json} --directory ${PWD} --profile snakes/slurm --cluster-config snakes/cluster/config_slurm.yamlx```

Further adaptions like grouping of jobs and advanced configs for rule based performance increase will follow.

##SGE(outdated)

Define the cluster config file and for SGE support simply append ```--cluster "qsub -q ${QUEUENAME} -e ${PWD}/sgeerror -o ${PWD}/sgeout -N ${JOBNAME}" --jobscript snakes/cluster/sge.sh```
