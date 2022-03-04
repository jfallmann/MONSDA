==============
RUN ON SLURM
==============

Snakemake
=========

You can either use the slurm profile adapted from
[Snakemake-Profiles](https://github.com/Snakemake-Profiles/slurm) that comes with this repository, or go
through the process of manually creating one, either using the cookiecutter example in the
``Snakemake-Profiles`` repository or on your own. You can also adapt the example that comes with this
repository and execute

```
monsda -j ${cpus} --configfile ${config.json} --directory ${PWD} --profile ${path_to_slurm_profile}
```

Further adaptions like grouping of jobs and advanced configs for rule
based performance increase will be tackled in future releases.

Nextflow
========

Cluster config for Nextflow follows the description [Nextflow-Executors](https://www.nextflow.io/docs/latest/executor.html) and [Nextflow-Profiles](https://www.nextflow.io/docs/latest/config.html#config-profiles). To use ```slurm``` as executor you can adapt the profile that comes with this repository and simply append 
```export NXF_EXECUTOR=slurm```
to the call to ```MONSDA```.
