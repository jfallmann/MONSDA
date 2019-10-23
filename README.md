# snakes

Collection of snakemake workflows for NGS analysis from mapping to featurecount and track generation
Contains sumodule so clone with ```git clone --recursive``` or if already cloned with pull submodules with 
```
git submodule init
git submodule update 
```

See [stackoverflow](https://stackoverflow.com/questions/25200231/cloning-a-git-repo-with-all-submodules) and 
[git docs](https://git-scm.com/book/en/v2/Git-Tools-Submodules) for details.

For details on ```snakemake``` and it's features please refer to the [snakemake documentation](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html).

In general it is necessary to write a configuration file containing information on paths, files to process and settings beyond default for mapping tools and others.
For examples on this please have a look into the ```config``` directory.

For ```snakemake``` to be fully FAIR, one needs to use ```conda``` or similar environment management systems. For details on ```conda``` please refer to the [conda manual](https://docs.conda.io/en/latest/).

This workflow collection makes heavy use of ```conda``` and especially the [bioconda](https://bioconda.github.io) channel.

For distribution of jobs one can either rely on local hardware, use scheduling software like the [SGE](https://docs.oracle.com/cd/E19957-01/820-0699/chp1-1/index.html) or [Slurm](https://slurm.schedmd.com/documentation.html).

This manual will only show examples on local and SGE usage, but more information on how to use other scheduling software is available elsewhere.

## What is happening

This repository hosts the executable ```RunSnakemake.py``` which acts a wrapper around ```snakemake``` and the ```config.json``` file. 
The ```config.json``` holds all the information that is needed to run the jobs and will be parsed by ```RunSnakemake.py``` and split into sub-configs that can later be found in the directory ```SubSnakes```.

To successfully run an analysis pipeline, a few steps have to be followed:
  * Clone this repository to you working directory, or symlink it there. *DO NOT CHANGE THE NAME OF THE REPO!!!* This is needed for subworkflow generation.
  * Directory structure: The structure for the directories is dictated by the config file
  * Config file: This is the central part of the analysis. Depending on this file ```RunSnakemake.py``` will determine processing steps and generate according config and ```snakemake``` workflow files to run each subworkflow until all processing steps are done.

## Create config.json

The configuration file follows standard ```json``` format and contains all the information needed to run the analysis.
It consists of multiple sections which will be explained in detail.
For examples please refer to the ```config``` directory of the repo.

### Cluster config
This is separate from the main configuration, for details please follow the explanations in the [snakemake documentation](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html).
For an example have a look at the ```cluster``` directory.

## Run the pipeline
Simply run

```
python snakes/RunSnakemake.py
```

to see the help and available options that will be passed through to ```snakemake```.

To start a simple run call
```
python snakes/RunSnakemake.py -j NUMBER_OF_CORES --configfile YOUR_CONFIG.json --directory ${PWD}
```
or add additional arguments for ```snakemake``` as you see fit.

[//]: # (#### Cluster SGE)
[//]: #(```snakemake --cluster "qsub -q ${QUEUENAME} -e ${PWD}/sgeerror -o ${PWD}/sgeout -N ${JOBNAME}" --jobscript Workflow/cluster/sge.sh -j ${THREADS} --configfile config.json --snakefile workflows/mapping.smk --use-conda --directory ${PWD} --printshellcmds &> LOGS/mapping.log```)
