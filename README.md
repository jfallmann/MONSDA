# snakes

[![Documentation Status](https://readthedocs.org/projects/MONSDA/badge/?version=latest)](https://MONSDA.readthedocs.io/en/latest/?badge=latest)


Welcome to MONSDA, a modular assembler of snakemake and nexflow workflows
for HTS analysis from sra download, preprocessing and mapping to
postprocessing/analysis and ucsc track generation centered on a single config file

Simply clone with ```git clone```.

For details on ```snakemake``` and it's features please refer to the [snakemake documentation](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html).

In general it is necessary to write a configuration file containing information on paths, files to process and settings beyond default for mapping tools and others.
The template on which analysis is based can be found in the ```config``` directory.

For ```snakemake``` to be fully FAIR, one needs to use ```conda``` or similar environment management systems. For details on ```conda``` please refer to the [conda manual](https://docs.conda.io/en/latest/).

This workflow collection makes heavy use of ```conda``` and especially the [bioconda](https://bioconda.github.io) channel.

To create a working environment for this repository please install the snakemake environment as found in the ```envs``` directory like so:

```
conda env create -n snakemake -f snakes/envs/snakemake.yaml
```

The ```envs``` directory holds all the environments needed to run the pipelines in the ```workflows``` directory, these will be installed automatically when needed.

More information can be found in the official [documentation](https://MONSDA.readthedocs.io/en/latest/?badge=latest)


## What is happening

This repository hosts the executable ```MONSDA.py``` which acts a wrapper around ```snakemake``` and the ```config.json``` file.
The ```config.json``` holds all the information that is needed to run the jobs and will be parsed by ```MONSDA.py``` and split into sub-configs that can later be found in the directory ```SubSnakes```.

To successfully run an analysis pipeline, a few steps have to be followed:
  * Clone this repository to you working directory, or symlink it there. *DO NOT CHANGE THE NAME OF THE REPO!!!* This is needed for subworkflow generation.
  * Directory structure: The structure for the directories is dictated by the ICS (IdentifierConditionSetting relationship) in the config file
  * Config file: This is the central part of the analysis. Depending on this file ```MONSDA.py``` will determine processing steps and generate according config and ```snakemake``` workflow files to run each subworkflow until all processing steps are done.

## Run the pipeline
Simply run

```
python snakes/MONSDA.py
```

to see the help and available options that will be passed through to ```snakemake```.

To start a simple run call
```
python snakes/MONSDA.py -j NUMBER_OF_CORES --configfile YOUR_CONFIG.json --directory ${PWD}
```
or add additional arguments for ```snakemake``` as you see fit.

### Run on cluster

####SLURM

You can either use the slurm profile adapted from [Snakemake-Profiles](https://github.com/Snakemake-Profiles/slurm) that comes with this repository, or go through the process of manually creating one, either using the cookiecutter example in the ```Snakemake-Profiles``` repository or on your own. To use the preconfigured example that comes with this repository simply adapt the call below to your needs.

```python snakes/MONSDA.py -j ${cpus} --configfile ${config.json} --directory ${PWD} --profile snakes/slurm --cluster-config snakes/cluster/config_slurm.yamlx```

Further adaptions like grouping of jobs and advanced configs for rule based performance increase will follow.

####SGE(outdated)

Define the cluster config file and for SGE support simply append ```--cluster "qsub -q ${QUEUENAME} -e ${PWD}/sgeerror -o ${PWD}/sgeout -N ${JOBNAME}" --jobscript snakes/cluster/sge.sh```

## Contribute
If you like this project, are missing features, want to contribute or
file bugs please leave an issue or contact me directly.

To contribute new tools feel free to adopt existing ones,
there should be a number of examples available that cover
implementation details for almost all sorts of tools. If you need to
add new python/groovy functions for processing of options or
parameters add them to the corresponding file in the lib directory.
New environments go into the envs directory, new subworkflows into the
workflows directory. Do not forget to also extend the template.json
and add some documentation.

PRs always welcome.

