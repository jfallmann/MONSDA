# MONSDA

[![Documentation Status](https://readthedocs.org/projects/monsda/badge/?version=latest)](https://monsda.readthedocs.io/en/latest/?badge=latest)


Welcome to MONSDA, Modular Organizer of Nextflow and Snakemake driven hts Data Analysis

Automating HTS analysis from data download, preprocessing and mapping to
postprocessing/analysis and track generation centered on a single config file.
MONSDA can create ```Snakemake``` and ```Nextflow``` workflows based on user defined configuration.
These workflows can either be saved to disk for manual inspection and execution or automatically executed.

For details on ```Snakemake``` and ```Nextflow``` and their features please refer to the corresponding [Snakemake](https://Snakemake.readthedocs.io/en/stable/tutorial/tutorial.html)  or [Nextflow](https://www.Nextflow.io/docs/latest/index.html) documentation.

In general it is necessary to write a configuration file containing information on paths, files to process and settings beyond default for mapping tools and others.
The template on which analysis is based can be found in the ```config``` directory.

For MONSDA to be as FAIR as possible, one needs to use ```conda``` or the faster drop-in replacement ```mamba```. For details on either please refer to the corresponding [conda](https://docs.conda.io/en/latest/) or [mamba](https://mamba.readthedocs.io/en/latest/) manual.

This workflow collection makes heavy use of ```conda``` and especially the [bioconda](https://bioconda.github.io) channel.

## Install MONSDA via ```conda``` or ```pip```

To install via ```conda/mamba``` simply run

```
mamba install -c bioconda -c conda-forge monsda
```

To install via ```pip``` you first need to create the ```MONSDA``` environment as found in the ```envs``` directory of this repository like so:

```
mamba env create -n monsda -f MONSDA/envs/monsda.yaml
```

The ```envs``` directory holds all the environments needed to run the pipelines in the ```workflows``` directory, these will be installed automatically alongside ```MONSDA```.

For that activate the ```monsda``` environment and run ```pip```

```
conda activate monsda
pip install MONSDA
```

More information can be found in the official [documentation](https://monsda.readthedocs.io/en/latest/?badge=latest)


## How does it work

This repository hosts the executable ```MONSDA.py``` which acts a wrapper around ```Snakemake``` and the ```config.json``` file.
The ```config.json``` holds all the information that is needed to run the jobs and will be parsed by ```MONSDA.py``` and split into sub-configs that can later be found in the directory ```SubSnakes``` or ```SubFlows``` respectively.

To successfully run an analysis pipeline, a few steps have to be followed:
  * Directory structure: The structure for the directories is dictated by the condition-tree in the config file
  * Config file: This is the central part of the analysis. Depending on this file ```MONSDA.py``` will determine processing steps and generate according config and ```Snakemake/Nextflow``` workflow files to run each subworkflow until all processing steps are done.

## Run the pipeline

Run

```
monsda
```
to see the help and available options that will be passed through to ```Snakemake``` or ```Nextflow```.

and 

```
monsda_configure
```

To spin up the configurator that guides you through the creation of config.json files.

Once a config.json is available you can start a ```Snakemake``` run with

```
monsda -j ${THREADS} --configfile ${CONFIG}.json --directory ${PWD} --conda-frontend mamba --conda-prefix ${PATH_TO_conda_envs}
```
and add additional arguments for ```Snakemake``` as you see fit.


For a ```Nextflow``` run use
```
monsda --Nextflow -j ${THREADS} --configfile ${CONFIG}.json --directory ${PWD}
```
and add additional arguments for ```Nextflow``` as you see fit.


### Run on workload manager

####SLURM

You can either use the slurm profile adapted from [Snakemake-Profiles](https://github.com/Snakemake-Profiles/slurm) that can be found in the ```profile_Snakemake``` directory, or go through the process of manually creating one, either using the cookiecutter example in the ```Snakemake-Profiles``` repository or on your own. 
For ```Nextflow``` a minimalist's example profile can be found under ```profile_Nextflow```.

Then run
```
monsda -j ${THREADS} --configfile ${CONFIG}.json --directory ${PWD} --conda-frontend mamba --profile ${SLURMPROFILE} --conda-prefix ${PATH_TO_conda_envs}
```
or
```
export NXF_EXECUTOR=slurm; monsda --Nextflow -j ${THREADS} --configfile ${CONFIG}.json --directory ${PWD}
```
respectively.


For other workload managers please refer to the documentation of ```Snakemake``` and ```Nextflow```.


## Contribute
If you like this project, are missing features, want to contribute or
file bugs please leave an issue or contact me directly.

To contribute new tools feel free to adopt existing ones,
there should be a number of examples available that cover
implementation details for almost all sorts of tools. If you need to
add new python/groovy functions for processing of options or
parameters add them to the corresponding file in the ```MONSDA``` directory.
New environments go into the envs directory, new subworkflows into the
workflows directory. Do not forget to also extend the template.json
and add some documentation.

PRs always welcome.

