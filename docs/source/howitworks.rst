================
HOW DOES IT WORK
================

For details on **Snakemake** and it's features please refer to the
[snakemake
documentation](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html).
For details on **Nextflow** and it's features please refer to the
[nextflow
documentation](https://www.nextflow.io/docs/latest/index.html).

In general it is necessary to write a configuration file containing
information on paths, files to process and settings beyond default for
mapping tools and others.  The template on which analysis is based can
be found in the **config** directory.


**MONSDA** acts a wrapper around **Snakemake** or **Nextflow** and the corresponding**config.json** file.  The **config.json** holds all the information that is needed to run the jobs and will be parsed by **MONSDA** and split into sub-configs that can later be found in the directory **SubSnakes** or **SubFlows** respectively. Command line execution calls are storen in the directory '''JOBS''', so users can manipulate those or rerun them manually as needed.

To successfully run an analysis pipeline, a few steps have to be followed:
  * Install MONSDA either via **bioconda** or **pip**` folling the instruction in INSTALL
  * Directory structure: The structure for the directories is dictated by the ICS (Identifier-Condition-Setting relationship) in the config file
  * Config file: This is the central part of the analysis. Depending on this file **MONSDA** will determine processing steps and generate according config and workflow files to run each subworkflow until all processing steps are done.


In general it is necessary to write a configuration file containing
information on paths, files to process and settings beyond default for
mapping tools and others.  The template on which analysis is based can
be found in the **config** directory.

For **MONSDA** to be as FAIR as possible, we rely on
**conda** as environment management system. For details on
**conda** please refer to the [conda
manual](https://docs.conda.io/en/latest/).

This workflow collection makes heavy use of **conda** and especially
the [bioconda](https://bioconda.github.io) channel.

To create a working environment for this repository please install the
**MONSDA.yaml** environment (if not installed via **bioconda**) as found in the **envs** directory
like so:

**
conda env create -n monsda -f envs/MONSDA.yaml
**

The **envs** directory holds all the environments needed to run the pipelines in the **workflows** directory,
these will be installed automatically when needed.

For fast resolve of conda packages, we recommend **mamba**
[mamba](https://github.com/TheSnakePit/mamba)

For distribution of jobs one can either rely on local hardware, use
scheduling software like
[Slurm](https://slurm.schedmd.com/documentation.html) or the
[SGE](https://docs.oracle.com/cd/E19957-01/820-0699/chp1-1/index.html)
or follow any other integration of
[Snakemake](https://snakemake.readthedocs.io/en/stable/executing/cluster-cloud.html)
or
[Nextflow](https://www.nextflow.io/docs/latest/awscloud.html#aws-batch)
but be aware that most of these have not been tested for this
repository and usually require additional system dependent setup and
configuration.

This manual will only show examples on local and SLURM usage, but more
information on how to use other scheduling software is available under
the link above.  We also provide an example for SGE integration, this
however dates back to the times before **snakemake profiles**.
