============
First Steps
============

**MONSDA** acts a wrapper around **Snakemake** or **Nextflow** based on a user defined **config.json** file.  This **config.json** holds all the information that is needed to run the jobs and will be parsed by **MONSDA** and split into independent sub-configs that can later be found in the directory **SubSnakes** or **SubFlows** respectively. Command line execution calls are stored in the directory *JOBS*, so users can manipulate those or rerun them manually as needed. By default, however, **MONSDA** will run those jobs automatically either locally, or through **Snakemake's** or **Nextflow's** integrated cluster interfaces.

To successfully run an analysis pipeline, a few steps have to be followed:
  * Install MONSDA either via **bioconda** or **pip** following the instruction in :ref:`install`
  * Directory structure: The structure for the directories is dictated by :ref:`condition-tree` in the config file
  * Config file: This is the central part of a **MONSDA** run. Depending on :ref:`config-file` **MONSDA** will determine processing steps and generate corresponding config and workflow files to run each subworkflow until all processing steps are done.


In general it is necessary to write a configuration file containing
information on paths, files to process and settings beyond default for
mapping tools and others.  The template on which analysis is based on can
be found in the **config** directory and will be explained in detail later.

To create a working environment for this repository please install the
**MONSDA.yaml** environment (if not installed via **bioconda**) as found in the **envs** directory
like so:

::

  conda env create -n monsda -f envs/MONSDA.yaml

The **envs** directory holds all the environments needed to run the pipelines in the **workflows** directory,
these will be installed automatically when needed.

For fast resolve of conda packages, we recommend conda-libmamba-solver_ which is a new solver for the conda package manager and speeds up conda without the need to install mamba  and is shipped with **MONSDA**. However, the user if free to use mamba_ which is currently also the standard conda-frontend for Snakemake_. 

.. _mamba: https://mamba.readthedocs.io/en/latest/
.. _conda-libmamba-solver: https://github.com/conda-incubator/conda-libmamba-solver

For distribution of jobs one can either rely on local hardware, use
scheduling software like
Slurm_ or the SGE_
or follow any other integration in
Snakemake_ or Nextflow_
but be aware that most of these have not been tested for this
repository and usually require additional system dependent setup and
configuration.

.. _Slurm: https://slurm.schedmd.com/documentation.html
.. _SGE: https://docs.oracle.com/cd/E19957-01/820-0699/chp1-1/index.html
.. _Snakemake: https://Snakemake.readthedocs.io/en/stable/executing/cluster-cloud.html
.. _Nextflow: https://www.Nextflow.io/docs/latest/awscloud.html#aws-batch

This manual will only show examples on local and SLURM usage.
