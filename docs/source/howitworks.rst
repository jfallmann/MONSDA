================
HOW DOES IT WORK
================

For details on ``Snakemake`` and it's features please refer to the
[snakemake
documentation](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html).
For details on ``Nextflow`` and it's features please refer to the
[nextflow
documentation](https://www.nextflow.io/docs/latest/index.html).

In general it is necessary to write a configuration file containing
information on paths, files to process and settings beyond default for
mapping tools and others.  The template on which analysis is based can
be found in the ``config`` directory.


``MONSDA`` acts a wrapper around ``snakemake`` or ``nextflow`` and the corresponding``config.json`` file.  The ``config.json`` holds all the information that is needed to run the jobs and will be parsed by ``MONSDA`` and split into sub-configs that can later be found in the directory ``SubSnakes`` or ``SubFlows`` respectively.

To successfully run an analysis pipeline, a few steps have to be followed:
  * Install MONSDA either via ``bioconda`` or ``pip``` folling the instruction in INSTALL
  * Directory structure: The structure for the directories is dictated by the ICS (Identifier-Condition-Setting relationship) in the config file
  * Config file: This is the central part of the analysis. Depending on this file ``MONSDA`` will determine processing steps and generate according config and ``snakemake`` subworkflow files to run each subworkflow until all processing steps are done.
