================
HOW DOES IT WORK
================

This repository hosts the executable ``MONSDA.py`` which acts a wrapper around ``snakemake`` and the
``config.json`` file.  The ``config.json`` holds all the information that is needed to run the jobs and will
be parsed by ``MONSDA.py`` and split into sub-configs that can later be found in the directory
``SubSnakes``.

To successfully run an analysis pipeline, a few steps have to be followed:
  * Clone this repository to you working directory, or symlink it there. *DO NOT CHANGE THE NAME OF THE REPO!!!* This is needed for subworkflow generation.
  * Directory structure: The structure for the directories is dictated by the ICS (Identifier-Condition-Setting relationship) in the config file
  * Config file: This is the central part of the analysis. Depending on this file ``MONSDA.py`` will determine processing steps and generate according config and ``snakemake`` subworkflow files to run each subworkflow until all processing steps are done.
