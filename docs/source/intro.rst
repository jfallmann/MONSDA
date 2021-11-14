==========
MONSDA
==========

Welcome to '''MONSDA''', Modular Organizer of Nextflow and Snakemake driven hts Data Analysis

Automizing hts analysis from data download, preprocessing and mapping to postprocessing/analysis and track generation centered on a single config file. '''MONSDA''' can create ```snakemake``` and ```nextflow``` workflows centered on a user friendly, sharable ``Json`` config file and reproducible subworkflows. These workflows can either be saved to disk for manual inspection and execution or automatically executed.

For details on ```snakemake``` and ```nextflow``` and their features please refer to the corresponding [snakemake](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html)  or [nextflow](https://www.nextflow.io/docs/latest/index.html) documentation.

In general it is necessary to write a configuration file containing workflows to execute, information on paths, files to process and settings beyond default for mapping tools and others.
The template on which ```MONSDA``` is based on can be found in the ```config``` directory.

For '''MONSDA''' to be as FAIR as possible, one needs to use ```conda``` or the faster drop-in replacement ```mamba```. For details on either please refer to the corresponding [conda](https://docs.conda.io/en/latest/) or [mamba](https://mamba.readthedocs.io/en/latest/) manual.

This workflow organizer makes heavy use of ```conda``` and especially the [bioconda](https://bioconda.github.io) channel.
