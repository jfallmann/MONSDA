MONSDA
======

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: GETTING STARTED

   source/installation
   source/first

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: PREPARING YOUR PROJECT

   source/preparation
   source/configurator

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: EXECUTING MONSDA

   source/runsmk
   source/cluster

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: WORKFLOW AND TOOL OVERVIEW

   source/workflows

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: TUTORIAL

   source/tutorial

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: DETAILS

   source/wrapper
   source/conditiontree
   source/config

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: CONTRIBUTE

   source/integrate
   source/contribute


Welcome to ``MONSDA``, Modular Organizer of Nextflow and Snakemake driven hts Data Analysis

Automizing HTS analysis from data download, preprocessing and mapping to postprocessing/analysis and track generation centered on a single config file. ``MONSDA`` can create ``snakemake`` and ``nextflow`` workflows centered on a user friendly, sharable ``Json`` config file and reproducible subworkflows. These workflows can either be saved to disk for manual inspection and execution or automatically executed.

For details on ``snakemake`` and ``nextflow`` and their features please refer to the corresponding snakemake_  or nextflow_ documentation.

.. _snakemake: https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html
.. _nextflow: https://www.nextflow.io/docs/latest/index.html

In general it is necessary to write a configuration file containing workflows to execute, information on paths, files to process and settings beyond default for mapping tools and others.
The template on which ``MONSDA`` is based on can be found in the ``config`` directory.

For ``MONSDA`` to be as FAIR as possible, one needs to use ``conda`` or the faster drop-in replacement ``mamba``. For details on either please refer to the corresponding conda_ or mamba_ manual.

.. _conda: https://docs.conda.io/en/latest/
.. _mamba: https://mamba.readthedocs.io/en/latest/

This workflow organizer makes heavy use of ``conda`` and especially the bioconda_ channel.

.. _bioconda: https://bioconda.github.io

==================
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

MONSDA.rst (END)

