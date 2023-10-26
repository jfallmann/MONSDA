===================
Start a pipline run
===================


Snakemake
---------

Activate the **MONSDA** conda environment and run


.. code-block::
    
    monsda --help


to see the help and available options that will be passed through to **Snakemake**.

To start a job with **Snakemake**, which is the default, run

.. code-block::

    monsda -j NUMBER_OF_CORES -c YOUR_CONFIG.json --directory ${PWD}


or add additional arguments for **Snakemake** as you see fit,
**Snakemake** currently defaults to mamba as conda frontend. Please be aware that for that to work one should follow the recommendations at MAMBA_. However, using conda-libmamba-solver, the conda frontend can lead to an even better and more stable experience. We currently recommend to set a fixed directory to store environments (here conda_envs) and run the conda frontend.  

.. _MAMBA: https://mamba.readthedocs.io/en/latest/mamba-installation.html

.. code-block::
    
    --conda-frontend conda --conda-prefix path_to_conda_envs


Nextflow
--------

To run **MONSDA** in **Nextflow** mode just add '--nextflow'

.. code-block::

    monsda --nextflow -j NUMBER_OF_CORES -c YOUR_CONFIG.json --directory ${PWD}


As with **Snakemake** additional arguments for **Nextflow** can be added and will be passed through.
