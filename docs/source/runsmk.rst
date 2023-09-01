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
we currently recommend to set mamba as conda frontend and set a fixed directory to store environments (here conda_envs). However, using conda-libmamba-solver, the conda frontend can lead to an even better and more stable experience.

.. code-block::
    
    --conda-frontend mamba --conda-prefix path_to_conda_envs


Nextflow
--------

To run **MONSDA** in **Nextflow** mode just add '--nextflow'

.. code-block::

    monsda --nextflow -j NUMBER_OF_CORES -c YOUR_CONFIG.json --directory ${PWD}


As with **Snakemake** additional arguments for **Nextflow** can be added and will be passed through.
