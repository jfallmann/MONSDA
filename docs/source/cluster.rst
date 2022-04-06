.. _slurm:

============
Run on Slurm
============

Snakemake
=========

You can either use the slurm profile adapted from
Snakemake-Profiles_ that comes with this repository, or go
through the process of manually creating one, either using the cookiecutter example in the
``Snakemake-Profiles`` repository or on your own. You can also adapt the example that comes with this repository and execute

.. _Snakemake-Profiles: https://github.com/Snakemake-Profiles/slurm

.. code-block::

    monsda -j ${cpus} --configfile ${config.json} --directory ${PWD} --profile ${path_to_slurm_profile}


Further adaptions like grouping of jobs and advanced configs for rule
based performance increase will be tackled in future releases of ``MONSDA``.

Nextflow
========

Cluster config for Nextflow follows the description Nextflow-Executors_ and Nextflow-Profiles_. To use ``SlURM`` as executor you can adapt the profile that comes with this repository and simply append 

.. code-block::
    
    export NXF_EXECUTOR=slurm
    
to the call to ``MONSDA``.

.. _Nextflow-Executors: https://www.nextflow.io/docs/latest/executor.html
.. _Nextflow-Profiles: https://www.nextflow.io/docs/latest/config.html#config-profiles 