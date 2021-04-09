===================
Start a pipline run
===================


Snakemake
---------

Run

``
python NextSnakes/NextSnakes.py
``

to see the help and available options that will be passed through to ``snakemake``.

To start a simple run call

``
python NextSnakes/NextSnakes.py -j NUMBER_OF_CORES --configfile YOUR_CONFIG.json --directory ${PWD}
``
or add additional arguments for ``Snakemake`` as you see fit.


Nextflow
--------

Run

``
python NextSnakes/NextSnakes.py
``

to see the help and available options that will be passed through to ``Nextflow``.

To start a simple run call

``
python NextSnakes/NextSnakes.py -c YOUR_CONFIG.json -d ${PWD}
``
or add additional arguments for ``Nextflow`` as you see fit.
