=======
Install
=======

Install MONSDA via ```conda``` or ```pip```
-------------------------------------------

To install via ```conda/mamba``` simply run

```
mamba install -c bioconda -c conda-forge monsda
```

To install via ```pip``` you first need to create the ```MONSDA``` environment as found in the ```envs``` directory of this repository like so:

```
mamba env create -n monsda -f MONSDA/envs/monsda.yaml
```

The ```envs``` directory holds all the environments needed to run the pipelines in the ```workflows``` directory, these will be installed automatically alongside ```MONSDA```.

For that activate the ```MONSDA``` environment and run ```pip```

```
conda activate monsda
pip install MONSDA
```

Install from source
-------------------

Simply clone this repository with ``git clone``.

You can then install dependencies as described for ```pip``` installation and manually run ```setup.py```.
Be aware that ```MONSDA``` is version dependent, so config files can only be run with the specified version of '''MONSDA''' in order to guarantee reproducibility by conserving dependencies and environments.