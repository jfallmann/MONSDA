START THE PIPELINE
==================

Simply run

```
python snakes/RunSnakemake.py
```

to see the help and available options that will be passed through to ```snakemake```.

To start a simple run call
```
python snakes/RunSnakemake.py -j NUMBER_OF_CORES --configfile YOUR_CONFIG.json --directory ${PWD}
```
or add additional arguments for ```snakemake``` as you see fit.
