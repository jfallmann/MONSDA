INSTALL
=======

Contains sumodule so clone with ```git clone --recursive``` or if already cloned pull submodules with
```
git submodule init
git submodule update
```

See [stackoverflow](https://stackoverflow.com/questions/25200231/cloning-a-git-repo-with-all-submodules) and
[git docs](https://git-scm.com/book/en/v2/Git-Tools-Submodules) for details.

For details on ```snakemake``` and it's features please refer to the [snakemake documentation](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html).

In general it is necessary to write a configuration file containing information on paths, files to process and settings beyond default for mapping tools and others.
For examples on this please have a look into the ```config``` directory.

For ```snakemake``` to be fully FAIR, one needs to use ```conda``` or similar environment management systems. For details on ```conda``` please refer to the [conda manual](https://docs.conda.io/en/latest/).

This workflow collection makes heavy use of ```conda``` and especially the [bioconda](https://bioconda.github.io) channel.

To create a working environment for this repository please install the snakemake environment as found in the ```envs``` directory like so:

```
conda env create -n snakemake -f snakes/envs/snakemake.yaml
```

The ```envs``` directory holds all the environments needed to run the pipelines in the ```workflows``` directory, these will be installed automatically when needed.

For distribution of jobs one can either rely on local hardware, use scheduling software like [Slurm](https://slurm.schedmd.com/documentation.html) or the [SGE](https://docs.oracle.com/cd/E19957-01/820-0699/chp1-1/index.html) or follow any other integration of [Snakemake](https://snakemake.readthedocs.io/en/stable/executing/cluster-cloud.html) although most of these are not tested for this repository.

This manual will only show examples on local and SLURM usage, but more information on how to use other scheduling software is available under the link above.
We also provide an example for SGE integration, this however dates back to the times before ```snakemake``` profiles.
