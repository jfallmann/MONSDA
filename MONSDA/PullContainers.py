#!/usr/bin/env python3
"""Pre-pull all apptainer/singularity containers referenced by a generated
MONSDA (sub-)Snakefile, without creating conda environments.

This mirrors what ``snakemake --conda-create-envs-only`` does for conda, but
for containers. It exists because snakemake's own
``--conda-create-envs-only`` unconditionally implies the conda deployment
method and then always calls ``DAG.create_conda_envs()`` even when only
``--use-apptainer``/``--use-singularity`` was requested. If the rules also
declare a ``container:`` directive (as all MONSDA rules do), snakemake then
tries to create the per-job conda environment *inside* the container via
``singularity/apptainer exec ... conda info --json``, which fails with
``conda: not found`` unless conda happens to be on PATH inside that
container (see snakemake/snakemake#3088, #3202).

This script builds the same DAG snakemake would, but calls
``DAG.pull_container_imgs()`` directly and neutralizes
``DAG.create_conda_envs()`` for the duration of the call, so only the
containers are pulled.
"""

import argparse
import sys
from pathlib import Path

from snakemake.api import SnakemakeApi
from snakemake.dag import DAG
from snakemake.settings.types import (
    ConfigSettings,
    DAGSettings,
    DeploymentMethod,
    DeploymentSettings,
    OutputSettings,
    ResourceSettings,
)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Pull apptainer/singularity containers for a MONSDA (sub-)Snakefile without creating conda envs."
    )
    parser.add_argument("-s", "--snakefile", required=True, type=Path)
    parser.add_argument("--configfile", required=True, type=Path)
    parser.add_argument("-j", "--cores", "--jobs", dest="cores", type=int, default=1)
    parser.add_argument("-d", "--directory", type=Path, default=None)
    parser.add_argument("--apptainer-prefix", "--singularity-prefix", dest="apptainer_prefix", default=None)
    parser.add_argument("--apptainer-args", "--singularity-args", dest="apptainer_args", default="")
    args, _unknown = parser.parse_known_args(argv)
    return args


def pull_containers(snakefile, configfile, cores=1, directory=None, apptainer_prefix=None, apptainer_args=""):
    """Build the DAG for snakefile/configfile and pull all referenced containers."""
    orig_create_conda_envs = DAG.create_conda_envs
    DAG.create_conda_envs = lambda self, *a, **k: None
    try:
        with SnakemakeApi(OutputSettings()) as snakemake_api:
            workflow_api = snakemake_api.workflow(
                resource_settings=ResourceSettings(cores=cores),
                config_settings=ConfigSettings(configfiles=[configfile]),
                deployment_settings=DeploymentSettings(
                    deployment_method={DeploymentMethod.APPTAINER},
                    apptainer_prefix=apptainer_prefix,
                    apptainer_args=apptainer_args,
                ),
                snakefile=snakefile,
                workdir=directory,
            )
            dag_api = workflow_api.dag(dag_settings=DAGSettings())
            dag_api.conda_create_envs()
    finally:
        DAG.create_conda_envs = orig_create_conda_envs


def main(argv=None):
    args = parse_args(argv)
    pull_containers(
        snakefile=args.snakefile,
        configfile=args.configfile,
        cores=args.cores,
        directory=args.directory,
        apptainer_prefix=args.apptainer_prefix,
        apptainer_args=args.apptainer_args,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
