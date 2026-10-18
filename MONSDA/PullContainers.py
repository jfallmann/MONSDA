#!/usr/bin/env python3
"""Pre-pull all apptainer/singularity containers referenced by a generated
MONSDA (sub-)Snakefile, without creating conda environments.

This mirrors what ``snakemake --conda-create-envs-only`` does for conda, but
for containers. It exists because snakemake's own
``--conda-create-envs-only`` (``DAGApi.conda_create_envs()``) unconditionally
implies the conda deployment method *before* the DAG is even built. Once
conda is implied, building the DAG itself (``Job.get_wait_for_files()``)
accesses each job's conda env address/hash, which for a rule that also
declares a ``container:`` directive (as all MONSDA rules do) lazily runs
``conda info --json`` *inside* the container via
``singularity/apptainer exec``. That fails with ``conda: not found`` unless
conda happens to be on PATH inside the container (see
snakemake/snakemake#3088, #3202) - regardless of whether
``DAG.create_conda_envs()`` itself ever actually runs.

This script therefore never lets CONDA enter the deployment method at all:
it builds the DAG by hand with only APPTAINER set (the same private
``Workflow._prepare_dag()``/``_build_dag()`` calls snakemake itself uses for
e.g. ``--unlock``) and then calls ``DAG.pull_container_imgs()`` directly.
"""

import argparse
import sys
from pathlib import Path

from snakemake.api import SnakemakeApi
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
    """Build the DAG for snakefile/configfile and pull all referenced containers.

    Deliberately avoids DAGApi.conda_create_envs()/WorkflowApi.dag(...).conda_create_envs(),
    both of which unconditionally imply DeploymentMethod.CONDA before building the DAG
    (see snakemake.api.DAGApi.conda_create_envs). With CONDA implied, DAG construction
    itself (Job.get_wait_for_files()) accesses each rule's conda env address/hash, which
    for a non-externally-managed env lazily instantiates Conda(container_img=...) and
    runs 'conda info --json' *inside* the container via apptainer/singularity exec -
    exactly the failure --apptainer-create-envs-only exists to avoid. So here the DAG is
    built by hand with only APPTAINER in the deployment method, the same private
    _prepare_dag()/_build_dag() calls snakemake itself uses for e.g. --unlock, and then
    only DAG.pull_container_imgs() is invoked.
    """
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
        workflow_api.dag(dag_settings=DAGSettings())
        workflow = workflow_api._workflow
        workflow._prepare_dag(forceall=False, ignore_incomplete=True, lock_warn_only=True)
        workflow._build_dag()
        workflow.dag.pull_container_imgs()


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
