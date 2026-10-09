# CLI and validation recipe

Full docs: `docs/source/runsmk.rst`, `docs/source/cluster.rst`. This covers
what an agent needs to generate, dry-run, and then execute a config safely.

## Environment

```bash
conda activate monsda   # or whatever env MONSDA was installed into
```

`monsda` (Snakemake backend, default) and `monsda --nextflow` both need
`snakemake`/`nextflow` on PATH respectively; MONSDA aborts immediately if the
selected engine's binary is missing.

## `monsda` flags that matter for an agent

```
monsda -c CONFIG.json -d DIRECTORY -j PROCS [options]
```

- `-c/--config/--configfile CONFIG [CONFIG ...]`: one JSON, or JSON + an
  extra Nextflow-specific config.
- `-d/--directory`: the run/project root (where `FASTQ/`, `GENOMES/`,
  `SubSnakes|SubFlows/`, `JOBS/`, `LOGS/` live or will be created).
- `-j/--procs`: parallelism, capped by the config's `MAXTHREADS`.
- `--save`: **generate only** — write the sub-workflow files and the CLI
  calls to `JOBS/` without executing anything. This is the key step for
  validating a config cheaply.
- `-s/--skeleton`: just create the minimal directory hierarchy, no config
  processing.
- `--snakemake` (default) / `--nextflow`: pick the execution engine.
- `-u/--use-conda`: use conda envs (default true).
- `-l/--unlock`: unlock a Snakemake working directory after an interrupted run.
- `--clean`: Nextflow work-dir cleanup (`-n` to preview, `-f` to actually delete).
- `--loglevel {WARNING,ERROR,INFO,DEBUG}`.
- `--samplesheet FILE`: populate `SETTINGS` from a CSV/TSV when the config
  lacks one (see `reference/config-guide.md`); writes
  `<config>_with_settings.json` on first use.
- `--oras-registry HOST` / `--oras-namespace NAMESPACE`: override the
  container registry/namespace for ORAS image pulls (defaults
  `ghcr.io` / `jfallmann/monsda`).
- `--version`: print and exit.

Additional unrecognized arguments are passed through to Snakemake/Nextflow
verbatim, e.g. `-n` (Snakemake dry-run), `--omit-from <rule>` (skip a rule and
everything downstream without removing it from `WORKFLOWS`), `--profile
<dir>` (cluster profile).

## Validate-before-run recipe (do this before a real execution)

1. Run the repo's lightweight structural check first:
   `python scripts/validate_config.py CONFIG.json` (see this skill's
   `scripts/` directory) — catches invalid JSON, missing step blocks, a
   `VERSION` mismatch, and `SETTINGS`/`REFERENCE`/`ANNOTATION` paths that
   don't exist on disk, before anything talks to Snakemake/Nextflow.
2. Generate without executing:
   ```bash
   monsda -c CONFIG.json -d DIRECTORY -j 2 --save
   ```
   This writes the per-condition sub-workflows (`SubSnakes/*.smk` or
   `SubFlows/*.nf`) and the exact CLI invocations MONSDA would run into
   `JOBS/MONSDA.commands`.
3. Dry-run every generated sub-workflow using the commands recorded in
   `JOBS/MONSDA.commands`, adding `-n` (Snakemake) — e.g.:
   ```bash
   snakemake --snakefile SubSnakes/<name>.smk --configfile <subconfig>.json -n
   ```
   For Nextflow, `nextflow run <subflow>.nf -preview` only validates syntax,
   not the script bodies — a real `-n`-equivalent dry run is not available,
   so inspect the generated `.nf`/its `params` instead, or do a tiny real run
   on a minimal subset of samples first.
4. Fix whatever the dry run reports (commonly `MissingInputException` for
   files that don't exist yet, or a bad tool/env name) and repeat from step 2
   before launching the full, real run.

## Running for real

```bash
monsda -c CONFIG.json -d DIRECTORY -j THREADS --conda-frontend mamba \
  --conda-prefix PATH_TO_CONDA_ENVS
```

or, for Nextflow:

```bash
monsda --nextflow -c CONFIG.json -d DIRECTORY -j THREADS
```

### On a cluster (SLURM example)

```bash
monsda -c CONFIG.json -d DIRECTORY -j THREADS --conda-frontend mamba \
  --profile profile_snakemake --conda-prefix PATH_TO_CONDA_ENVS
```

```bash
export NXF_EXECUTOR=slurm
monsda --nextflow -c CONFIG.json -d DIRECTORY -j THREADS
```

`profile_snakemake/` and `profile_nextflow/` in this repo ship example
profiles to adapt.

### Resuming a run while skipping a step

Keep the full original `WORKFLOWS` list (see `reference/config-guide.md`),
and pass Snakemake's own skip flag instead of editing the config:

```bash
monsda -c CONFIG.json -d DIRECTORY -j THREADS --omit-from multiqc
```
