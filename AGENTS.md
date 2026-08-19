# Notes for AI agents

## Keep the user updated

Report progress **while** working, not only at the end. Before/while running a
step, say what is being done and why, and show the relevant command output or
diff. A silent stretch of work followed by one large summary makes the change
impossible to judge and is not acceptable.

State clearly what was actually verified and what was not.

## Test environment

`snakemake` and `nextflow` live in the `monsda` conda env, they are *not* on the
default PATH:

    export PATH=/home/fall/anaconda3/envs/monsda/bin:$PATH

Confirmed there: snakemake 9.20.0, nextflow 26.04.0.12031. `multiqc` is **not**
installed in that env, MultiQC behaviour cannot be verified locally.

That env contains an installed copy of MONSDA in `site-packages`, *not* an
editable install. It only picks up the working tree when the cwd is the repo, so
runs from anywhere else need `PYTHONPATH=<repo>` or they silently test the
installed release instead of the change.

MONSDA aborts early if `nextflow` is missing, so `--nextflow` runs need that
PATH.

### Gotchas when running MONSDA from a checkout

- `envpath`/`workflowpath` resolve to `<repo>/MONSDA/MONSDA/{envs,workflows}`,
  which does not exist in a checkout. Create the symlinks before running and
  remove them afterwards, they must not be committed:

      mkdir -p MONSDA/MONSDA
      ln -sfn ../../workflows MONSDA/MONSDA/workflows
      ln -sfn ../../envs MONSDA/MONSDA/envs
      ln -sfn ../../scripts MONSDA/MONSDA/scripts

- `VERSION` in the config has to match the installed version exactly, otherwise
  MONSDA refuses to run. Get it from `MONSDA.__version__` of the interpreter
  used for the run, it differs between envs.

- Generating without executing: run from the target directory with
  `python -m MONSDA.RunMONSDA -c <config> -d <dir> -j 2 --save`. The declared
  console script is `monsda`, not `MONSDA`, and it is not on PATH in the
  `monsda` env, so use the module form and set `PYTHONPATH=<repo>` if the
  interpreter does not already import MONSDA from the tree. Then
  dry run every generated workflow via the calls in `JOBS/MONSDA.commands`
  (`snakemake ... -n`). Note that `tests/data` ships no BAMs, so most sub
  workflows fail with `MissingInputException` even on an unmodified tree,
  compare the failure set against a pristine worktree instead of expecting zero
  failures.

- Workflow templates are snapshot tested. On intended changes regenerate with
  `UPDATE_GOLDEN=1 python -m pytest` in `tests/` and check the golden diff only
  contains the intended lines.

- `${HOME}/MONSDA` is a symlink to this repo, the same tree under a second path.
  Editor diagnostics may report paths there, they refer to these files.

## Test data and configs

`${HOME}/Work2/Tests/MONSDATest` is the user's test directory. It holds real
input data (`FASTQ/{Ecoli,SIMPLE,FGUMI}`, `GENOMES/Ecoli` incl. prebuilt
indices, 122 BAMs in `MAPPED/`), their own configs at the top level and outputs
of earlier runs. Use it instead of `tests/data`, which ships no BAMs. It is not
a git repo, so nothing there is protected by version control.

Put configs written for a test into `TESTCONFIGS/` so they can be reused:

- `config_multiqc_fastqc.json`, `config_multiqc_rustqc.json`: QC, TRIMMING,
  DEDUP, MAPPING on `Ecoli/KO` with a `MULTI` option set. The rustqc one selects
  the `multiqc_rustqc` template, the fastqc one `multiqc`.
- `monsda_gen.sh <config> <outdir> [--nextflow]`: generates only (`--save`).
  Creates the `MONSDA/MONSDA` symlinks, stamps `VERSION` with the tree's
  `MONSDA.__version__` and links `FASTQ`/`GENOMES` into a sandbox `outdir`, so
  the outputs in the test root are not overwritten. MONSDA resolves input data
  relative to its working directory, hence the links.
- `mqc_script_test.sh <generated_subflow.nf> [combo] [condition] [params]`: runs
  the generated nextflow `mqc` process with a fake `multiqc` that records its
  argv, and checks the module list from `LOGS/versions.txt`, the scanned
  directories and the published files. `nextflow lint` and `-preview` do not
  render script bodies, so this is the only local check of that bash.

Do not run tests in the test root itself, generate into a subdirectory.
