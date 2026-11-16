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

    export PATH=${HOME}/anaconda3/envs/monsda/bin:$PATH

Confirmed there: snakemake 9.20.0, nextflow 26.04.0.12031. `multiqc` is **not**
installed in that env, MultiQC behaviour cannot be verified locally.

That env contains an installed copy of MONSDA in `site-packages`, *not* an
editable install. It only picks up the working tree when the cwd is the repo, so
runs from anywhere else need `PYTHONPATH=<repo>` or they silently test the
installed release instead of the change.

MONSDA aborts early if `nextflow` is missing, so `--nextflow` runs need that
PATH.

### Running the test suite (always: against the installed wheel)

This is the only way to run the tests. Never run pytest in the working tree or
in `tests/` of the checkout, it pollutes the repo with caches, generated
workflows and logs. Build the package, install it into the minimal env
(`envs/monsda_min.yaml`) and run the tests outside the repo. Use

    bash tests/run_installed_tests.sh [pytest args]     # or: make test-installed

It builds a wheel from the checkout, installs it into a conda env created from
`envs/monsda_min.yaml` (default `/tmp/monsda_installed_test/env`, reused on later
runs), copies `tests/` to a stage dir whose `workflows/ envs/ scripts/ configs/`
are symlinks into the installed `share/MONSDA`, and runs pytest there with
`PYTHONPATH` unset. The tests therefore exercise the packaged code and data
files, which also catches files missing from `setup.py` `data_files`.

Env vars: `MONSDA_TEST_HOME` (work dir), `MONSDA_TEST_ENV` (existing env prefix,
skips creation), `MONSDA_TEST_ENVFILE`, `MONSDA_TEST_FRESH=1` (recreate env).

This is a public repository: do not commit personal infrastructure here — no
host aliases, user names, or cluster paths. If a large-data run on a personal
HPC cluster is needed, keep that recipe in the *user's* global agent notes, not
in this file, and run the suite there inside a batch job (never on a login
node). The R-dependent PostDE tests skip when no `Rscript` with GSVA/decoupleR
is on PATH; if an unrelated `Rscript` without those packages is found, they
fail — that is an environment problem, not a code problem. The wheel build may
warn about `MONSDA.lib` package-data discovery; harmless.

### Gotchas when running MONSDA from a checkout

(Only for generating/dry-running configs from the tree, not for the test suite.)

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

- Workflow templates are snapshot tested. A golden mismatch is a test failure,
  not something to refresh on every run: `UPDATE_GOLDEN` is opt-in and must
  only be used after the diff has been reviewed and shown to be an intended
  template change (the version string is already normalized to `VERSION`, so
  version bumps never need it). Never set it by default or to make a red run
  green. Regenerate with
  `UPDATE_GOLDEN=1 bash tests/run_installed_tests.sh -k golden`. That rewrites
  the golden files in the stage copy, so copy them back with
  `rsync -a ${MONSDA_TEST_HOME:-/tmp/monsda_installed_test}/stage/tests/golden/ tests/golden/`
  and check `git diff tests/golden` only contains the intended lines.
