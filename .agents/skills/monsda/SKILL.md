---
name: monsda
description: Use MONSDA (Modular Organizer of Nextflow and Snakemake driven hts Data Analysis) to plan and run HTS/NGS pipelines — QC, trimming, mapping, dedup, counting, tracks, peak calling, DE/DEU/DAS/DTU, circRNA, fusion detection — from a single config.json that MONSDA turns into Snakemake or Nextflow workflows. Use whenever the user wants to analyze FASTQ/BAM sequencing data with this repo, needs a MONSDA config built or edited, needs an experimental design mapped to a condition-tree, or wants to generate/dry-run/execute the resulting workflows instead of doing it by hand.
---

# MONSDA

MONSDA is a config-driven wrapper: one JSON config describes samples, genomes
and tool choices; MONSDA compiles it into per-condition Snakemake or Nextflow
sub-workflows and can run them. Your job as the agent is to turn the user's
experiment description into a correct config, validate it cheaply before
spending compute, and only then run it — instead of the user hand-editing
JSON and guessing at tool names.

Read the reference files below on demand (do not load all of them up front):

- `reference/config-guide.md` — config.json anatomy, condition-tree, OPTIONS/TOOLS, samplesheet shortcut, and the `SAMPLES` naming rules for `single`/`paired`/`singlecell` data including Illumina BCL Convert names. Read before writing or editing any config.
- `reference/workflows.md` — available `WORKFLOWS` steps and the tools each one accepts, plus the `RAW/<condition-path>` raw signal layout for `BASECALL` (Guppy FAST5, Dorado POD5). Read when choosing which steps/tools to enable.
- `reference/cli.md` — `monsda` / `monsda_configure` CLI flags, the generate-then-dry-run validation recipe, cluster profiles. Read before invoking any command.

## Decision flow

1. **Clarify the experiment, not the tool syntax.** Ask the user (or infer from
   what they already said): experimental design/conditions, where the
   FASTQ/BAM files and genome/annotation live, sequencing type (single/paired/
   single-cell, stranded?), and which analysis steps they want
   (QC/TRIMMING/MAPPING/DEDUP/COUNTING/TRACKS/PEAKS/DE/DEU/DAS/DTU/CIRCS/
   FUSIONS). Don't ask about MONSDA internals they aren't expected to know.
2. **Map the design to a condition-tree** (`ID -> CONDITION -> SETTING`,
   see `reference/config-guide.md`). This fixes the expected input directory
   layout (`FASTQ/<ID>/<CONDITION>/...`) and the output tool-key directories.
3. **Build the config.** Prefer, in order:
   - If the user has (or can export) a sample table, use the built-in
     `--samplesheet file.csv` support (`configs/samplesheet_template.csv`
     shows the columns) instead of hand-writing the nested `SETTINGS` block.
   - Otherwise copy the closest starting point — `configs/template_clean.json`
     (empty skeleton) or `configs/tutorial_quick.json` (small worked example)
     — and edit it. Several other files under `configs/` (`template.json`,
     `tutorial_toolmix.json`, `tutorial_exhaustive.json`,
     `tutorial_postprocess.json`) carry trailing `#comments` for human
     documentation and are **not valid JSON**; read them for explanations but
     don't copy them as-is — see `reference/config-guide.md` for the full list.
   - Keep `WORKFLOWS` as a single comma-separated string in the **order**
     steps should run, and give every listed step a corresponding top-level
     key (`QC`, `MAPPING`, ...), each with a `TOOLS` map and per-condition
     `OPTIONS` — MONSDA errors out otherwise.
   - Set `VERSION` to the exact installed MONSDA version the run will use
     (`python -c "import MONSDA; print(MONSDA.__version__)"` with the same
     interpreter/env as the run) — a mismatch aborts the run.
4. **Validate before executing anything expensive.** Run
   `scripts/validate_config.py <config.json>` (does structural/JSON checks and
   cross-checks referenced sample/reference files) and then the
   generate-only + dry-run recipe in `reference/cli.md` (`monsda --save` then
   `snakemake -n` / nextflow equivalent on the generated sub-workflows). This
   catches missing files, bad tool names and condition-tree mistakes without
   running a single real job.
5. **Run it.** Only after a clean dry run, execute for real (local, or via a
   cluster profile — see `reference/cli.md`). Report the actual command used
   and where outputs land (tool-key directories under e.g. `MAPPING/<id>/
   <condition>/<combo>/`).

## Preprocessing

UMIs that sit in a separate read instead of inline in R1/R2 have to be
prepended to the mate before MONSDA sees the data; MONSDA has no step for it.

- From the repo root:
  `python scripts/Preprocessing/addUMI.py --umi-fastq sample_UMI.fastq.gz sample_R1.fastq.gz`
- Choose **one** mate according to the library chemistry. If the protocol
  calls for the UMI on R2, pass `sample_R2.fastq.gz` instead. Configure
  extraction for that same mate; the script does not infer the chemistry.
- The script prepends the full UMI sequence *and* its base qualities at the 5'
  end and keeps the original read header unchanged. UMI and read files must be
  in the same order with matching names (dnaio validates this and aborts on a
  mismatch); requires `dnaio >= 0.9`.
- `--umi-fastq` is not optional for real UMIs: without it the script generates
  random UMIs, which are synthetic and useless for deduplication.
- Output lands next to the input as `<name>_umi.fastq.gz`; raw files are never
  overwritten.
- Feed it to MONSDA as a separate run: stage the modified mate under
  `FASTQ/<ID>/<condition>/` with the expected name (`sample_R1.fastq.gz` or
  `sample_R2.fastq.gz`) and the untouched other mate under its original
  matched basename. Do not stage the UMI FASTQ as a mate.

Run `addUMI.py` before MONSDA; do not add it to `WORKFLOWS`. For the MONSDA
`DEDUP` step, select `umitools` (`umi_tools`) and set `OPTIONS.EXTRACT` to the
`umi_tools extract` barcode pattern for the mate and UMI length you prepended
(`--bc-pattern` with one `N` per UMI base for R1, `--bc-pattern2` for R2), then
map and dedup the BAM as usual. With `EXTRACT` empty MONSDA links the FASTQs
through unchanged.

## What this buys the user over doing it by hand

- No need to memorize the condition-tree <-> directory-layout mapping or the
  per-tool `OPTIONS` schema — you look it up in the reference files.
- `scripts/validate_config.py` catches config mistakes before a conda env
  solve or a multi-hour mapping job starts.
- You generate-and-dry-run (`--save` + `snakemake/nextflow -n`) before
  committing to a real run, which is the step people skip when working by
  hand and MONSDA doesn't enforce.

## Notes

- `monsda_configure` is an interactive terminal UI meant for humans; it is
  not practical to drive from an agent. Build/edit the JSON directly (or via
  `--samplesheet`) instead of trying to script the TUI.
- When testing changes to MONSDA itself, do not run pytest in the working tree.
  By default run `bash tests/run_installed_tests.sh [pytest args]` (or
  `make test-installed`): it builds the wheel, installs it into a conda env
  made from `envs/monsda_min.yaml`, and runs the whole suite against the
  installed package from a staged copy of `tests/`. Set `MONSDA_TEST_ENV` to
  reuse an existing env, `MONSDA_TEST_FRESH=1` to rebuild it. On a cluster
  submit it as a Slurm job with `--time`. Report which tests failed because of
  missing R packages (GSVA, decoupleR) separately from real failures.
- If you are working inside a MONSDA git checkout rather than an installed
  release, see this repo's `AGENTS.md` for the dev-only gotchas (symlinks
  under `MONSDA/MONSDA`, PATH for `snakemake`/`nextflow`, golden-file tests).
  None of that applies to an end user running an installed `monsda`.
