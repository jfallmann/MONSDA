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

- `reference/config-guide.md` — config.json anatomy, condition-tree, OPTIONS/TOOLS, samplesheet shortcut. Read before writing or editing any config.
- `reference/workflows.md` — available `WORKFLOWS` steps and the tools each one accepts. Read when choosing which steps/tools to enable.
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
- If you are working inside a MONSDA git checkout rather than an installed
  release, see this repo's `AGENTS.md` for the dev-only gotchas (symlinks
  under `MONSDA/MONSDA`, PATH for `snakemake`/`nextflow`, golden-file tests).
  None of that applies to an end user running an installed `monsda`.
