# config.json anatomy

Full human docs: `docs/source/config.rst`, `docs/source/conditiontree.rst`,
`docs/source/preparation.rst`. This is the condensed, agent-facing version.

## Top level

```json
{
  "WORKFLOWS": "QC,TRIMMING,MAPPING,DEDUP",
  "BINS": "",
  "MAXTHREADS": "20",
  "VERSION": "1.5.0",
  "SETTINGS": { ... },
  "QC": { ... },
  "TRIMMING": { ... },
  "MAPPING": { ... },
  "DEDUP": { ... }
}
```

- `WORKFLOWS`: comma-separated step names, in run order. Every step listed
  here needs a matching top-level key with its own settings, or MONSDA errors.
  Available steps: see `reference/workflows.md`.
- `VERSION` must equal the installed `MONSDA.__version__` for the interpreter
  that will run the job, or the run aborts.
- `MAXTHREADS` caps `-j` passed on the CLI.
- `BINS`: path to custom user scripts; advanced/rare, leave `""` otherwise.

### Rerunning on an existing output directory

Keep the **full, original** `WORKFLOWS` list even if you only want to rerun a
later step (e.g. keep `QC,TRIMMING,MAPPING,DEDUP` even if only `DE` changed).
The tool-key/`combo` directory name is derived from which stages are listed;
dropping an earlier stage changes that name and makes existing outputs look
missing, forcing a full rerun. To skip actually re-executing one rule, add
Snakemake's own `--omit-from <rule>` on the command line instead of editing
`WORKFLOWS`.

## The condition-tree (`SETTINGS`)

```
'ID' -> 'CONDITION' -> 'SETTING' (optional third level)
```

```json
"SETTINGS": {
  "LabA": {
    "wild-type": {
      "SAMPLES": ["Sample_1", "Sample_2"],
      "SEQUENCING": "paired",
      "REFERENCE": "GENOMES/Dm6/dm6.fa.gz",
      "ANNOTATION": {
        "GTF": "GENOMES/Dm6/dm6.gtf.gz",
        "GFF": "GENOMES/Dm6/dm6.gff3.gz"
      }
    }
  }
}
```

- Each leaf of the tree is one independent analysis unit. The tree also
  dictates the **input** directory layout: `FASTQ/<ID>/<CONDITION>/...`
  (MONSDA looks one level above the deepest leaf for input, so adding a
  `SETTING` sub-level for e.g. benchmarking different tool params does not
  require duplicating the FASTQ files).
- `SAMPLES`: list of sample basenames, no file extension. For paired-end,
  list the name once without `_R1`/`_R2` — MONSDA appends those
  automatically and requires files to end in `_R1.fastq.gz`/`_R2.fastq.gz`
  (single-end: `<name>.fastq.gz`).
- `SEQUENCING`: `single`, `paired`, or `singlecell`; optionally append
  comma-separated strandedness (`rf` or `fr`, see RSeQC `infer_experiment.py`
  convention) e.g. `"paired,fr"`. Leave strandedness empty if unstranded —
  MONSDA configures per-tool strandedness flags automatically from this.
- `REFERENCE` / `ANNOTATION.GTF` / `ANNOTATION.GFF` / `INDEX` / `PREFIX` /
  `DECOY`: per-condition defaults, can be overridden per workflow step.
- `GROUPS` / `BATCHES` / `TYPES`: optional, same length as `SAMPLES`, used
  only by DE/DEU/DAS/DTU design matrices (batch correction, type covariate).
- `IP`: only for `PEAKS` with CLIP-type protocols (`CLIP`, `iCLIP`, `revCLIP`).

Output directories are named by a derived **tool-key/combo**, e.g.
`MAPPING/LabA/wild-type/fastqc-cutadapt-star-umitools/`, built from the
enabled workflow steps + tools so different tool combinations never collide.

## Per-workflow-step blocks

Each step listed in `WORKFLOWS` needs:

```json
"MAPPING": {
  "TOOLS": { "star": "STAR" },
  "LabA": {
    "wild-type": {
      "ENV": "star",
      "BIN": "STAR",
      "star": {
        "OPTIONS": {
          "INDEX": "--some-star-index-flag value",
          "MAP": "--some-star-mapping-flag value"
        }
      }
    }
  }
}
```

- `TOOLS`: map of `conda-env-name: executable-name` for every tool used in
  this step (env name usually equals the bioconda package, see `envs/`).
- `ENV`/`BIN` at the condition level override the step-level default (useful
  to benchmark two installs of the same tool, or point at a local binary for
  tools without a conda package, e.g. `guppy`, `dorado`).
- Inside each tool's block, `OPTIONS` holds one key per sub-step (e.g.
  `INDEX` vs `MAP` for a mapper) with a **raw command-line string** of extra
  flags for that tool — exactly as you'd type them, minus anything MONSDA
  derives automatically (paired/stranded flags, threads). Leave `OPTIONS: {}`
  if no extra flags are needed; never invent a flag that isn't in the tool's
  own `--help`.
- `COMPARABLE` (DE/DEU/DAS/DTU and MACS peak calling): if omitted/empty,
  MONSDA generates all-vs-all comparisons from `GROUPS`. Otherwise:
  ```json
  "COMPARABLE": { "comparison-name": ["Group1", "Group2"] }
  ```

## Fastest path: samplesheet instead of hand-written SETTINGS

If the user can provide a CSV/TSV, skip writing `SETTINGS` by hand. Columns
(see `configs/samplesheet_template.csv`):

```
CONDITION,SAMPLE,GROUP,SEQUENCING,REFERENCE,GTF,GFF,INDEX,PREFIX,DECOY,TYPE,BATCH,IP
```

`CONDITION` is a `/`-separated path that becomes the condition-tree branch
(e.g. `FGUMI/WT/dummylevel`). Rows sharing a `CONDITION` are grouped into one
leaf's `SAMPLES` list; only the first row of a group needs the
REFERENCE/GTF/GFF/etc. columns filled in. Pass it to `monsda` directly with
`--samplesheet file.csv` (it populates `SETTINGS` in the given config JSON on
the first run and writes `<config>_with_settings.json` for reuse), or let
`monsda_configure --samplesheet file.csv` do it when scaffolding a brand-new
project.

## PostDE (optional enrichment on top of DE)

Only for the `DE` workflow with `deseq2`/`edger` (not DEU/DAS/DTU). Enabled
via a top-level `POSTDE` key pointing at a separate analysis JSON (gProfiler/
clusterProfiler/GSVA/decoupleR/dream). See `docs/source/config.rst` section
"PostDE" and the shipped example `configs/postde_analysis.json` before
building one — it has its own schema and conda env (`postde`).

## Starting points in this repo

Valid JSON, safe to copy/parse directly:

- `configs/template_clean.json` — empty skeleton with every workflow step's
  key already scaffolded.
- `configs/tutorial_quick.json` — small worked example (FETCH+MAPPING on the
  bundled `tests/data` Ecoli set).
- `configs/template_umicollapse.json` — UMI-aware dedup example.
- `configs/postde_analysis.json` — example PostDE analysis JSON.

**Documentation-only, NOT valid JSON** — these have trailing `#comment`
annotations after values and will fail `json.load`; read them for the
explanations but copy structure from the valid files above instead:
`configs/template.json`, `configs/template_base_commented.json`,
`configs/tutorial_toolmix.json`, `configs/tutorial_exhaustive.json`,
`configs/tutorial_postprocess.json`. Strip the `#...` trailing comments
before using one of these as a base, or diff it against
`configs/template_clean.json` to see the annotation-free version of the same
keys.
