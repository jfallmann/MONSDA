# config.json anatomy

Full human docs: `docs/source/config.rst`, `docs/source/conditiontree.rst`,
`docs/source/preparation.rst`. This is the condensed, agent-facing version.

## Top level

```json
{
  "WORKFLOWS": "QC,TRIMMING,MAPPING,DEDUP",
  "BINS": "",
  "MAXTHREADS": "20",
  "VERSION": "1.6.0",
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
- `SAMPLES`: one entry per sample, no file extension. What the entry means
  depends on `SEQUENCING`:
  - `paired` and `singlecell`: the bare `<id>` of a pair
    `<id>_R1.fastq.gz`/`<id>_R2.fastq.gz`. MONSDA derives the `_R1`/`_R2`
    suffixes itself, so read tags never appear in the config.
  - `single`: the file stem, i.e. the name without `.fastq.gz`
    (`<id>.fastq.gz`). If that stem itself ends in a literal `_R1`
    (`sample_R1.fastq.gz`), the `_R1` is part of the sample name and is kept.
  - Lane-split files (`<id>_L001_R1.fastq.gz`, and the older reversed
    `<id>_R1_L001.fastq.gz`) are concatenated into the canonical
    `<id>_R1.fastq.gz`/`<id>_R2.fastq.gz` for `paired`/`singlecell` and into
    `<id>.fastq.gz` for `single`, so no manual renaming is needed.
  - Illumina BCL Convert writes `<Sample_ID>_S#_L00#_R#_001.fastq.gz`, or
    without the lane part under `--no-lane-splitting`. One entry covers all
    of it: `381915_S1_L001_R1_001.fastq.gz` and
    `381915_S1_L002_R2_001.fastq.gz` are both sample `381915_S1`. Keep the
    `S#` sample number as part of the name; the trailing `_001` is a file
    counter, not a sample or replicate.
  - Identifiers are not validated, but keep the metadata tokens `_L<digits>`,
    `_R1`/`_R2` and `_001` out of a biological label (`Tumor-L001`, not
    `Tumor_L001`), since they cannot be told apart from file-name syntax.
    Keep identifiers unique across the project; if the same basename is
    reused in different condition directories, the condition position in the
    tree qualifies it, so do not type the condition into `SAMPLES`.
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
- `SEQUENCING` can likewise be set at the step level (or nested inside a
  tool's own block) to override the global `SETTINGS.SEQUENCING` for just
  that step, same `"single"/"paired"/"singlecell"[,rf|fr]` syntax. Useful for
  e.g. QuantSeq-style 3' protocols: keep the global value `"paired,rf"` so
  TRIMMING/DEDUP see both reads (read2 carries the UMI), but set
  `"SEQUENCING": "single,rf"` under `MAPPING` (and `COUNTING` if needed) so
  only read1 is actually mapped/counted. Steps downstream of MAPPING (PEAKS,
  TRACKS, COUNTING, CIRCS, FUSIONS, DE, DEU, DAS) default to the resolved
  MAPPING pairedness/strandedness, not the global one, when they don't set
  their own — a MAPPING-level override already propagates to DE/DEU/DAS
  automatically, e.g. so their featureCounts calls use the right `-p`/`-s`
  flags. DTU is the exception: its default salmon path quantifies from the
  trimmed FASTQ directly, independent of MAPPING, so it falls back to the
  global value instead.
- Inside each tool's block, `OPTIONS` holds one key per sub-step (e.g.
  `INDEX` vs `MAP` for a mapper) with a **raw command-line string** of extra
  flags for that tool — exactly as you'd type them, minus anything MONSDA
  derives automatically (paired/stranded flags, threads). Leave `OPTIONS: {}`
  if no extra flags are needed; never invent a flag that isn't in the tool's
  own `--help`.
- `OPTIONS.QUANT` is the documented key for the quantification flags of
  `COUNTING`, `DE`, `DEU`, `DAS` and `DTU`, in both the Snakemake and the
  Nextflow engine. If `QUANT` is missing or blank, MONSDA falls back to a
  non-blank legacy `OPTIONS.COUNT` and logs a warning; the two are never
  merged, so a set `QUANT` simply wins and `COUNT` is dropped. `COUNTTOOL`
  (which counting tool to run, e.g. `featureCounts`/`salmon`) is a separate
  key and is not affected.
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

Only for the `DE` workflow with `deseq2`/`edger` (not DEU/DAS/DTU). Configured
in the main config under `"DE": {"POSTDE": {...}}` — one block, no separate
JSON file (methods: gProfiler/clusterProfiler/GSVA/decoupleR/dream; top-level
`POSTDE` and a `config` key are rejected). It has its own conda env
(`postde`). Defaults live in `MONSDA/extras/PostDE.py::POSTDE_DEFAULTS` and in
`configs/template_base_commented.json`; the summary keys `plots`, `report`
and `shiny` are rejected until a reporting module ships.

Key semantics (validated at startup, before any workflow is generated):

- One global block next to the DE tool options, applied to every contrast.
- `background` (optional, single block replacing per-method background keys):
  - `file`: one gene id per line, relative to the **run directory** (not the
    config file's directory); staged into the PostDE workdir as
    `background_file.<ext>`.
  - `expressed`: derive from the DE run — `metric` `tpm|cpm|expression`,
    `min` strict `>` cutoff, `min_samples` (default 1), `per_group`
    (default true). `expression` = the engine's own average-expression
    column: DESeq2 `baseMean`, edgeR `logCPM` — works with both tools but
    `min` is on the engine's scale (baseMean raw counts vs logCPM log2);
    no `min_samples`/`per_group`. `tpm` uses the full gene span from the
    GTF, so absolute cutoffs are approximate; `cpm` is the robust default.
    Counts have spike-ins removed and are raw (not normalized).
  - `combine`: `intersect`|`union` of `file` and `expressed`; the result is
    also intersected with the tested gene set.
  - `include_significant` (default true): significant genes below the
    cutoff are added back and a warning lists them.
- With a background set, **every** analysis (gprofiler, clusterProfiler/GSEA,
  GSVA, decoupleR, dream) runs on the full tested gene set AND on the
  background-restricted set: outputs land in `<entry>/full/` and
  `<entry>/background/` (plus `background_summary.json`), so results can be
  compared against the ground truth. Without a background, outputs land in
  `<entry>/` directly. Runtime roughly doubles with a background.
- Restricted gProfiler runs map `domain_scope` `annotated`→
  `custom_annotated` and `known`→`custom` automatically.

## Starting points in this repo

Valid JSON, safe to copy/parse directly:

- `configs/template_clean.json` — empty skeleton with every workflow step's
  key already scaffolded.
- `configs/tutorial_quick.json` — small worked example (FETCH+MAPPING on the
  bundled `tests/data` Ecoli set).
- `configs/template_umicollapse.json` — UMI-aware dedup example.

**Documentation-only, NOT valid JSON** — these have trailing `#comment`
annotations after values and will fail `json.load`; read them for the
explanations but copy structure from the valid files above instead:
`configs/template.json`, `configs/template_base_commented.json`,
`configs/tutorial_toolmix.json`, `configs/tutorial_exhaustive.json`,
`configs/tutorial_postprocess.json`. Strip the `#...` trailing comments
before using one of these as a base, or diff it against
`configs/template_clean.json` to see the annotation-free version of the same
keys.
