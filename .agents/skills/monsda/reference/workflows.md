# Workflow steps and tools

Full tables with descriptions/links: `docs/source/workflows.rst`. Condensed
reference below: `WORKFLOWS` step name -> `conda-env: executable` pairs you
can put under that step's `TOOLS` key. Env names match a `.yaml` in `envs/`.

| Step | Purpose | `TOOLS` env : bin examples |
|---|---|---|
| `FETCH` | Download from SRA | `sra: fasterq-dump` |
| `BASECALL` | FAST5/POD5 -> FASTQ | `guppy: $PATH_TO_LOCAL_BIN` (no conda pkg), `dorado: $PATH_TO_LOCAL_BIN` (no conda pkg) |
| `QC` | QC of FASTQ/BAM (runs standalone pre-processing, or again post-processing for any later step) | `fastqc: fastqc` (FASTQ/BAM, includes MultiQC), `rustqc: rustqc` (BAM only, includes MultiQC, has a strandedness guard — aborts the run if inferred strandedness conflicts with `SEQUENCING`) |
| `TRIMMING` | Adapter/quality trimming | `trimgalore: trim_galore`, `cutadapt: cutadapt`, `bbduk: bbduk`, `fastp: fastp` |
| `MAPPING` | Align to reference | `hisat2: hisat2`, `star: STAR` (also used for `STARsolo` single-cell mode), `rustar: rustar` (Rust STAR reimplementation), `segemehl2`/`segemehl3`/`segemehl2bisulfite`/`segemehl3bisulfite: segemehl.x`, `bwa: bwa mem`, `bwa2: bwa2-mem`, `bwameth: bwameth.py`, `bowtie2: bowtie2` (no native stranded-protocol flag, so the strandedness part of `SEQUENCING` is not forwarded to it), `minimap: minimap2` (long-read, use `-x map-ont`/`-x map-pb` in OPTIONS), `salmonalign: salmon` (alignment mode), `rammap: rammap` (long-read, not on bioconda), `piscem: piscem` (single-cell only, emits RAD for `alevinfry`) |
| `DEDUP` | Remove/mark duplicates | `umitools: umi_tools`, `fgumi: fgumi` (UMI extraction + dedup), `picard: picard` (position-based) |
| `COUNTING` | Quantify reads/transcripts | `countreads: featureCounts`, `salmon: salmon` (FASTQ or, via `OPTIONS.BAM = "transcriptome"\|"genome"`, directly from BAM), `kallisto: kallisto`, `oarfish: oarfish` (long-read, from transcriptome BAM), `alevinfry: alevin-fry` (needs a piscem RAD file), `simpleaf: simpleaf` (self-contained single-cell index+map+quant) |
| `TRACKS` | trackdb.txt + BIGWIG for UCSC/genome browsers | `ucsc: ucsc` |
| `PEAKS` | Peak calling (ChIP/RIP/CLIP) | `piranha: piranha` (CLIP/RIP), `macs: macs` (ChIP TF binding; use `COMPARABLE` for signal-vs-background pairs), `scyphy: piranha` (cyPhyRNA-Seq), `peaks: peaks` (quick/unpublished window scan) |
| `DE` | Differential gene expression | `deseq2`, `edger` (R scripts under `Analysis/DE/`; both accept a `COUNTTOOL` option selecting the counter: `featureCounts` for short reads, default, or `oarfish` for long-read BAMs — oarfish quantifies from the mapped genome BAM and transcript counts are aggregated to gene level; tool options come from `OPTIONS.QUANT`, a legacy `OPTIONS.COUNT` is still accepted as a fallback and `QUANT` wins on conflict) |
| `DEU` | Differential exon usage | `dexseq`, `edger` (both accept a `COUNTTOOL` option; only `featureCounts`, the default, is currently implemented) |
| `DAS` | Differential alternative splicing | `edger`, `diego` (both accept a `COUNTTOOL` option; only `featureCounts`, the default, is currently implemented) |
| `DTU` | Differential transcript usage | `dexseq`, `drimseq`, `spit` (all three accept an optional `TERMINUS` option to pre-group transcripts, and a `COUNTTOOL` option selecting the quantifier: `salmon` for short reads, default, or `oarfish` for long-read BAMs; `TERMINUS` only works with `COUNTTOOL: salmon`) |
| `CIRCS` | circRNA detection | `ciri2: $Path_to_CIRI2.pl` (install locally, no conda pkg) |
| `FUSIONS` | Fusion gene detection | `starfusion: STAR-Fusion` — needs `MAPPING` with `star` and chimeric output enabled (`--chimSegmentMin 12 --chimOutType Junctions` in the STAR `MAP` OPTIONS), or set `FUSIONS.OPTIONS.FASTQ = "true"` to run directly on trimmed FASTQ instead. Needs a CTAT genome resource dir in `INDEX` (auto-built from `REFERENCE`+`ANNOTATION` if missing — slow). |

## BASECALL raw signal (`RAW/<condition-path>`)

One basecalling job is generated per biological sample and every raw signal
file of that sample is a tracked input. The selected caller fixes the format:
`guppy` reads `.fast5`, `dorado` reads `.pod5`; files of the other format are
ignored and never mixed into a sample. Each sample owns its raw signal through
exactly one of three sources below `RAW/<condition-path>` (the condition path
of its leaf, e.g. `RAW/Test` for `Test`):

- the single file `<id>.fast5` (`<id>.pod5` for dorado)
- the directory `<id>/`, scanned recursively for the caller's extension
- the explicit list `<id>.fast5.list` (`<id>.pod5.list`), one path per line,
  relative to `RAW/<condition-path>`, for loose files outside a per-sample
  directory

```text
RAW/Test/S1.fast5              # single file, guppy
RAW/Test/S2/run_1.fast5        # recursive directory, guppy
RAW/Test/S3.fast5.list         # explicit list, guppy
RAW/Test/loose/loose_a.fast5
```

with `S3.fast5.list` holding `loose/loose_a.fast5`.

Matching is exact, there is no prefix guessing: `S1` never picks up
`S10.fast5` or `S1_extra.fast5`. MONSDA stops with an actionable error when a
sample has no source, when more than one source exists (ambiguous), when a
listed path is absolute, leaves the condition directory, is missing, is not a
regular file, has the wrong extension or is empty, when a list repeats an
entry, and when one raw file would be claimed by two samples. Output is one
FASTQ per sample, `FASTQ/<condition-path>/<id>.fastq.gz`, in both the
Snakemake and the Nextflow engine.

Notes:

- QC runs automatically again after any processing step if enabled, in
  addition to the optional pre-processing QC pass.
- Single-cell: `piscem` (MAPPING) + `alevinfry` (COUNTING), or the
  self-contained `simpleaf` (COUNTING only); generic post-mapping QC/DEDUP is
  skipped for the RAD-file path.
- Long-read (Nanopore/PacBio): `minimap`/`rammap` (MAPPING) + `oarfish`
  (COUNTING, from transcriptome BAM). `DTU` can also run on long reads by
  setting `COUNTTOOL: oarfish` in each DTU tool's `OPTIONS` (it reuses the
  already-mapped genome BAM instead of a salmon index/quant step). `DE`
  works the same way via `COUNTTOOL: oarfish` in each DE tool's `OPTIONS`
  (oarfish quantifies from the genome BAM, transcript counts are summed to
  gene level into a featureCounts-compatible table).
- Bisulfite: `bwameth` or `segemehl{2,3}bisulfite` for MAPPING.
- Tools marked "no conda pkg" need `ENV`/`BIN` pointing at a locally
  installed executable instead of a bioconda environment name.
- featureCounts count units: the featureCounts calls in `COUNTING`, `DE`,
  `DEU` and `DAS` run with `-p --countReadPairs` when the resolved
  pairedness of the step is `paired`, so counts are fragments (read pairs);
  single-end data is counted per read. No config key and no extra `OPTIONS`
  flag is needed.
- `PostDE` (gProfiler/clusterProfiler/GSVA/decoupleR/dream enrichment on top
  of `DE` with `deseq2`/`edger`) lives under `DE.POSTDE` in the main config,
  not a `WORKFLOWS` entry — see `reference/config-guide.md`. With a
  background configured every analysis runs twice (full gene set +
  background-restricted), outputs in `<entry>/full/` and
  `<entry>/background/`.
