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
| `MAPPING` | Align to reference | `hisat2: hisat2`, `star: STAR` (also used for `STARsolo` single-cell mode), `rustar: rustar` (Rust STAR reimplementation), `segemehl2`/`segemehl3`/`segemehl2bisulfite`/`segemehl3bisulfite: segemehl.x`, `bwa: bwa mem`, `bwa2: bwa2-mem`, `bwameth: bwameth.py`, `minimap: minimap2` (long-read, use `-x map-ont`/`-x map-pb` in OPTIONS), `salmonalign: salmon` (alignment mode), `rammap: rammap` (long-read, not on bioconda), `piscem: piscem` (single-cell only, emits RAD for `alevinfry`) |
| `DEDUP` | Remove/mark duplicates | `umitools: umi_tools`, `fgumi: fgumi` (UMI extraction + dedup), `picard: picard` (position-based) |
| `COUNTING` | Quantify reads/transcripts | `countreads: featureCounts`, `salmon: salmon` (FASTQ or, via `OPTIONS.BAM = "transcriptome"\|"genome"`, directly from BAM), `kallisto: kallisto`, `oarfish: oarfish` (long-read, from transcriptome BAM), `alevinfry: alevin-fry` (needs a piscem RAD file), `simpleaf: simpleaf` (self-contained single-cell index+map+quant) |
| `TRACKS` | trackdb.txt + BIGWIG for UCSC/genome browsers | `ucsc: ucsc` |
| `PEAKS` | Peak calling (ChIP/RIP/CLIP) | `piranha: piranha` (CLIP/RIP), `macs: macs` (ChIP TF binding; use `COMPARABLE` for signal-vs-background pairs), `scyphy: piranha` (cyPhyRNA-Seq), `peaks: peaks` (quick/unpublished window scan) |
| `DE` | Differential gene expression | `deseq2`, `edger` (R scripts under `Analysis/DE/`) |
| `DEU` | Differential exon usage | `dexseq`, `edger` |
| `DAS` | Differential alternative splicing | `edger`, `diego` |
| `DTU` | Differential transcript usage | `dexseq`, `drimseq`, `spit` (all three accept an optional `TERMINUS` option to pre-group transcripts) |
| `CIRCS` | circRNA detection | `ciri2: $Path_to_CIRI2.pl` (install locally, no conda pkg) |
| `FUSIONS` | Fusion gene detection | `starfusion: STAR-Fusion` — needs `MAPPING` with `star` and chimeric output enabled (`--chimSegmentMin 12 --chimOutType Junctions` in the STAR `MAP` OPTIONS), or set `FUSIONS.OPTIONS.FASTQ = "true"` to run directly on trimmed FASTQ instead. Needs a CTAT genome resource dir in `INDEX` (auto-built from `REFERENCE`+`ANNOTATION` if missing — slow). |

Notes:

- QC runs automatically again after any processing step if enabled, in
  addition to the optional pre-processing QC pass.
- Single-cell: `piscem` (MAPPING) + `alevinfry` (COUNTING), or the
  self-contained `simpleaf` (COUNTING only); generic post-mapping QC/DEDUP is
  skipped for the RAD-file path.
- Long-read (Nanopore/PacBio): `minimap`/`rammap` (MAPPING) + `oarfish`
  (COUNTING, from transcriptome BAM).
- Bisulfite: `bwameth` or `segemehl{2,3}bisulfite` for MAPPING.
- Tools marked "no conda pkg" need `ENV`/`BIN` pointing at a locally
  installed executable instead of a bioconda environment name.
- `PostDE` (gProfiler/clusterProfiler/GSVA/decoupleR/dream enrichment on top
  of `DE` with `deseq2`/`edger`) is a separate top-level `POSTDE` config, not
  a `WORKFLOWS` entry — see `reference/config-guide.md`.
