# Tool option audit

Purpose: cross-check hardcoded CLI flags in workflows/*.smk and *.nf against the
actual `--help` output of each tool inside its container in
~/Work2/Container/MONSDA (apptainer conda env), to catch flag drift like the
oarfish `--reference` -> `--genome-fasta` rename fixed in workflows/oarfish.smk
and workflows/oarfish.nf (commit pending on dev_1.5.0).

Status: IN PROGRESS. This file is updated incrementally as tools are checked.

## Automated audit

The manual process below has been turned into a repeatable script,
`audit_tool_flags.py`, driven by the declarative `tool_manifest.yaml`
(one entry per tool, transcribing every finding in this file). Run it
whenever containers or tool versions are updated:

    source ~/anaconda3/etc/profile.d/conda.sh && conda activate apptainer
    python tests/tool_audit/audit_tool_flags.py            # default: all tools
    python tests/tool_audit/audit_tool_flags.py --tool oarfish --tool salmon
    python tests/tool_audit/audit_tool_flags.py --strict    # also fail known gaps
    python tests/tool_audit/audit_tool_flags.py --json-out /tmp/report.json

It classifies every check as PASS, DRIFT (a hardcoded flag no longer
matches the tool's own `--help`/behaviour - a real regression, always
fails the run), GAP (container or binary missing - informational if
already listed with `known_gap: true`, otherwise a hard failure), or SKIP
(explicitly not applicable, e.g. Piranha's config-only params). Exit code
is non-zero on any DRIFT or newly-appeared GAP.

It is also collected by `pytest` via `tests/test_tool_audit.py`, which
self-skips with a visible reason on any machine without `apptainer` and
the container directory (mirrors the `multiqc`-not-installed note in
AGENTS.md - this cannot run in a bare CI checkout).

`tool_manifest.yaml` is the reviewable source of truth and must be updated
by hand whenever a workflow's hardcoded flags change. Use
`python tests/tool_audit/extract_workflow_commands.py --refresh` to
regenerate `smk_shell_commands.txt`/`nf_script_commands.txt` as a reference
while doing so, and `--check-manifest` to sanity-check that every
`source_refs` entry in the manifest is still found in those dumps (catches
the extractor itself silently breaking on new syntax, not shell semantics).
The script never edits this README automatically - promote genuinely new
findings into the prose sections below by hand after fixing them.

## Method (manual, superseded by the script above for routine re-checks)
    source ~/anaconda3/etc/profile.d/conda.sh && conda activate apptainer
    apptainer exec ~/Work2/Container/MONSDA/<tool>-1.5.0.sif <binary> --help

## Findings

### oarfish (fixed, commit 2f58130)
`--reference` renamed to `--genome-fasta` in the current oarfish binary
(container oarfish-1.5.0.sif). `--genome-alignments` mode error was:
`error: unexpected argument '--reference' found`.
Fixed in workflows/oarfish.smk and workflows/oarfish.nf.

### salmonalign (fixed)
salmon in the container (salmon-1.5.0.sif) is salmon 2.5.1, a full Rust
rewrite. `--writeMappings` now REQUIRES an explicit filename value
(`-z, --writeMappings <WRITE_MAPPINGS>`), unlike the old C++ salmon where
`--writeMappings` with no argument defaulted to writing SAM to stdout.
workflows/salmonalign.smk and .nf relied on that old stdout-default behavior
(`--writeMappings -o {outdir} ... | tee >(samtools ...)`), so the shell
piped nothing into `tee` and the whole aligned/unmapped split silently
produced empty streams (confirmed by reproducing with a minimal salmon index
+ single read: without a value clap errors
`a value is required for '--writeMappings <WRITE_MAPPINGS>' but none was
supplied`).
Fix: use `--writeMappings=/dev/stdout`, confirmed to write valid SAM to
stdout in salmon 2.5.1. Applied in workflows/salmonalign.smk (both PE/SE
rules) and workflows/salmonalign.nf (both PE/SE processes).

### salmon.smk / salmon_trim.smk (checked, no change needed)
`-d/--decoy` for `salmon index`, `-a/--alignments`, `--annotation`,
`--genome`, `--deterministic` for `salmon quant` (genome/transcriptome
counting modes used by workflows/salmon.smk) all exist as-is in salmon
2.5.1 `salmon quant --help`. No `--writeMappings` used there since it does
not stream SAM to stdout (counts-only), so unaffected by the above.

### kallisto (checked, no change needed)
kallisto 0.51.0 in container. `-d/--d-list` (decoy), `-i`, `-t`, `-o`,
`--single`, `-1/-2`/positional fastqs used in workflows/kallisto.smk and
kallisto_trim.smk all match `kallisto index --help` / `kallisto quant --help`.

### bwameth (checked, no change needed)
`--reference` is the correct, still-current flag for bwameth.py (confirmed
via `bwameth.py --help`), unlike oarfish/salmon. Not a false positive.

### fgumi (checked, no change needed)
Container binary is `fgumi` (not `fgbio`). `extract --inputs/--sample
/--library/--output` all present in `fgumi extract --help`; matches
workflows/fgumi.smk.

### Tools whose containers do NOT actually contain the expected binary
(cannot verify CLI flags locally, flagged as a gap, not necessarily a bug):
- scyphy-1.5.0.sif: no `scyphy` binary in /opt/conda/envs/monsda/bin
- dorado-1.5.0.sif: no `dorado` binary
- guppy-1.5.0.sif: no `guppy_basecaller` binary (env only has gettext tools)
- ciri2-1.5.0.sif: no `CIRI2.pl` binary
- rustar-1.5.0.sif: no rustar-named binary found (has rust toolchain only)
These may be intentionally thin/build-in-progress containers (guppy/dorado
often require separate ONT-licensed installers), or the .sif build failed to
install the tool - not confirmed which. Not re-tested further this session.

### piscem / alevinfry (checked, no change needed)
`piscem build` and `piscem map-sc` help match workflows/piscem.smk exactly
(`-o`, `--decoy-paths`, `-k`, `-t`; `map-sc -t -i -1 -2 -o`, `--geometry` is a
required user config value, not hardcoded). `alevin-fry generate-permit-list
/collate/quant` flags (`-i`, `-d`, `-o`, `-t`, `-r`, `-m`, `--use-mtx`) all
match workflows/alevinfry.smk.

### star / hisat2 (checked, no change needed)
STAR help (`STAR --help`) confirms `--runThreadN`, `--runMode
genomeGenerate`, `--genomeDir`, `--genomeFastaFiles`, `--outFileNamePrefix`,
`--outTmpDir`, `--readFilesCommand`, `--readFilesIn`, `--outReadsUnmapped
Fastx`, `--outSAMtype BAM SortedByCoordinate`, and every tag in
`--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM` used in
workflows/star.smk. hisat2 `--un-conc-gz`/`-p`/`--summary-file` etc. all
present in `hisat2 --help`.

### macs2 / macs3 (checked, no change needed)
`callpeak -t -c --outdir -n -f` all present in both `macs2 callpeak --help`
and `macs3 callpeak --help`, matching workflows/macs2.smk and macs3.smk.

### umicollapse (checked, no change needed)
workflows/umicollapse.smk uses two different binaries: `umi_tools` for the
`whitelist`/`extract` rules (flags `--stdin/--stdout/--read2-in/--read2-out
/--error-correct-cell/--whitelist/--temp-dir/--log` all confirmed present in
`umi_tools whitelist --help` / `umi_tools extract --help`), and the Java
`umicollapse` tool only for the final `bam` dedup rule. Ran `umicollapse bam
-i tiny.bam -o out.bam --paired` against a real tiny BAM in the container
(umicollapse 1.1.0-2) and confirmed `-i/-o/--paired` work and produce output.

### umitools (checked, no change needed)
`umi_tools dedup --stdin/--stdout/--paired/--log/--temp-dir` all present in
`umi_tools dedup --help`, matching workflows/umitools.smk.

### rustqc (checked, no change needed)
`rustqc rna --gtf/-g -t -s -p -o -j --skip-dup-check` etc. all present in
`rustqc rna --help`, matching workflows/rustqc.smk.

### trimgalore / cutadapt / fastp / bbduk (checked, no change needed)
All hardcoded flags (`--cores/--paired/--gzip/-o`; `--cores -o -p`;
`--thread --in1/--in2/--out1/--out2`; `t=/in1=/in2=/out1=/out2=`) confirmed
present in each tool's own `--help`.

### segemehl / segemehl3 (checked, no change needed, naming note)
Binary is `segemehl.x` in both the segemehl-1.5.0.sif and segemehl3-1.5.0.sif
containers (segemehl3 does not ship a separate `segemehl3` executable). Not a
bug: workflows/segemehl3.smk resolves the binary via config (`MAPPERBIN`),
never hardcodes the name `segemehl3`. `-d/-x/-q/-p/--threads` all present in
`segemehl.x --help`.

### minimap2 / sra (checked, no change needed)
`minimap2 -d/-t` and `fasterq-dump -O/-e/-t/--split-files` confirmed present
in their own `--help`, matching workflows/minimap.smk and sra.smk.

### picard (checked, no change needed)
`MarkDuplicates --REMOVE_DUPLICATES/--ASSUME_SORT_ORDER/--TMP_DIR/--INPUT
/--OUTPUT/--METRICS_FILE` all confirmed present in `picard MarkDuplicates
--help`, matching workflows/picard_dedup.smk.

### piranha (not applicable)
workflows/piranha.smk never hardcodes Piranha CLI flags itself - the actual
peak-calling parameters come entirely from user config (`{params.ppara}`).
`Piranha --help` requires a positional input and errored differently when
called with only `--help`, but this does not affect MONSDA since no fixed
flags are baked into the template.

## Summary
Two real bugs found and fixed this session (both were caused by the
containers having been rebuilt against newer tool versions with renamed/
stricter CLI flags):
1. `oarfish`: `--reference` -> `--genome-fasta` (workflows/oarfish.smk,
   workflows/oarfish.nf).
2. `salmon` (salmonalign only): `--writeMappings` now requires an explicit
   path argument in salmon 2.5.1; changed to `--writeMappings=/dev/stdout`
   (workflows/salmonalign.smk, workflows/salmonalign.nf).

Every other tool with hardcoded flags in workflows/*.smk (piscem, alevinfry,
star, hisat2, macs2, macs3, umicollapse, umitools, rustqc, trimgalore,
cutadapt, fastp, bbduk, kallisto, bwameth, salmon counting modes, segemehl,
minimap2, sra, picard, fgumi) was checked against the actual container's
`--help` (or, where `--help` errors before printing usage, a real minimal
run) and found consistent - no other hardcoded-flag drift detected.

Not checked this session (no hardcoded external-tool CLI flags were found in
their .smk/.nf, or they are R/perl scripts calling internal MONSDA scripts
rather than external tool binaries): dexseq_DEU/DTU, edger_DAS/DE/DEU,
deseq2_DE, drimseq_DTU, diego_DAS, spit_DTU, annotatebed, ucsc, ciri2,
scyphy, dorado, guppy, rustar, multiqc/multiqc_rustqc, countreads,
manipulate_genome, bwa/bwa2.

### Tools whose containers do NOT actually contain the expected binary
(cannot verify CLI flags locally, flagged as a gap, not necessarily a bug):
- scyphy-1.5.0.sif: no `scyphy` binary in /opt/conda/envs/monsda/bin
- dorado-1.5.0.sif: no `dorado` binary
- guppy-1.5.0.sif: no `guppy_basecaller` binary (env only has gettext tools)
- ciri2-1.5.0.sif: no `CIRI2.pl` binary
- rustar-1.5.0.sif: no rustar-named binary found (has rust toolchain only)
These may be intentionally thin/build-in-progress containers (guppy/dorado
often require separate ONT-licensed installers), or the .sif build failed to
install the tool - not confirmed which. Not re-tested further this session.

