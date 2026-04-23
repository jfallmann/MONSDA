# fgumi UMI smoke-test data

This directory contains a tiny synthetic paired-end fixture generator for `fgumi`:

- `make_fgumi_umi_fixture.py` writes
  - `FASTQ/Test/umi/FGUMI01_R1.fastq.gz`
  - `FASTQ/Test/umi/FGUMI01_R2.fastq.gz`
- `config_fgumi_test.json` is a minimal tutorial-style MONSDA config that enables only `DEDUP` with `fgumi`.

## Why synthetic data?

Most publicly referenced UMI tutorial datasets are not truly small enough for quick CI-style smoke tests.
The synthetic fixture keeps runtime tiny and deterministic while still exercising UMI extraction behavior.

## External UMI datasets (if you want real data)

- Galaxy Training (CEL-Seq2 UMI tutorial data on Zenodo):
  - https://zenodo.org/record/2573177
  - Files:
    - `test_barcodes_celseq2_R1.fastq.gz` (~243.7 MB)
    - `test_barcodes_celseq2_R2.fastq.gz` (~594.9 MB)
- Galaxy tutorial page:
  - https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-umis/tutorial.html

For nf-core test datasets, use branch discovery/search tooling first:
- https://github.com/nf-core/test-datasets
- https://github.com/nf-core/test-datasets/blob/master/docs/USE_EXISTING_DATA.md
