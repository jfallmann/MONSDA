#!/usr/bin/env python3
"""Create a tiny paired-end FASTQ fixture with synthetic UMIs for fgumi tests.

This script writes gzipped FASTQ files into:
  FASTQ/Test/umi/FGUMI01_R1.fastq.gz
  FASTQ/Test/umi/FGUMI01_R2.fastq.gz

R1 layout is intentionally compatible with a simple fgumi read structure like:
  --read-structures 6M+T +T
where the first 6 bases in R1 are a mock UMI.
"""

from __future__ import annotations

import gzip
from pathlib import Path

OUTDIR = Path("FASTQ") / "Test" / "umi"
SAMPLE = "FGUMI01"

# (umi6, r1_suffix, r2_sequence)
READS = [
    ("ACGTAA", "TTTTTTTTTTTTAACCGGTTCCAA", "GATTACAGATTACAGATTACAGATTACA"),
    ("ACGTAA", "TTTTTTTTTTTTAACCGGTTCCAA", "GATTACAGATTACAGATTACAGATTACA"),  # duplicate
    ("TGCACT", "TTTTTTTTTTTTGGCCAATTGGCC", "CGTACGTACGTACGTACGTACGTACGTAC"),
    (
        "TGCACT",
        "TTTTTTTTTTTTGGCCAATTGGCC",
        "CGTACGTACGTACGTACGTACGTACGTAC",
    ),  # duplicate
    ("GGAACC", "TTTTTTTTTTTTCCGGAATTCCGG", "TTGCAATTGCAATTGCAATTGCAATTGCA"),
    ("CCTTGG", "TTTTTTTTTTTTAATTCCGGAATT", "AACCGGTTAACCGGTTAACCGGTTAACCG"),
    ("TTAAGG", "TTTTTTTTTTTTGGAATTCCGGAA", "GGCCAATTGGCCAATTGGCCAATTGGCCA"),
    ("CCGGTT", "TTTTTTTTTTTTAACCTTGGAACC", "ATATCGCGATATCGCGATATCGCGATATC"),
]


def _fq_record(name: str, seq: str, qual_char: str = "I") -> str:
    return f"@{name}\n{seq}\n+\n{qual_char * len(seq)}\n"


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    r1_path = OUTDIR / f"{SAMPLE}_R1.fastq.gz"
    r2_path = OUTDIR / f"{SAMPLE}_R2.fastq.gz"

    with gzip.open(r1_path, "wt") as r1h, gzip.open(r2_path, "wt") as r2h:
        for i, (umi, r1_suffix, r2_seq) in enumerate(READS, start=1):
            read_id = f"{SAMPLE}:{i}"
            r1_seq = umi + r1_suffix
            r1h.write(_fq_record(read_id + " 1:N:0:TEST", r1_seq))
            r2h.write(_fq_record(read_id + " 2:N:0:TEST", r2_seq))

    print(f"Wrote {r1_path}")
    print(f"Wrote {r2_path}")


if __name__ == "__main__":
    main()
