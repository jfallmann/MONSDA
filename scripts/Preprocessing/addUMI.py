#!/usr/bin/env python3
"""
Prepend a UMI to every read in one or more FASTQ files.

UMIs can be:
  - randomly generated (default, see --umi-length / --qual-char)
  - extracted from a paired UMI FASTQ file (--umi-fastq)

When --umi-fastq is given the UMI sequence *and* its base qualities are taken
from the corresponding read in that file.  Read pairing is validated by dnaio,
which compares read names at every position (tolerating /1 /2 suffixes and
trailing descriptions).  A mismatch raises an error immediately instead of
silently producing wrong results.

Output files are written next to the originals with an "_umi" suffix,
e.g. input_umi.fastq.gz.  Gzip I/O uses igzip/ISA-L when available.

Usage (random UMIs):
    python addUMI.py [--umi-length 12] [--qual-char F] input.fastq[.gz] ...

Usage (UMIs from file):
    python addUMI.py --umi-fastq umi.fastq[.gz] input.fastq[.gz]

Requires:
    dnaio >= 0.9   (pip install dnaio)
"""

import argparse
import os
import sys
from pathlib import Path

try:
    import dnaio
except ImportError:
    sys.exit("dnaio is required: pip install dnaio")

NUC = b"ACGT"


def random_umi(length: int) -> str:
    """Return a random nucleotide string using os.urandom (CSPRNG)."""
    raw = os.urandom(length)
    return "".join(chr(NUC[b & 0x03]) for b in raw)


def output_path(inp: str) -> str:
    p = Path(inp)
    if p.suffixes[-2:] == [".fastq", ".gz"]:
        stem = p.with_suffix("").with_suffix("")  # strip .fastq.gz
        return str(stem) + "_umi.fastq.gz"
    elif p.suffix == ".gz":
        stem = p.with_suffix("")
        return str(stem.with_suffix("")) + "_umi" + "".join(p.suffixes)
    else:
        return str(p.with_suffix("")) + "_umi" + p.suffix


def process_random(inpath: str, umi_len: int, qual_char: str) -> None:
    """Prepend a randomly generated UMI to every read in *inpath*."""
    outpath = output_path(inpath)
    fixed_qual = qual_char * umi_len
    n_reads = 0

    with dnaio.open(inpath) as fin, dnaio.open(outpath, mode="w") as fout:
        for rec in fin:
            umi = random_umi(umi_len)
            out_rec = dnaio.SequenceRecord(
                rec.name,
                umi + rec.sequence,
                fixed_qual + rec.qualities,
            )
            fout.write(out_rec)
            n_reads += 1

    print(f"  {inpath}  →  {outpath}  ({n_reads:,} reads)")


def process_from_umi_fastq(inpath: str, umi_fastq: str, qual_char: str) -> None:
    """Prepend UMI sequences from *umi_fastq* to reads in *inpath*.

    dnaio's paired reader validates that read names match at every position
    (modulo /1 /2 suffixes and space-separated descriptions), so mismatched
    files are caught immediately rather than silently producing wrong output.
    The UMI base qualities from the UMI FASTQ are preserved; *qual_char* is
    only used as a fallback if a UMI record has no quality string.
    """
    outpath = output_path(inpath)
    n_reads = 0

    with dnaio.open(inpath, umi_fastq) as paired, dnaio.open(outpath, mode="w") as fout:
        for main_rec, umi_rec in paired:
            umi_qual = (
                umi_rec.qualities
                if umi_rec.qualities
                else qual_char * len(umi_rec.sequence)
            )
            out_rec = dnaio.SequenceRecord(
                main_rec.name,
                umi_rec.sequence + main_rec.sequence,
                umi_qual + main_rec.qualities,
            )
            fout.write(out_rec)
            n_reads += 1

    print(f"  {inpath}  →  {outpath}  ({n_reads:,} reads)")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Prepend a UMI to every read in FASTQ file(s).  "
            "By default a random UMI is generated per read.  "
            "Use --umi-fastq to take UMI sequences from a matching FASTQ file."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("fastq", nargs="+", help="Input FASTQ file(s) (.fastq or .fastq.gz)")
    ap.add_argument(
        "--umi-fastq",
        metavar="FILE",
        help=(
            "FASTQ file whose read sequences are used as UMIs (and whose "
            "qualities are prepended to the main read qualities).  "
            "Read names must match the main input file at every position; "
            "exactly one input FASTQ must be provided when this option is used."
        ),
    )
    ap.add_argument(
        "--umi-length",
        type=int,
        default=12,
        help="Length of randomly generated UMIs (ignored with --umi-fastq).",
    )
    ap.add_argument(
        "--qual-char",
        default="F",
        help=(
            "Quality character for randomly generated UMI bases "
            "(Phred+33; 'F' = Phred 37).  "
            "Also used as fallback when a UMI FASTQ record has no quality string."
        ),
    )
    args = ap.parse_args()

    if args.umi_fastq and len(args.fastq) != 1:
        ap.error("--umi-fastq requires exactly one input FASTQ file.")

    if args.umi_fastq:
        print(f"Mode       : UMI from file")
        print(f"UMI source : {args.umi_fastq}")
        print()
        process_from_umi_fastq(args.fastq[0], args.umi_fastq, args.qual_char)
    else:
        print(f"Mode       : random UMI")
        print(f"UMI length : {args.umi_length}")
        print(f"UMI quality: {args.qual_char * args.umi_length}  ('{args.qual_char}' × {args.umi_length})")
        print()
        for fq in args.fastq:
            process_random(fq, args.umi_length, args.qual_char)

    print("\nDone.")


if __name__ == "__main__":
    main()

