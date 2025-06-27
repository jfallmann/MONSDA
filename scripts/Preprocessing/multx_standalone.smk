# snakemake --snakefile multx.smk --config samples=path/to/samples.txt outdir=output/dir envfile=env.yaml container="oras://your/container:tag" paired=paired
import os

# Config parameters
SAMPLES = config.get("samples", "samples.txt")
OUTDIR = config.get("outdir", "output")
ENVFILE = config.get("envfile", None)
CONTAINER = config.get("container", None)
PAIRED = config.get("paired", "paired")  # "paired" or "single"

# Read sample table (expects columns: sample, r1, r2)
sample_table = [line.strip().split() for line in open(SAMPLES) if line.strip() and not line.startswith("#")]
SAMPLE_LIST = [row[0] for row in sample_table]
SAMPLE_R1 = {row[0]: row[1] for row in sample_table}
SAMPLE_R2 = {row[0]: row[2] for row in sample_table} if PAIRED == "paired" else {}

# Whitelist path (assume one for all samples, can be customized)
WHITELIST = config.get("whitelist", "Multx_whitelist.txt")

# Helper for env/container
def env_directive():
    if ENVFILE:
        return f'conda: "{ENVFILE}"'
    elif CONTAINER:
        return f'container: "{CONTAINER}"'
    else:
        return ""

env_line = env_directive()

if PAIRED == "paired":
    rule all:
        input:
            expand("{outdir}/{sample}_R1.final.fastq.gz", outdir=OUTDIR, sample=SAMPLE_LIST),
            expand("{outdir}/{sample}_R2.final.fastq.gz", outdir=OUTDIR, sample=SAMPLE_LIST)

    rule multx_demux_first:
        input:
            r1 = lambda wildcards: SAMPLE_R1[wildcards.sample],
            r2 = lambda wildcards: SAMPLE_R2[wildcards.sample],
            whitelist = WHITELIST
        output:
            o1 = temp("{outdir}/{sample}_R1.demux.fastq.gz"),
            o2 = temp("{outdir}/{sample}_R2.demux.fastq.gz"),
            unmatched_r1 = temp("{outdir}/{sample}_unmatched_R1.fastq.gz"),
            unmatched_r2 = temp("{outdir}/{sample}_unmatched_R2.fastq.gz")
        log:
            "{outdir}/{sample}_multx_first.log"
        threads: 4
        # Dynamically set conda or container
        __env__ = env_line
        shell:
            """
            multx-fastq --threads {threads} --barcodes {input.whitelist} --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.o1} --out2 {output.o2} --unmatched1 {output.unmatched_r1} --unmatched2 {output.unmatched_r2} \
            > {log} 2>&1
            """

    rule fastp_trim_unmatched:
        input:
            r1 = rules.multx_demux_first.output.unmatched_r1
        output:
            r1_trimmed = temp("{outdir}/{sample}_unmatched_R1_trimmed.fastq.gz")
        threads: 1
        __env__ = env_line
        shell:
            """
            fastp --in1 {input.r1} --out1 {output.r1_trimmed} --trim_front1 1 --disable_adapter_trimming \
            --disable_quality_filtering --disable_trim_poly_g --disable_length_filtering --dont_eval_duplication \
            --thread {threads}
            """

    rule multx_demux_second:
        input:
            r1 = rules.fastp_trim_unmatched.output.r1_trimmed,
            r2 = rules.multx_demux_first.output.unmatched_r2,
            whitelist = WHITELIST
        output:
            o1 = temp("{outdir}/{sample}_R1.demux_second.fastq.gz"),
            o2 = temp("{outdir}/{sample}_R2.demux_second.fastq.gz")
        log:
            "{outdir}/{sample}_multx_second.log"
        threads: 4
        __env__ = env_line
        shell:
            """
            multx-fastq --threads {threads} --barcodes {input.whitelist} --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.o1} --out2 {output.o2} > {log} 2>&1
            """

    rule concat_final:
        input:
            o1_first = rules.multx_demux_first.output.o1,
            o1_second = rules.multx_demux_second.output.o1,
            o2_first = rules.multx_demux_first.output.o2,
            o2_second = rules.multx_demux_second.output.o2
        output:
            o1 = "{outdir}/{sample}_R1.final.fastq.gz",
            o2 = "{outdir}/{sample}_R2.final.fastq.gz"
        shell:
            """
            cat {input.o1_first} {input.o1_second} > {output.o1}
            cat {input.o2_first} {input.o2_second} > {output.o2}
            """

else:
    rule all:
        input:
            expand("{outdir}/{sample}.final.fastq.gz", outdir=OUTDIR, sample=SAMPLE_LIST)

    rule multx_demux_first:
        input:
            r1 = lambda wildcards: SAMPLE_R1[wildcards.sample],
            whitelist = WHITELIST
        output:
            o1 = temp("{outdir}/{sample}.demux.fastq.gz"),
            unmatched_r1 = temp("{outdir}/{sample}_unmatched.fastq.gz")
        log:
            "{outdir}/{sample}_multx_first.log"
        threads: 4
        __env__ = env_line
        shell:
            """
            multx-fastq --threads {threads} --barcodes {input.whitelist} --in1 {input.r1} \
            --out1 {output.o1} --unmatched1 {output.unmatched_r1} > {log} 2>&1
            """

    rule fastp_trim_unmatched:
        input:
            r1 = rules.multx_demux_first.output.unmatched_r1
        output:
            r1_trimmed = temp("{outdir}/{sample}_unmatched_trimmed.fastq.gz")
        threads: 1
        __env__ = env_line
        shell:
            """
            fastp --in1 {input.r1} --out1 {output.r1_trimmed} --trim_front1 1 --disable_adapter_trimming \
            --disable_quality_filtering --disable_trim_poly_g --disable_length_filtering --dont_eval_duplication \
            --thread {threads}
            """

    rule multx_demux_second:
        input:
            r1 = rules.fastp_trim_unmatched.output.r1_trimmed,
            whitelist = WHITELIST
        output:
            o1 = temp("{outdir}/{sample}.demux_second.fastq.gz")
        log:
            "{outdir}/{sample}_multx_second.log"
        threads: 4
        __env__ = env_line
        shell:
            """
            multx-fastq --threads {threads} --barcodes {input.whitelist} --in1 {input.r1} \
            --out1 {output.o1} > {log} 2>&1
            """

    rule concat_final:
        input:
            o1_first = rules.multx_demux_first.output.o1,
            o1_second = rules.multx_demux_second.output.o1
        output:
            o1 = "{outdir}/{sample}.final.fastq.gz"
        shell:
            """
            cat {input.o1_first} {input.o1_second} > {output.o1}
            """