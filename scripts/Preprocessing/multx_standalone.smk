# snakemake --snakefile multx.smk --use-conda --config samples=path/to/samples.txt outdir=output/dir envfile=env.yaml container="oras://your/container:tag" paired=paired
import os
import re

# Config parameters
SAMPLES = config.get("samples", "samples.txt")
OUTDIR = config.get("outdir", "output")
ENVFILE = config.get("envfile", None)
CONTAINER = config.get("container", None)
PAIRED = config.get("paired", "paired")  # "paired" or "single"
THREADS = config.get("threads", 4)
# Read sample table (expects columns: sample, r1, r2)
sample_table = [line.strip().split() for line in open(SAMPLES) if line.strip() and not line.startswith("#")]
SAMPLE_LIST = [row[0] for row in sample_table]
SAMPLE_R1 = [row[1] for row in sample_table]
SAMPLE_R2 = [row[2] for row in sample_table]

# Whitelist path (assume one for all samples, can be customized)
WHITELIST = config.get("whitelist", "Multx_whitelist.txt")
demux_table = [line.strip().split() for line in open(WHITELIST) if line.strip() and not line.startswith("#")]
DEMUX_LIST = [row[0] for row in demux_table]

wildcard_constraints:
    outdir = OUTDIR,
    sample = '|'.join([re.escape(x) for x in DEMUX_LIST]),
    sample_mux = '|'.join([re.escape(x) for x in SAMPLE_LIST])

if PAIRED == "paired":
    rule all:
        input:
            i1 = expand("{outdir}/final/{sample}_R1.fastq.gz", outdir=OUTDIR, sample=DEMUX_LIST),
            i2 = expand("{outdir}/final/{sample}_R2.fastq.gz", outdir=OUTDIR, sample=DEMUX_LIST)

    rule multx_demux_first:
        input:
            r1 = lambda wildcards: SAMPLE_R1,
            r2 = lambda wildcards: SAMPLE_R2,
            whitelist = WHITELIST
        output:
            o1 = temp(expand("{outdir}/first/{sample}_R1_demux.fastq.gz", outdir=OUTDIR, sample=DEMUX_LIST)),
            o2 = temp(expand("{outdir}/first/{sample}_R2_demux.fastq.gz", outdir=OUTDIR, sample=DEMUX_LIST)),
            unmatched_r1 = temp(expand("{outdir}/first/{sample_mux}_unmatched_R1.fastq.gz", outdir=OUTDIR, sample_mux=SAMPLE_LIST)),
            unmatched_r2 = temp(expand("{outdir}/first/{sample_mux}_unmatched_R2.fastq.gz", outdir=OUTDIR, sample_mux=SAMPLE_LIST))       
        threads: THREADS
        conda: ENVFILE if ENVFILE else None
        container: CONTAINER if CONTAINER else None
        params:
            outdir = OUTDIR,            
        shell:
            """
            fastq-multx -B {input.whitelist} -b <(zcat {input.r1}) <(zcat {input.r2}) -o {params.outdir}/first/%_R1_demux.fastq.gz -o {params.outdir}/first/%_R2_demux.fastq.gz &> LOGS/multx_first.log &&
            mv -f {params.outdir}/first/unmatched_R1_demux.fastq.gz {output.unmatched_r1} &&
            mv -f {params.outdir}/first/unmatched_R2_demux.fastq.gz {output.unmatched_r2} &&
            touch {output.o1} {output.o2} {output.unmatched_r1} {output.unmatched_r2}
            """

    rule fastp_trim_unmatched:
        input:
            r1 = "{outdir}/first/{sample_mux}_unmatched_R1.fastq.gz",
            r2 = "{outdir}/first/{sample_mux}_unmatched_R2.fastq.gz"
        output:
            r1_trimmed = temp("{outdir}/first/{sample_mux}_unmatched_R1_trimmed.fastq.gz"),
            r2_trimmed = temp("{outdir}/first/{sample_mux}_unmatched_R2_trimmed.fastq.gz")        
        threads: THREADS
        conda: ENVFILE if ENVFILE else None
        container: CONTAINER if CONTAINER else None
        shell:
            """
            fastp --in1 {input.r1} --in2 {input.r2} --out1 {output.r1_trimmed} --out2 {output.r2_trimmed} --trim_front1 1 --disable_adapter_trimming --disable_quality_filtering --disable_trim_poly_g --disable_length_filtering --dont_eval_duplication --thread {threads} &> LOGS/fastp.log
            """

    rule multx_demux_second:
        input:
            r1 = expand("{outdir}/first/{sample_mux}_unmatched_R1_trimmed.fastq.gz", outdir=OUTDIR, sample_mux=SAMPLE_LIST),
            r2 = expand("{outdir}/first/{sample_mux}_unmatched_R2_trimmed.fastq.gz", outdir=OUTDIR, sample_mux=SAMPLE_LIST),
            whitelist = WHITELIST
        output:
            o1 = temp("{outdir}/second/{sample}_R1_demux.fastq.gz"),
            o2 = temp("{outdir}/second/{sample}_R2_demux.fastq.gz")        
        threads: THREADS
        conda: ENVFILE if ENVFILE else None
        container: CONTAINER if CONTAINER else None
        params:
            outdir = OUTDIR
        shell:
            """
            fastq-multx -B {input.whitelist} -b <(zcat {input.r1}) <(zcat {input.r2}) -o {params.outdir}/second/%_R1_demux_second.fastq.gz -o {params.outdir}/second/%_R2_demux_second.fastq.gz &> LOGS/multx_second.log &&
            touch {output.o1} {output.o2}
            """

    rule concat_final:
        input:
            i1_first = "{outdir}/first/{sample}_R1_demux.fastq.gz",
            i1_second = "{outdir}/first/{sample}_R2_demux.fastq.gz",
            i2_first = "{outdir}/second/{sample}_R1_demux.fastq.gz",
            i2_second = "{outdir}/second/{sample}_R2_demux.fastq.gz"
        output:
            o1 = "{outdir}/final/{sample}_R1.fastq.gz",
            o2 = "{outdir}/final/{sample}_R2.fastq.gz"
        shell:
            """
            cat {input.i1_first} {input.i1_second} > {output.o1}
            cat {input.i2_first} {input.i2_second} > {output.o2}
            """
else:
    rule all:
        input:
            i1 = expand("{outdir}/final/{sample}.fastq.gz", outdir=OUTDIR, sample=DEMUX_LIST)           

    rule multx_demux_first:
        input:
            r1 = lambda wildcards: SAMPLE_R1,            
            whitelist = WHITELIST
        output:
            o1 = temp(expand("{outdir}/first/{sample}_demux.fastq.gz", outdir=OUTDIR, sample=DEMUX_LIST)),            
            unmatched_r1 = temp(expand("{outdir}/first/{sample_mux}_unmatched.fastq.gz", outdir=OUTDIR, sample_mux=SAMPLE_LIST))        
        threads: THREADS
        conda: ENVFILE if ENVFILE else None
        container: CONTAINER if CONTAINER else None
        params:
            outdir = OUTDIR         
        shell:
            """
            fastq-multx -B {input.whitelist} -b <(zcat {input.r1}) -o {params.outdir}/first/%_demux.fastq.gz  &> LOGS/multx_first.log &&
            mv -f {params.outdir}/first/unmatched_demux.fastq.gz {output.unmatched_r1} &&            
            touch {output.o1} {output.unmatched_r1}
            """

    rule fastp_trim_unmatched:
        input:
            r1 = "{outdir}/first/{sample_mux}_unmatched.fastq.gz"          
        output:
            r1_trimmed = temp("{outdir}/first/{sample_mux}_unmatched_trimmed.fastq.gz")           
        threads: THREADS
        conda: ENVFILE if ENVFILE else None
        container: CONTAINER if CONTAINER else None
        shell:
            """
            fastp --in1 {input.r1} --out1 {output.r1_trimmed} --trim_front1 1 --disable_adapter_trimming --disable_quality_filtering --disable_trim_poly_g --disable_length_filtering --dont_eval_duplication --thread {threads} &> LOGS/fastp.log
            """

    rule multx_demux_second:
        input:
            r1 = expand("{outdir}/first/{sample_mux}_unmatched_trimmed.fastq.gz", outdir=OUTDIR, sample_mux=SAMPLE_LIST),            
            whitelist = WHITELIST
        output:
            o1 = temp("{outdir}/second/{sample}_demux.fastq.gz")            
        threads: THREADS
        conda: ENVFILE if ENVFILE else None
        container: CONTAINER if CONTAINER else None
        params:
            outdir = OUTDIR
        shell:
            """
            fastq-multx -B {input.whitelist} -b <(zcat {input.r1}) -o {params.outdir}/second/%_R1_demux_second.fastq.gz &> LOGS/multx_second.log &&
            touch {output.o1}
            """

    rule concat_final:
        input:
            i1_first = "{outdir}/first/{sample}_demux.fastq.gz",            
            i1_second = "{outdir}/second/{sample}_demux.fastq.gz"
        output:
            o1 = "{outdir}/final/{sample}.fastq.gz"
        shell:
            """
            cat {input.i1_first} {input.i1_second} > {output.o1}
            """