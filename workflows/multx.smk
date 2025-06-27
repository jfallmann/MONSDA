include: "fastp.smk"

MULTXBIN, MULTXENV = env_bin_from_config(config, 'DEMULTIPLEXING')

if paired == 'paired':
    rule multx_demux_first:
        input:
            r1 = lambda wildcards: "PREFASTQ/{rawfile}_R1.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
            r2 = lambda wildcards: "PREFASTQ/{rawfile}_R2.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
            whitelist = "Multx_whitelist.txt"
        output:
            o1 = temp("FASTQ/{combo}/{file}_R1_demux.fastq.gz"),
            o2 = temp("FASTQ/{combo}/{file}_R2_demux.fastq.gz"),
            unmatched_r1 = temp("FASTQ/{combo}/{file}_unmatched_R1.fastq.gz"),
            unmatched_r2 = temp("FASTQ/{combo}/{file}_unmatched_R2.fastq.gz")
        log:
            "LOGS/{combo}/{file}_multx_first.log"
        conda: ""+MULTXENV+".yaml"
        container: "oras://jfallmann/monsda:"+MULTXENV+""
        threads: MAXTHREAD
        params:
            odir = lambda wildcards, output: os.path.dirname(output.o1),
            tpara = lambda wildcards: tool_params(wildcards.file, None, config, "DEMULTIPLEXING", MULTXENV).get('DEMUX', ""),
            multx = MULTXBIN
        shell:
            "{params.multx} --threads {threads} --barcodes {input.whitelist} --in1 {input.r1} --in2 {input.r2} "
            "--out1 {output.o1} --out2 {output.o2} --unmatched1 {output.unmatched_r1} --unmatched2 {output.unmatched_r2} {params.tpara} > {log} 2>&1"

    rule fastp_trim_unmatched:
        input:
            r1 = rules.multx_demux_first.output.unmatched_r1
        output:
            r1_trimmed = "DEMUX_FASTQ/{combo}/{file}_unmatched_R1_trimmed.fastq.gz"
        conda: ""+MULTXENV+".yaml"
        container: "oras://jfallmann/monsda:"+MULTXENV+""
        threads: 1
        shell:
            "fastp --in1 {input.r1} --out1 {output.r1_trimmed} --trim_front1 1 --disable_adapter_trimming --disable_quality_filtering --disable_trim_poly_g --disable_length_filtering --dont_eval_duplication --thread {threads}"

    rule multx_demux_second:
        input:
            r1 = rules.fastp_trim_unmatched.output.r1_trimmed,
            r2 = rules.multx_demux_first.output.unmatched_r2,
            whitelist = "Multx_whitelist.txt"
        output:
            o1 = temp("FASTQ/{combo}/{file}_R1_demux_second.fastq.gz"),
            o2 = temp("FASTQ/{combo}/{file}_R2_demux_second.fastq.gz")
        log:
            "LOGS/{combo}/{file}_multx_second.log"
        conda: ""+MULTXENV+".yaml"
        container: "oras://jfallmann/monsda:"+MULTXENV+""
        threads: MAXTHREAD
        params:
            odir = lambda wildcards, output: os.path.dirname(output.o1),
            tpara = lambda wildcards: tool_params(wildcards.file, None, config, "DEMULTIPLEXING", MULTXENV).get('DEMUX', ""),
            multx = MULTXBIN
        shell:
            "{params.multx} --threads {threads} --barcodes {input.whitelist} --in1 {input.r1} --in2 {input.r2} --out1 {output.o1} --out2 {output.o2} {params.tpara} > {log} 2>&1"

    rule concat_final:
        input:
            o1_first = rules.multx_demux_first.output.o1,
            o1_second = rules.multx_demux_second.output.o1,
            o2_first = rules.multx_demux_first.output.o2,
            o2_second = rules.multx_demux_second.output.o2
        output:
            o1 = "FASTQ/{combo}/{file}_R1.fastq.gz",
            o2 = "FASTQ/{combo}/{file}_R2.fastq.gz"
        shell:
            "cat {input.o1_first} {input.o1_second} > {output.o1} && cat {input.o2_first} {input.o2_second} > {output.o2}"

# You can add similar logic for single-end if needed.