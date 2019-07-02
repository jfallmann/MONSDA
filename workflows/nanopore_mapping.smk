include: "header.smk"

#include: "porechop.smk" #We lack a proper adapter trimming tool
include: "fastqc.smk"           #Need alternative for really long reads

rule all:
    input: expand("DONE/{file}_mapped",file=SAMPLES),
           "QC/Multi/DONE"

rule mapping:
    input:  rules.fastqc_trimmed.input, rules.fastqc_trimmed.output if "ON" in config["QC"] else expand("TRIMMED_FASTQ/{file}_trimmed.fastq.gz",file=SAMPLES)
    output: report("MAPPED/{file}{runstate}_mapped.sam", category="MAPPING"),
            "UNMAPPED/{file}{runstate}_unmapped.fastq"
    log:    "LOGS/{file}{runstate}/mapping.log"
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: 20
    params: mpara = lambda wildcards: mapping_params(wildcards.file, None ,config),
            index = lambda wildcards: "{ref}/{gen}{name}.idx".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=NAME),
            ref = lambda wildcards: check_ref("{ref}/{gen}{name}.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=NAME['genomic'])),
            mapp=MAPPERBIN,
            of = lambda wildcards: "{out}".format(out="MAPPED/"+str(wildcards.file))
    shell: "arr=({params.mpara});alen=$({{#arr[@]}});for i in \"${{!arr[@]}}\";do {params.mapp} {params.mpara[$i][1]} {params.ref} {input[0]} | tee -a >(grep -v -P '\t4\t' >{params.of}_{params.mpara[$i][1]}_mapped.sam) | grep -P '\t4\t' > UN{params.of}_{params.mpara[$i][1]}_unmapped.fastq) 2>> {log};done"
#           arr=({input}); alen=${{#arr[@]}}; orr=({output}); for i in \"${{!arr[@]}}\";do a=$(zcat ${{arr[$i]}}|wc -l ); echo $((a/4)) > ${{orr[$i]}};done

rule gzipsam:
    input: "MAPPED/{file}_mapped.sam",
           "UNMAPPED/{file}_unmapped.fastq",
    output: report("MAPPED/{file}_mapped.sam.gz", category="ZIPIT"),
            "UNMAPPED/{file}_unmapped.fastq.gz"
    log:    "LOGS/{file}/gzipsam.log"
    conda:  "../envs/base.yaml"
    threads: 20
    shell: "pigz -k -p {threads} -f {input} > {output} 2> {log}"

rule sortsam:
    input:  "MAPPED/{file}_mapped.sam.gz"
    output: report("SORTED_MAPPED/{file}_mapped_sorted.sam.gz", category="SORTING"),
            temp("SORTED_MAPPED/{file}_mapped_header.gz"),
            temp("SORTTMP/{file}")
    log:    "LOGS/{file}/sortsam.log"
    conda: "../envs/samtools.yaml"
    threads: 20
    shell: "set +o pipefail;samtools view -H {input[0]}|grep -P '^@HD' |pigz -p {threads} -f > {output[1]} ; samtools view -H {input[0]}|grep -P '^@SQ'|sort -t$'\t' -k1,1 -k2,2V |pigz -p {threads} -f >> {output[1]} ; samtools view -H {input[0]}|grep -P '^@RG'|pigz -p {threads} -f >> {output[1]} ; samtools view -H {input[0]}|grep -P '^@PG'|pigz -p {threads} -f >> {output[1]} ; export LC_ALL=C;zcat {input[0]} | grep -v \"^@\"|sort --parallel={threads} -S 25% -T SORTTMP -t$'\t' -k3,3V -k4,4n - |pigz -p {threads} -f > {output[2]} ; cat {output[1]} {output[2]} > {output[0]} 2> {log}"

rule sam2bam:
    input:  "SORTED_MAPPED/{file}_mapped_sorted.sam.gz", "QC/{file}_mapped_sorted_qc.zip" if  "ON" in config["QC"] else "SORTED_MAPPED/{file}_mapped_sorted.sam.gz"
    output: report("SORTED_MAPPED/{file}_mapped_sorted.bam", category="2BAM"),
            "SORTED_MAPPED/{file}_mapped_sorted.bam.bai"
    log:    "LOGS/{file}/sam2bam.log"
    conda: "../envs/samtools.yaml"
    threads: 20
    params: bins = BINS,
#            sizes = lambda wildcards: "{ref}/{gen}.{name}.chrom.sizes".format(ref=REFERENCE, gen=genomepath(wildcards.file,config), name=NAME),
            fn = lambda wildcards: "{fn}".format(fn=str(sample_from_path(wildcards.file)))
    shell: "zcat {input[0]} | samtools view -bS - | samtools sort -T {params.fn} -o {output[0]} --threads {threads} && samtools index {output[0]} 2> {log}"
        #"{params.bins}/Shells/Sam2Bam.sh {input[0]} {params.sizes} {output}"

rule uniqsam:
    input: "SORTED_MAPPED/{file}_mapped_sorted.sam.gz",
           "SORTED_MAPPED/{file}_mapped_sorted.bam",
           "SORTED_MAPPED/{file}_mapped_sorted.bam.bai"
    output: report("UNIQUE_MAPPED/{file}_mapped_sorted_unique.sam.gz", category="UNIQUE")
    log: "LOGS/{file}/uniqsam.log"
    conda: "../envs/base.yaml"
    threads: 20
    params: bins=BINS
    shell:  "{params.bins}/Shells/UniqueSam_woPicard.sh {input[0]} {output[0]} {threads} 2> {log}"

rule sam2bamuniq:
    input:   "UNIQUE_MAPPED/{file}_mapped_sorted_unique.sam.gz",
             "SORTED_MAPPED/{file}_mapped_sorted.bam",
             "SORTED_MAPPED/{file}_mapped_sorted.bam.bai"
    output:  report("UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam", category="2BAM"),
             "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam.bai"
    log:     "LOGS/{file}/sam2bamuniq.log"
    conda:   "../envs/samtools.yaml"
    threads: 20
    params: bins=BINS,
            #           sizes = lambda wildcards: "{ref}/{gen}.{name}.chrom.sizes".format(ref=REFERENCE, gen=genomepath(wildcards.file,config), name=NAME),
            fn = lambda wildcards: "{fn}".format(fn=sample_from_path(wildcards.file))
    shell: "zcat {input[0]} | samtools view -bS - | samtools sort -T {params.fn} -o {output[0]} --threads {threads} && samtools index {output[0]} 2> {log}"
           #"{params.bins}/Shells/Sam2Bam.sh {input[0]} {params.sizes} {output}"

rule multiqc:
    input:  expand("QC/{file}_qc.zip", file=SAMPLES),
            expand("QC/{file}_trimmed_qc.zip", file=SAMPLES),
            expand("QC/{file}_mapped_sorted_qc.zip", file=SAMPLES),
            expand("QC/{file}_mapped_sorted_unique_qc.zip", file=SAMPLES)
    output: report("QC/Multi/DONE", category="QC")
    log:    "LOGS/multiqc.log"
    conda:  "../envs/qc.yaml"
    shell:  "OUT=$(dirname {output}); multiqc -k json -z -o $OUT $PWD 2> {log} && touch QC/Multi/DONE"

rule themall:
    input:  "QC/Multi/DONE", "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam" if "ON" in config["QC"] else "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
    output: "DONE/{file}_mapped"
    run:
        for f in output:
            with open(f, "w") as out:
                out.write("DONE")
