include: "header.smk"

#include: "porechop.smk" #We lack a proper adapter trimming tool
include: "fastqc.smk"           #Need alternative for really long reads

rule all:
    input: expand("DONE/{file}_mapped",file=SAMPLES),#samplecond(SAMPLES,config)),
           "QC/Multi/DONE"
        #input: rules.multiqc.output #not working, rule multiqc not yet defined

rule mapping:
    input:  rules.qc_trimmed.input, rules.qc_trimmed.output if "ON" in config["QC"] else "TRIMMED_FASTQ/{file}_trimmed.fastq.gz"
    output: report("MAPPED/{cond}/{file}_mapped.sam", category="MAPPING"),
            "UNMAPPED/{cond}/{file}_unmapped.fastq"
    log:    "LOGS/{cond}/{file}/mapping.log"
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: 20
    params: mpara = lambda wildcards: mapping_params(wildcards.file, None ,config),
            index = lambda wildcards: "{ref}/{gen}{name}.idx".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=NAME),
            ref = lambda wildcards: check_ref("{ref}/{gen}{name}.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=NAME['genomic'])),
            mapp=MAPPERBIN,
            of = lambda wildcards: "{out}".format(out="MAPPED/"+samplecond(wildcards.file,config))
    shell: "arr=({params.mpara});alen=$({{#arr[@]}});for i in \"${{!arr[@]}}\";do {params.mapp} {params.mpara[$i][1]} {params.ref} {input[0]} | tee -a >(grep -v -P '\t4\t' >{params.of}_mapped.sam) | grep -P '\t4\t' > UN{params.of}_unmapped.fastq) 2>> {log};done"
#           arr=({input}); alen=${{#arr[@]}}; orr=({output}); for i in \"${{!arr[@]}}\";do a=$(zcat ${{arr[$i]}}|wc -l ); echo $((a/4)) > ${{orr[$i]}};done

rule gzipsam:
    input: rules.mapping.output
    output: report("MAPPED/{cond}/{file}_mapped.sam.gz", category="ZIPIT"),
            "UNMAPPED/{cond}/{file}_unmapped.fastq.gz"
    log:    "LOGS/{cond}/{file}/gzipsam.log"
    conda:  "../envs/base.yaml"
    threads: 20
    shell: "pigz -k -p {threads} -f {input} > {output} 2> {log}"

rule sortsam:
    input:  rules.gzipsam.output
    output: report("SORTED_MAPPED/{cond}/{file}_mapped_sorted.sam.gz", category="SORTING"),
            temp("SORTED_MAPPED/{cond}/{file}_mapped_header.gz"),
            temp("SORTTMP/{cond}/{file}")
    log:    "LOGS/{cond}/{file}/sortsam.log"
    conda: "../envs/samtools.yaml"
    threads: 20
    shell: "set +o pipefail;samtools view -H {input[0]}|grep -P '^@HD' |pigz -p {threads} -f > {output[1]} ; samtools view -H {input[0]}|grep -P '^@SQ'|sort -t$'\t' -k1,1 -k2,2V |pigz -p {threads} -f >> {output[1]} ; samtools view -H {input[0]}|grep -P '^@RG'|pigz -p {threads} -f >> {output[1]} ; samtools view -H {input[0]}|grep -P '^@PG'|pigz -p {threads} -f >> {output[1]} ; export LC_ALL=C;zcat {input[0]} | grep -v \"^@\"|sort --parallel={threads} -S 25% -T SORTTMP -t$'\t' -k3,3V -k4,4n - |pigz -p {threads} -f > {output[2]} ; cat {output[1]} {output[2]} > {output[0]} 2> {log}"

rule sam2bam:
    input:  rules.sortsam.output
    output: report("SORTED_MAPPED/{cond}/{file}_mapped_sorted.bam", category="2BAM"),
            "SORTED_MAPPED/{cond}/{file}_mapped_sorted.bam.bai"
    log:    "LOGS/{cond}/{file}/sam2bam.log"
    conda: "../envs/samtools.yaml"
    threads: 20
    params: bins = BINS,
#            sizes = lambda wildcards: "{ref}/{gen}.{name}.chrom.sizes".format(ref=REFERENCE, gen=genomepath(wildcards.file,config), name=NAME),
            fn = lambda wildcards: "{fn}".format(fn=str(sample_from_path(wildcards.file)))
    shell: "zcat {input[0]} | samtools view -bS - | samtools sort -T {params.fn} -o {output[0]} --threads {threads} && samtools index {output[0]} 2> {log}"
        #"{params.bins}/Shells/Sam2Bam.sh {input[0]} {params.sizes} {output}"

rule uniqsam:
    input: rules.sortsam.output,
           rules.sam2bam.output
    output: report("UNIQUE_MAPPED/{cond}/{file}_mapped_sorted_unique.sam.gz", category="UNIQUE")
    log: "LOGS/{cond}/{file}/uniqsam.log"
    conda: "../envs/base.yaml"
    threads: 20
    params: bins=BINS
    shell:  "{params.bins}/Shells/UniqueSam_woPicard.sh {input[0]} {output[0]} {threads} 2> {log}"

rule sam2bamuniq:
    input: rules.uniqsam.output,
           rules.sam2bam.output
    output:  report("UNIQUE_MAPPED/{cond}/{file}_mapped_sorted_unique.bam", category="2BAM"),
             "UNIQUE_MAPPED/{cond}/{file}_mapped_sorted_unique.bam.bai"
    log:     "LOGS/{cond}/{file}/sam2bamuniq.log"
    conda:   "../envs/samtools.yaml"
    threads: 20
    params: bins=BINS,
            #           sizes = lambda wildcards: "{ref}/{gen}.{name}.chrom.sizes".format(ref=REFERENCE, gen=genomepath(wildcards.file,config), name=NAME),
            fn = lambda wildcards: "{fn}".format(fn=sample_from_path(wildcards.file))
    shell: "zcat {input[0]} | samtools view -bS - | samtools sort -T {params.fn} -o {output[0]} --threads {threads} && samtools index {output[0]} 2> {log}"
           #"{params.bins}/Shells/Sam2Bam.sh {input[0]} {params.sizes} {output}"

rule multiqc:
    input:  snakemake.utils.listfiles("QC/{file}*_gc.zip", restriction=None, omit_value=None)
#    input:  expand("QC/{file}_qc.zip", file=SAMPLES),
#            expand("QC/{file}_trimmed_qc.zip", file=SAMPLES),
#            expand("QC/{file}_mapped_sorted_qc.zip", file=SAMPLES),
#            expand("QC/{file}_mapped_sorted_unique_qc.zip", file=SAMPLES)
#    input: rules.qc_raw.output,
#           rules.qc_trimmed.output,
#           rules.qc_mapped.output,
#           rules.qc_uniquemapped.output
    output: report("QC/Multi/DONE", category="QC")
    log:    "LOGS/multiqc.log"
    conda:  "../envs/qc.yaml"
    shell:  "OUT=$(dirname {output}); multiqc -k json -z -o $OUT $PWD 2> {log} && touch QC/Multi/DONE"

onsuccess:
    print("Workflow finished, no error")

rule themall:
    input:  rules.multiqc.output, rules.sam2bamuniq.output if "ON" in config["QC"] else rules.sam2bamuniq.output
    output: "DONE/{cond}/{file}_mapped"
    run:
        for f in output:
            with open(f, "w") as out:
                out.write("DONE")
