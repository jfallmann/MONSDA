include: "header.smk"

rule all:
    input: expand("DONE/{file}_mapped", file=samplecond(SAMPLES,config)),#file=SAMPLES),#samplecond(SAMPLES,config)),
           expand("QC/{rawfile}_qc.zip", rawfile=SAMPLES) if "ON" in config["QC"] else [],
           "QC/Multi/DONE"

#include: "porechop.smk" #We lack a proper adapter trimming tool
if 'ON' in config['QC']:
    include: "fastqc.smk"           #Need alternative for really long reads

rule mapping:
    input:  expand("TRIMMED_FASTQ/{rawfile}_trimmed.fastq.gz",rawfile=SAMPLES)
    output: report("MAPPED/{file}_mapped.sam", category="MAPPING"),
            "UNMAPPED/{file}_unmapped.fastq.gz"
    log:    "LOGS/{file}/mapping.log"
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: 20
    params: mpara = lambda wildcards: ''.join(mapping_params(wildcards.file, None ,config)[0]),
            index = lambda wildcards: "{ref}/{gen}{name}.idx".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=''.join(mapping_params(wildcards.file, None ,config)[1])),
            ref = lambda wildcards: check_ref("{ref}/{gen}{name}.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=NAME['genomic'])),
            mapp=MAPPERBIN
    shell: "{params.mapp} -t {threads} {params.mpara} {params.index} {params.ref} {input[0]} | tee >(grep -v -P '\t4\t' > {output[0]}) >(grep -P '\t4\t' |samtools fastq -n - | pigz > {output[1]}) 1>/dev/null 2>> {log}"

rule gzipsam:
    input:  rules.mapping.output
    output: report("MAPPED/{file}_mapped.sam.gz", category="ZIPIT")
    log:    "LOGS/{file}/gzipsam.log"
    conda:  "../envs/base.yaml"
    threads: 20
    shell: "pigz -k -p {threads} -f {input} > {output} 2> {log}"

rule sortsam:
    input:  rules.gzipsam.output
    output: report("SORTED_MAPPED/{file}_mapped_sorted.sam.gz", category="SORTING"),
            temp("SORTED_MAPPED/{file}_mapped_header.gz"),
            temp("SORTTMP/{file}")
    log:    "LOGS/{file}/sortsam.log"
    conda: "../envs/samtools.yaml"
    threads: 20
    shell: "set +o pipefail;samtools view -H {input[0]}|grep -P '^@HD' |pigz -p {threads} -f > {output[1]} ; samtools view -H {input[0]}|grep -P '^@SQ'|sort -t$'\t' -k1,1 -k2,2V |pigz -p {threads} -f >> {output[1]} ; samtools view -H {input[0]}|grep -P '^@RG'|pigz -p {threads} -f >> {output[1]} ; samtools view -H {input[0]}|grep -P '^@PG'|pigz -p {threads} -f >> {output[1]} ; export LC_ALL=C;zcat {input[0]} | grep -v \"^@\"|sort --parallel={threads} -S 25% -T SORTTMP -t$'\t' -k3,3V -k4,4n - |pigz -p {threads} -f > {output[2]} ; cat {output[1]} {output[2]} > {output[0]} 2> {log}"

rule sam2bam:
    input:  rules.sortsam.output
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
    input: rules.sortsam.output,
           rules.sam2bam.output
    output: report("UNIQUE_MAPPED/{file}_mapped_sorted_unique.sam.gz", category="UNIQUE")
    log: "LOGS/{file}/uniqsam.log"
    conda: "../envs/base.yaml"
    threads: 20
    params: bins=BINS
    shell:  "{params.bins}/Shells/UniqueSam_woPicard.sh {input[0]} {output[0]} {threads} 2> {log}"

rule sam2bamuniq:
    input: rules.uniqsam.output,
           rules.sam2bam.output
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
#    input:  snakemake.utils.listfiles("QC/{file}*_gc.zip", restriction=None, omit_value=None)
    input:  expand("QC/{qcfile}_qc.zip", qcfile=SAMPLES),
            expand("QC/{qcfile}_trimmed_qc.zip", qcfile=SAMPLES),
            expand("QC/{file}_mapped_sorted_qc.zip", file=samplecond(SAMPLES,config)),
            expand("QC/{file}_mapped_sorted_unique_qc.zip", file=samplecond(SAMPLES,config))
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
    input:  rules.sam2bamuniq.output
    output: "DONE/{file}_mapped"
    run:
        for f in output:
            with open(f, "w") as out:
                out.write("DONE")
