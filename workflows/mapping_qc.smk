include: "header.smk"

rule all:
    input: expand("DONE/{file}_mapped", file=samplecond(SAMPLES,config)),#file=SAMPLES),#samplecond(SAMPLES,config)),
           expand("QC/{rawfile}_qc.zip", rawfile=SAMPLES), "QC/Multi/multiqc_report.html" if "OFF" not in config["QC"] else []

if config['MAPPING'] is 'paired':
    paired = '_paired'
else:
    paired = ''

if 'OFF' not in config['QC']:
    include: str(config['QC'])+paired+'.smk'

if 'OFF' not in config['TRIMMENV']:
    include: str(config['TRIMMENV'])+paired+'.smk'

include: str(config['MAPPERENV'])+paired+'.smk'

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
            fn = lambda wildcards: "{fn}".format(fn=str(sample_from_path(wildcards.file)))
    shell: "zcat {input[0]} | samtools view -bS - | samtools sort -T {params.fn} -o {output[0]} --threads {threads} && samtools index {output[0]} 2> {log}"

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
            fn = lambda wildcards: "{fn}".format(fn=sample_from_path(wildcards.file))
    shell: "zcat {input[0]} | samtools view -bS - | samtools sort -T {params.fn} -o {output[0]} --threads {threads} && samtools index {output[0]} 2> {log}"

onsuccess:
    print("Workflow finished, no error")

rule themall:
    input:  rules.sam2bamuniq.output
    output: "DONE/{file}_mapped"
    run:
        for f in output:
            with open(f, "w") as out:
                out.write("DONE")
