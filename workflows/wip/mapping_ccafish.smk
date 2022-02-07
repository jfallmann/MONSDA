import glob, os, sys, inspect, snakemake #, datetime

### snakemake --use-conda --ri --latency-wait 120 -j 16 -s mapping.smk --configfile config_dicty.json --directory ${PWD} -k
### optional with date
###snakemake --use-conda --ri --latency-wait 120 -j 16 -s mapping.smk --configfile config_dicty.json --directory ${PWD}/`date +%d-%m-%Y` -k

cmd_subfolder = os.path.join(os.path.dirname(os.path.realpath(os.path.abspath( inspect.getfile( inspect.currentframe() )) )),"../lib")
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from Collection import *

QC=config["QC"]
ADAPTERS=config["ADAPTERS"]
REFERENCE=config["REFERENCE"]
GENOME=config["GENOME"]
NAME=config["NAME"]
BINS=config["BINS"]
SOURCE=sources(config)
SAMPLES=samples(config)

include: "cutadapt.smk"

rule all:
    input:  expand("DONE/{file}_processed",file=SAMPLES),
            "QC/Multi/multiqc_report.html"

#rule sra2fastq:
#       input:  "RAW/{source}/{file}.sra"
#       output: "FASTQ/{source}/{file}.fastq.gz"
#       log:    "LOGS/sra2fastq_{file}.log"
#       shell:  "fastq-dump -Z {input} | gzip > {output}"

rule fastqc_raw:
    input:  "FASTQ/{file}.fastq.gz"
    output: report("QC/{file}_fastqc.zip", category="QC")
    conda: "../envs/qc.yaml"
    threads: 20
    params: dir= expand("QC/{source}",source=SOURCE)
    shell: "for i in {input}; do OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f fastq {input};done"

rule fastqc_trimmed:
    input:  "TRIMMED_FASTQ/{file}_trimmed.fastq.gz"
    output: report("QC/{file}_trimmed_fastqc.zip", category="QC")
    conda: "../envs/qc.yaml"
    threads: 20
    params: dir= expand("QC/{source}",source=SOURCE)
    shell: "for i in {input}; do OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f fastq {input};done"

rule segemehl_map:
    input:  "TRIMMED_FASTQ/{file}_trimmed.fastq.gz", "QC/{file}_trimmed_fastqc.zip" if "ON" in config["QC"]  else "TRIMMED_FASTQ/{file}_trimmed.fastq.gz"
    output: report("MAPPED/{file}_mapped.sam", category="MAPPING"),
            "UNMAPPED/{file}_unmapped.fastq"
    conda: "../envs/segemehl.yaml"
    threads: 20
    params: p=lambda wildcards:mapping_params(wildcards.file, "cluster", config),
            index = lambda wildcards: "{ref}/{gen}.{name}_cluster.idx".format(ref=REFERENCE, gen=genomepath(wildcards.file, config), name=NAME['tRNA']),
            ref = lambda wildcards: "{ref}/{gen}.{name}_cluster.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file, config), name=NAME['tRNA'])
    shell:  "segemehl.x {params.p} -d {params.ref} -i {params.index} -q {input[0]} --threads {threads} -o {output[0]} -u {output[1]}"

rule filter_tRNA:
    input:  "MAPPED/{file}_mapped.sam",
            "TRIMMED_FASTQ/{file}_trimmed.fastq.gz"
    output: "FILTERED/{file}_in.sam",
            "FILTERED/{file}_filtered.sam",
            report("FILTERED/{file}_filtered.fastq.gz", category="FILTER")
    conda: "../envs/perl.yaml"
    threads: 1
    params: bins = BINS,
            ref = lambda wildcards: "{ref}/{gen}.{name}_cluster.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file, config), name=NAME['tRNA'])
    shell:  "perl {params.bins}/Universal/Readfilter.pl {input[0]} {params.ref} {output[0]} {output[1]} && perl {params.bins}/Universal/sam2fastq.pl {output[1]} {input[1]} |gzip > {output[2]}"

rule segemehl_map_artificial:
    input:  "FILTERED/{file}_filtered.fastq.gz", "MAPPED/{file}_mapped.sam", "QC/{file}_trimmed_fastqc.zip" if "ON" in config["QC"]  else "FILTERED/{file}_filtered.fastq.gz"
    output: report("MAPPED/{file}_mapped_artificial.sam",category="MAPPING"),
            "UNMAPPED/{file}_unmapped_artificial.fastq"
    conda: "../envs/segemehl.yaml"
    threads: 20
    params: p=lambda wildcards: mapping_params(wildcards.file, "artificial", config),
            index = lambda wildcards: "{ref}/{gen}.{name}_artificial.idx".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=NAME['genomic']),,
            ref = lambda wildcards: "{ref}/{gen}.{name}_artificial.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=NAME['genomic'])
    shell:  "segemehl.x {params.p} -d {params.ref} -i {params.index} -q {input[0]} --threads {threads} -o {output[0]} -u {output[1]}"

rule remove_reads:
    input:   "MAPPED/{file}_mapped_artificial.sam",
             "FILTERED/{file}_filtered.sam"
    output:  "FILTERED/{file}_Gfiltered.sam",
             "FILTERED/{file}_GPfiltered.sam",
             report("FILTERED/{file}_almostFinal.sam", category="FILTER")
    conda: "../envs/perl.yaml"
    threads: 1
    params: bins = BINS,
            ref= lambda wildcards: "{ref}/{gen}.{name}_pre-tRNAs.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=NAME['tRNA']),
            bed = lambda wildcards: "{ref}/{gen}.{name}_pre-tRNAs.bed".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=NAME['tRNA'])
    shell:  "perl {params.bins}/Universal/removeGenomeMapper_2.pl {params.ref} {input[0]} {output[0]} && perl {params.bins}/Mapping/1a_removePrecursor2sam_Dicty.pl {params.bed} {output[0]} 30 > {output[1]} && perl {params.bins}/Universal/getFinalSam.pl {output[1]} {input[1]} > {output[2]}"

rule gzipsam:
    input:  "MAPPED/{file}_mapped.sam",
            "UNMAPPED/{file}_unmapped.fastq",
            "MAPPED/{file}_mapped_artificial.sam",
            "UNMAPPED/{file}_unmapped_artificial.fastq",
            "FILTERED/{file}_in.sam",
            "FILTERED/{file}_filtered.sam",
            "FILTERED/{file}_Gfiltered.sam",
            "FILTERED/{file}_GPfiltered.sam",
            "FILTERED/{file}_almostFinal.sam"
    output: report("MAPPED/{file}_mapped.sam.gz", category="ZIPIT"),
            "UNMAPPED/{file}_unmapped.fastq.gz",
            "MAPPED/{file}_mapped_artificial.sam.gz",
            "UNMAPPED/{file}_unmapped_artificial.fastq.gz",
            "FILTERED/{file}_in.sam.gz",
            "FILTERED/{file}_filtered.sam.gz",
            "FILTERED/{file}_Gfiltered.sam.gz",
            "FILTERED/{file}_GPfiltered.sam.gz",
            "FILTERED/{file}_almostFinal.sam.gz"
    conda:  "../envs/base.yaml"
    threads: 20
    shell:  "pigz -k -p {t} -f {f} > {o}".format(t=threads,f=input[i],o=output[i])

rule combine_reads:
    input: "FILTERED/{file}_in.sam.gz",
           "FILTERED/{file}_almostFinal.sam.gz"
    output: report("FINAL/{file}_final.sam.gz", category="COMBINE")
    conda: "../envs/base.yaml"
    threads: 1
    shell:  "cat {input[0]} {input[1]} > {output[0]}"

rule sortsam:
    input:  "FINAL/{file}_final.sam.gz",
    output: report("SORTED_MAPPED/{file}_mapped_sorted.sam.gz", category="SORTING"),
            temp("SORTED_MAPPED/{file}_mapped_header.gz"),
            temp("SORTTMP/{file}")
    conda: "../envs/samtools.yaml"
    threads: MAXTHREAD
    params: sortmem = MAXTHREAD*2.5
    shell:  "samtools view -H {input[0]}|grep '@HD' |pigz -p {threads} -f > {output[1]} && samtools view -H {input[0]}|grep '@SQ'|sort -t$'\t' -k1,1 -k2,2V |pigz -p {threads} -f >> {output[1]} && samtools view -H {input[0]}|grep '@RG'|pigz -p {threads} -f >> {output[1]} && samtools view -H {input[0]}|grep '@PG'|pigz -p {threads} -f >> {output[1]} && export LC_ALL=C;zcat {input[0]} | grep -v \"^@\"|sort --parallel={threads} -S {params.sortmem}% -T SORTTMP -t$'\t' -k3,3V -k4,4n - |pigz -p {threads} -f > {output[2]} && cat {output[1]} {output[2]} > {output[0]}"

### Fixing header order for Picardtools
# for i in *.sam.gz;do echo $i;zcat $i|samtools view -H|grep '@HD'|pigz > tmp.gz && zcat $i|samtools view -H|grep '@SQ'|sort -t$'\t' -k1,1 -k2,2V |pigz >> tmp.gz && zcat $i|samtools view -H|grep '@PG' |pigz >> tmp.gz && zcat $i|grep -v "^@"|pigz > temp.gz && cat tmp.gz temp.gz > $i\_rehead.gz;done

rule fastqc_mapped:
    input:  "SORTED_MAPPED/{file}_mapped_sorted.sam.gz"
    output: report("QC/{file}_mapped_sorted_fastqc.zip", category="QC")
    conda: "../envs/qc.yaml"
    threads: 20
    params: dir= expand("QC/{source}",source=SOURCE)
    shell: "for i in {input}; do OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f sam_mapped {input};done"

rule sam2bam:
    input:  "SORTED_MAPPED/{file}_mapped_sorted.sam.gz", "QC/{file}_mapped_sorted_fastqc.zip" if "ON" in config["QC"] else "SORTED_MAPPED/{file}_mapped_sorted.sam.gz"
    output: report("SORTED_MAPPED/{file}_mapped_sorted.bam", category="2BAM"),
            "SORTED_MAPPED/{file}_mapped_sorted.bam.bai"
    conda: "../envs/samtools.yaml"
    threads: 20
    params: bins = BINS
    #            ref = "{ref}/{gen}.{name}.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=NAME['tRNA'])
    run:  "zcat {input[0]} | samtools view -bS - | samtools sort -T {wildcards.file} -o {output[0]} --threads {threads} && samtools index {output[0]}"
#       "zcat {input[0]}|samtools view -bT {input[2]} -o {output[0]} --threads {threads} -"
#"java -jar -Xmx40g {params.bins}/picard.jar SortSam  VALIDATION_STRINGENCY=LENIENT  MAX_RECORDS_IN_RAM=7500000 TMP_DIR=SORTTMP  INPUT={input[0]}  OUTPUT={output}  SORTORDER=coordinate"

rule uniqsam:
    input:  "SORTED_MAPPED/{file}_mapped_sorted.sam.gz",
            "SORTED_MAPPED/{file}_mapped_sorted.bam",
            "SORTED_MAPPED/{file}_mapped_sorted.bam.bai"
    output: report("UNIQUE_MAPPED/{file}_mapped_sorted_unique.sam.gz", category="UNIQUE")
    conda: "../envs/base.yaml"
    threads: 20
    params: bins=BINS
    shell:  "{params.bins}/Shells/UniqueSam_woPicard.sh {input[0]} {output[0]} {threads}"
            #"samtools view -H <(zcat {input[0]})|grep '@HD' |pigz -p {threads} -f > {output[0]} && samtools view -H <(zcat {input[0]})|grep '@SQ'|sort -t$'\t' -k1,1 -k2,2V |pigz -p {threads} -f >> {output[0]} && samtools view -H <(zcat {input[0]})|grep '@RG'|pigz -p {threads} -f >> {output[0]} && samtools view -H <(zcat {input[0]})|grep '@PG'|pigz -p {threads} -f >> {output[0]} && {params.bins}/Shells/UniqueSam_woPicard.sh {input[0]} {params.bins} {output[0]}"

rule sam2bamuniq:
    input:  "UNIQUE_MAPPED/{file}_mapped_sorted_unique.sam.gz",
            "SORTED_MAPPED/{file}_mapped_sorted.bam",
            "SORTED_MAPPED/{file}_mapped_sorted.bam.bai"
    output: report("UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam", category="2BAM"),
            "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam.bai"
    conda: "../envs/samtools.yaml"
    threads: 20
    params: bins=BINS
#            ref = lambda wildcards: "{ref}/{gen}.{name}.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=NAME['tRNA'])
    shell:  "zcat {input[0]} | samtools view -bS - | samtools sort -T {wildcards.file} -o {output[0]} --threads {threads} && samtools index {output[0]}"
#"{params.bins}/Shells/Sam2Bam.sh {input[0]} {params.ref} {params.bins} {output[0]} {threads}"
#"zcat {input[0]}|samtools view -bT {input[1]} -o {output} --threads {threads} -"

rule phasesam:
    input:  "SORTED_MAPPED/{file}_mapped_sorted.sam.gz",
            "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
    output: report("PHASED_MAPPED/{file}_mapped_sorted_phased.sam.gz", category="PHASED"),
            temp("PHASED_MAPPED/{file}_mapped_sorted_phased.sam"),
            temp("SORTED_MAPPED/{file}_mapped_sorted.sam")
    conda: "../envs/perl.yaml"
    threads: 1
    params: bins=BINS
    shell: "zcat {input[0]} > {output[2]} && perl {params.bins}/Universal/multimapperPhasing.pl -ed 0 -id 0 -verbose 0 -sam {output[2]} -out {output[1]} && cat {output[1]}|pigz -p {threads} -f  > {output[0]}"

rule sam2bamphased:
    input:  "PHASED_MAPPED/{file}_mapped_sorted_phased.sam.gz",
            "SORTED_MAPPED/{file}_mapped_sorted.bam",
            "SORTED_MAPPED/{file}_mapped_sorted.bam.bai"
    output: report("PHASED_MAPPED/{file}_mapped_sorted_phased.bam", category="2BAM"),
            "PHASED_MAPPED/{file}_mapped_sorted_phased.bam.bai"
    conda: "../envs/samtools.yaml"
    threads: 20
    params: bins=BINS
#            ref = lambda wildcards: "{ref}/{gen}.{name}.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=NAME['tRNA'])
    shell:  "zcat {input[0]} | samtools view -bS - | samtools sort -T {wildcards.file} -o {output[0]} --threads {threads} && samtools index {output[0]}"

rule fastqc_uniquemapped:
    input:  "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam",
            "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam.bai"
    output: report("QC/{file}_mapped_sorted_unique_fastqc.zip", category="QC")
    conda: "../envs/qc.yaml"
    threads: 20
    params: dir=expand("QC/{source}",source=SOURCE)
    shell: "for i in {input[0]}; do OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f bam {input[0]};done"

rule fastqc_phasedmapped:
    input:  "PHASED_MAPPED/{file}_mapped_sorted_phased.bam",
            "PHASED_MAPPED/{file}_mapped_sorted_phased.bam.bai"
    output: report("QC/{file}_mapped_sorted_phased_fastqc.zip", category="QC")
    conda: "../envs/qc.yaml"
    threads: 20
    params: dir=expand("QC/{source}",source=SOURCE)
    shell: "for i in {input[0]}; do OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f bam {input[0]};done"

rule remove_dups:
    input:  "SORTED_MAPPED/{file}_mapped_sorted.bam", "QC/{file}_mapped_sorted_fastqc.zip" if "ON" in config["QC"] else "SORTED_MAPPED/{file}_mapped_sorted.bam"
    output: report("SORTED_MAPPED/{file}_mapped_sorted_dedup.bam", category="DEDUP"),
            "LOGS/{file}_dedup.txt",
            "SORTED_MAPPED/{file}_mapped_sorted_dedup.bam.bai"
    conda: "../envs/picardtools.yaml"
    threads: 1
    shell:  "picard MarkDuplicates I={input[0]} O={output[0]} M={output[1]} ASO=coordinate REMOVE_SEQUENCING_DUPLICATES=true REMOVE_DUPLICATES=true && samtools index {output[0]}"
#for newer picardtools version "picard MarkDuplicates -I {input[0]} -O {output[0]} -M {output[1]} -ASO coordinate -REMOVE_SEQUENCING_DUPLICATES true -REMOVE_DUPLICATES true && samtools index {output[0]}"

rule remove_dups_unique:
    input:  "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam",
            "QC/{file}_mapped_sorted_unique_fastqc.zip"
    output: report("UNIQUE_MAPPED/{file}_mapped_sorted_unique_dedup.bam", category="UNIQUE"),
            "LOGS/{file}_unique_dedup.txt",
            "UNIQUE_MAPPED/{file}_mapped_sorted_unique_dedup.bam.bai"
    conda: "../envs/picardtools.yaml"
    threads: 1
    shell:  "picard MarkDuplicates I={input[0]} O={output[0]} M={output[1]} ASO=coordinate REMOVE_SEQUENCING_DUPLICATES=true REMOVE_DUPLICATES=true && samtools index {output[0]}"
#"picard MarkDuplicates -I {input[0]} -O {output[0]} -M {output[1]} -ASO coordinate -REMOVE_SEQUENCING_DUPLICATES true -REMOVE_DUPLICATES true && samtools index {output[0]}"

rule remove_dups_phased:
    input:  "PHASED_MAPPED/{file}_mapped_sorted_phased.bam",
            "QC/{file}_mapped_sorted_phased_fastqc.zip"
    output: report("PHASED_MAPPED/{file}_mapped_sorted_phased_dedup.bam", category="PHASED"),
            "LOGS/{file}_phased_dedup.txt",
            "PHASED_MAPPED/{file}_mapped_sorted_phased_dedup.bam.bai"
    conda: "../envs/picardtools.yaml"
    threads: 1
    shell:  "picard MarkDuplicates I={input[0]} O={output[0]} M={output[1]} ASO=coordinate REMOVE_SEQUENCING_DUPLICATES=true REMOVE_DUPLICATES=true && samtools index {output[0]}"
#"picard MarkDuplicates -I {input[0]} -O {output[0]} -M {output[1]} -ASO coordinate -REMOVE_SEQUENCING_DUPLICATES true -REMOVE_DUPLICATES true && samtools index {output[0]}"

rule fastqc_dedup:
    input:  "SORTED_MAPPED/{file}_mapped_sorted_dedup.bam"
    output: report("QC/{file}_mapped_sorted_dedup_fastqc.zip", category="DEDUP")
    conda: "../envs/qc.yaml"
    threads: 20
    params: dir= expand("QC/{source}",source=SOURCE)
    shell: "for i in {input}; do OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f bam {input};done"

rule fastqc_uniqdedup:
    input:  "UNIQUE_MAPPED/{file}_mapped_sorted_unique_dedup.bam"
    output: report("QC/{file}_mapped_sorted_unique_dedup_fastqc.zip", category="QC")
    conda: "../envs/qc.yaml"
    threads: 20
    params: dir= expand("QC/{source}", source=SOURCE)
    shell: "for i in {input}; do OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f bam {input};done"

rule fastqc_phaseddedup:
    input:  "PHASED_MAPPED/{file}_mapped_sorted_phased_dedup.bam"
    output: report("QC/{file}_mapped_sorted_phased_dedup_fastqc.zip", category="QC")
    conda: "../envs/qc.yaml"
    threads: 20
    params: dir= expand("QC/{source}", source=SOURCE)
    shell: "for i in {input}; do OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f bam {input};done"

rule multiqc:
    input:  expand("QC/{file}_fastqc.zip",file=SAMPLES),
            expand("QC/{file}_mapped_sorted_fastqc.zip",file=SAMPLES),
            expand("QC/{file}_mapped_sorted_unique_fastqc.zip",file=SAMPLES),
            expand("QC/{file}_mapped_sorted_phased_fastqc.zip",file=SAMPLES),
            expand("QC/{file}_mapped_sorted_dedup_fastqc.zip",file=SAMPLES),
            expand("QC/{file}_mapped_sorted_unique_dedup_fastqc.zip",file=SAMPLES),
            expand("QC/{file}_mapped_sorted_phased_dedup_fastqc.zip",file=SAMPLES)
    output: report("QC/Multi/multiqc_report.html", category="QC")
    conda: "../envs/qc.yaml"
    params: dir=expand("QC/Multi",source=SOURCE)
    shell:  "multiqc -k json -z -o {params.dir} ."

rule themall:
    input:  "QC/Multi/multiqc_report.html", "PHASED_MAPPED/{file}_mapped_sorted_phased_dedup.bam" if "ON" in config["QC"] else "PHASED_MAPPED/{file}_mapped_sorted_phased_dedup.bam"
    output: "DONE/{file}_processed"
#    conda: "../envs/base.yaml"
    run:
        for f in output:
            with open(f, "w") as out:
                out.write("DONE")

#rule star_pe_multi:
#    input:
#        # use a list for multiple fastq files for one sample
#        # usually technical replicates across lanes/flowcells
#        fq1 = ["reads/{sample}_R1.1.fastq", "reads/{sample}_R1.2.fastq"],
#        # paired end reads needs to be ordered so each item in the two lists match
#        fq2 = ["reads/{sample}_R2.1.fastq", "reads/{sample}_R2.2.fastq"] #optional
#    output:
#        # see STAR manual for additional output files
#        "MAPPED/{file}_mapped.sam"
#    log:
#        "logs/star/pe/{sample}.log"
#    params:
#        # path to STAR reference genome index
#        index=expand("{ref}{gen}/{gen}{name}.idx",ref=REFERENCE,gen=GENOME, name=NAME),
#        # optional parameters
#        extra=""
#    threads: 8
#    wrapper:
#        "0.27.1/bio/star/align"
#
#rule star_se:
#    input:
#        fq1 = "reads/{sample}_R1.1.fastq"
#    output:
#        # see STAR manual for additional output files
#        "MAPPED/{file}_mapped.sam"
#    log:
#        "logs/star/{sample}.log"
#    params:
#        # path to STAR reference genome index
#        index=expand("{ref}{gen}/{gen}{name}.idx",ref=REFERENCE,gen=GENOME, name=NAME),
#        # optional parameters
#        extra=""
#    threads: 8
#    wrapper:
#        "0.27.1/bio/star/align"
