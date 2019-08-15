include: "header.smk"

rule all:
    input:  expand("COUNTS/{file}/Counts", file=samplecond(SAMPLES,config)),
            expand("COUNTS/{rawfile}/Counts",rawfile=SAMPLES),
            expand("COUNTS/Featurecounter/{file}_mapped_sorted.counts", file=samplecond(SAMPLES,config)),
            expand("COUNTS/Featurecounter/{file}_mapped_sorted_unique.counts", file=samplecond(SAMPLES,config)),
            "COUNTS/DONE"

if config['MAPPING'] is 'paired':
    rule count_fastq:
        input:  "FASTQ/{rawfile}_r1.fastq.gz",
                "FASTQ/{rawfile}_r2.fastq.gz",
                "TRIMMED_FASTQ/{rawfile}_r1_trimmed.fastq.gz",
                "TRIMMED_FASTQ/{rawfile}_r2_trimmed.fastq.gz"
        output: "COUNTS/{file}_raw_r1_fq.count",
                "COUNTS/{file}_raw_r2_fq.count",
                "COUNTS/{file}_trimmed_r1_fq.count",
                "COUNTS/{file}_trimmed_r2_fq.count"
        conda:  "../envs/base.yaml"
        threads: 1
        shell:  "arr=({input}); alen=${{#arr[@]}}; orr=({output}); for i in \"${{!arr[@]}}\";do a=$(zcat ${{arr[$i]}}|wc -l ); echo $((a/4)) > ${{orr[$i]}};done"

else:
    rule count_fastq:
        input:  "FASTQ/{rawfile}.fastq.gz",
                "TRIMMED_FASTQ/{rawfile}_trimmed.fastq.gz"
        output: "COUNTS/{file}_raw_fq.count",
                "COUNTS/{file}_trimmed_fq.count"
        conda:  "../envs/base.yaml"
        threads: 1
        shell:  "arr=({input}); alen=${{#arr[@]}}; orr=({output}); for i in \"${{!arr[@]}}\";do a=$(zcat ${{arr[$i]}}|wc -l ); echo $((a/4)) > ${{orr[$i]}};done"

rule count_mappers:
    input:  "SORTED_MAPPED/{file}_mapped_sorted.bam",
            "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
    output: "COUNTS/{file}_mapped.count",
            "COUNTS/{file}_mapped_unique.count"
    conda:  "../envs/samtools.yaml"
    threads: MAXTHREAD
    shell:  "export LC_ALL=C; arr=({input}); alen=${{#arr[@]}}; orr=({output}); for i in \"${{!arr[@]}}\";do samtools view -F 260 ${{arr[$i]}} | cut -d$'\t' -f1|sort --parallel={threads} -S 25% -T SORTTMP -u |wc -l > ${{orr[$i]}};done"

rule featurecount:
    input:  "SORTED_MAPPED/{file}_mapped_sorted.bam"
    output: "COUNTS/Featurecounter/{file}_mapped_sorted.counts"
    conda:  "../envs/featurecount.yaml"
    threads: MAXTHREAD
    shell:  "featureCounts -O -M --fraction -T {threads} -t exon -a {ANNOTATION} -o {output[0]} {input[0]}"

rule featurecount_uniq:
    input:  "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam",
            "COUNTS/MAPPING/{file}_readcounts.tsv"
    output: "COUNTS/Featurecounter/{file}_mapped_sorted_unique.counts"
    conda:  "../envs/featurecount.yaml"
    threads: MAXTHREAD
    shell:  "featureCounts -O -T {threads} -t exon -a {ANNOTATION} -o {output[0]} {input[0]}"

rule summarize_counts:
    input:  rules.count_fastq.output,
            rules.count_mappers.output
    output: "COUNTS/{file}/Counts",
            "COUNTS/{rawfile}/Counts"
    conda:  "../envs/base.yaml"
    threads: 1
    params: current = lambda w,input: os.path.dirname(w.input)
    shell:  "arr=({input}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> COUNTS/{params.current}/Counts && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> COUNTS/{params.current}/Counts; else echo '0' >> COUNTS/{params.current}/Counts;fi;done"

onsuccess:
    print("Workflow finished, no error")

rule themall:
    input:  rules.summarize_counts.output,
            rules.featurecount.output,
            rules.featurecount_uniq.output
    output: "COUNTS/DONE"
    conda:  "../envs/base.yaml"
    threads: 1
    params: bins = BINS
    shell:  "touch {output}"

###rnacounter and cufflinks are to be added later
#rule RNAcountReads:
#   input:  "SORTED_MAPPED/{file}_mapped_sorted.bam",
#           "COUNTS/Featurecounter/{file}_mapped_sorted.counts"
#   output: "COUNTS/RNAcounter/{file}_mapped_sorted.counts",
#           "COUNTS/RNAcounter/{file}_mapped_sorted.gene_counts",
#           "COUNTS/RNAcounter/{file}_mapped_sorted.transcript_counts",
#           "COUNTS/RNAcounter/{file}_mapped_sorted.exon_counts",
#           "COUNTS/RNAcounter/{file}_mapped_sorted.intron_counts"
#   shell:  "/usr/bin/rnacounter --nh -n 1 -t genes,transcripts,exons,introns {input[0]} {ANNOTATION} > {output[0]} && grep 'gene' {output[0]} > {output[1]} && grep 'transcript' {output[0]} > {output[2]} && grep 'exon' {output[0]} > {output[3]} && grep 'intron' {output[0]} > {output[4]}"
#
#rule RNAcountReads_uniq:
#   input:  "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam",
#           "COUNTS/Featurecounter/{file}_mapped_sorted_unique.counts"
#   output: "COUNTS/RNAcounter/{file}_mapped_sorted_unique.counts",
#           "COUNTS/RNAcounter/{file}_mapped_sorted_unique.gene_counts",
#           "COUNTS/RNAcounter/{file}_mapped_sorted_unique.transcript_counts",
#           "COUNTS/RNAcounter/{file}_mapped_sorted_unique.exon_counts",
#           "COUNTS/RNAcounter/{file}_mapped_sorted_unique.intron_counts"
#   shell:  "/usr/bin/rnacounter -n 1 -t genes,transcripts,exons,introns {input[0]} {ANNOTATION} > {output[0]} && grep 'gene' {output[0]} > {output[1]} && grep 'transcript' {output[0]} > {output[2]} && grep 'exon' {output[0]} > {output[3]} && grep 'intron' {output[0]} > {output[4]}"
#
#rule cufflinks:
#   input:  "SORTED_MAPPED/{file}_mapped_sorted.bam",
#           "COUNTS/RNAcounter/{file}_mapped_sorted.counts"
#   output: "QUANT/Cufflinks/{file}/transcripts.gtf"
#   log:    "LOGS/Cufflinks/{file}.log"
#   params: out="QUANT/Cufflinks/{file}"
#   threads: 20
#   shell:  "cufflinks -o {params.out} -p {threads} -G {ANNOTATION} {input[0]}"
#
#rule cufflinks_uniq:
#   input:  "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam",
#           "COUNTS/RNAcounter/{file}_mapped_sorted_unique.counts"
#   output: "QUANT/Cufflinks/{file}_unique/transcripts.gtf"
#   log:    "LOGS/Cufflinks/{file}.log"
#   params: out="QUANT/Cufflinks/{file}_unique"
#   threads: 20
#   shell:  "cufflinks -o {params.out} -p {threads} -G {ANNOTATION} {input[0]}"
