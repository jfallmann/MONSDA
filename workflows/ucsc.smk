include: "header.smk"

rule all:
    input:  expand("DONE/{file}_tracks",file=samplecond(SAMPLES,config))

rule bamtobed:
    input:  "SORTED_MAPPED/{file}_mapped_sorted.bam",
            "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
    output: "UCSC/{file}_mapped_sorted.bed.gz",
            "UCSC/{file}_mapped_unique.bed.gz"
    log:    "LOGS/{file}_ucscbamtobed"
    conda:  "../envs/bedtools.yaml"
    threads: 1
    shell:  "bedtools bamtobed -i {input[0]} |gzip > {output[0]} && bedtools bamtobed -i {input[1]} |gzip > {output[1]} &> {logs}"

rule index_fa:
    input:  expand("{ref}/{{org}}/{{gen}}{{name}}.fa",ref=REFERENCE),
    output: expand("{ref}/{{org}}/{{gen}}{{name}}.fa.fai",ref=REFERENCE)
    log:    "LOGS/{org}/{gen}{name}_ucscindexfa"
    conda:  "../envs/samtools.yaml"
    params: bins = BINS
    threads: 1
    shell:  "for i in {input};do {params.bins}/Preprocessing/indexfa.sh $i 2> {log};done"

rule get_chromsize_genomic:
    input: expand("{ref}/{{org}}/{{gen}}{{name}}.fa.fai",ref=REFERENCE)
    output: expand("{ref}/{{org}}/{{gen}}{{name}}.chrom.sizes",ref=REFERENCE)
    log:    "LOGS/{org}/{gen}{name}_ucscgetcrom"
    conda:  "../envs/samtools.yaml"
    threads: 1
    params: bins = BINS
    shell:  "cut -f1,2 {input} |sed 's/^chr//g' > {output} 2> {log}"

rule BedToBedg:
    input:  rules.bamtobed.output
#            lambda wildcards: "{ref}/{gen}{name}.idx".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=''.join(tool_params(wildcards.file, None ,config, 'MAPPING')[2]))
    output: "UCSC/{file}_mapped_sorted.fw.bedg.gz",
            "UCSC/{file}_mapped_sorted.re.bedg.gz",
            "UCSC/{file}_mapped_unique.fw.bedg.gz",
            "UCSC/{file}_mapped_unique.re.bedg.gz"
    log:    "LOGS/{file}_ucscbedtobedgraph"
    conda:  "../envs/perl.yaml"
    threads: 1
    params: out=lambda wildcards: expand("QC/{source}",source=source_from_sample(wildcards.file)),
            bins = BINS,
            sizes = lambda wildcards: "{ref}/{gen}{name}.chrom.sizes".format(ref=REFERENCE,gen=genomepath(wildcards.file,config),name=namefromfile(wildcards.file, config))
    shell:  "perl {params.bins}/Universal/Bed2Bedgraph.pl -f {input[0]} -c {params.sizes} -v on -x {output[0]} -y {output[1]} -a track 2> {log} && perl {params.bins}/Universal/Bed2Bedgraph.pl -f {input[1]} -c {params.sizes} -v on -x {output[2]} -y {output[3]} -a track 2>> {log}"

### This step generates bigwig files for peaks which can then be copied to a web-browsable directory and uploaded to UCSC via the track field
rule BedgToUCSC:
    input:  rules.BedToBedg.output
    output: "UCSC/{file}_mapped_sorted.fw.bw",
            "UCSC/{file}_mapped_sorted.re.bw",
            "UCSC/{file}_mapped_unique.fw.bw",
            "UCSC/{file}_mapped_unique.re.bw",
            temp("UCSC/{file}_fw_tmp"),
            temp("UCSC/{file}_re_tmp"),
            temp("UCSC/{file}_unique_fw_tmp"),
            temp("UCSC/{file}_unique_re_tmp")
    log:    "LOGS/{file}_bedgtoucsc"
    conda:  "../envs/ucsc.yaml"
    threads: 1
    params: sizes = lambda wildcards: "{ref}/{gen}{name}.chrom.sizes".format(ref=REFERENCE,gen=genomepath(wildcards.file,config),name=namefromfile(wildcards.file, config))
    shell:  "export LC_ALL=C; zcat {input[0]} |sort --parallel={threads} -S 25% -T SORTTMP -t$'\t' -k1,1 -k2,2n > {output[4]} 2> {log} && bedGraphToBigWig {output[4]} {params.sizes} {output[0]} 2>> {log} && zcat {input[1]} |sort --parallel={threads} -S 25% -T SORTTMP -t$'\t' -k1,1 -k2,2n > {output[5]} 2>> {log} && bedGraphToBigWig {output[5]} {params.sizes} {output[1]} 2>> {log} && zcat {input[2]} |sort --parallel={threads} -S 25% -T SORTTMP -t$'\t' -k1,1 -k2,2n > {output[6]} 2>> {log} && bedGraphToBigWig {output[6]} {params.sizes} {output[2]} 2>> {log} && zcat {input[3]} |sort --parallel={threads} -S 25% -T SORTTMP -t$'\t' -k1,1 -k2,2n > {output[7]} 2>> {log} && bedGraphToBigWig {output[7]} {params.sizes} {output[3]} 2>> {log}"

rule themall:
    input:  rules.BedgToUCSC.output
    output: "DONE/{file}_tracks"
    run:
        for f in output:
            with open(f, "w") as out:
                        out.write("DONE")
