include: "header.smk"

NAME=config["NAME"]["genomic"]

rule all:
    input:  expand("DONE/{file}_tracks",file=samplecond(SAMPLES,config))

rule bamtobed:
    input:  "SORTED_MAPPED/{file}_mapped_sorted.bam",
            "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
    output: "UCSC/{file}_mapped_sorted.bed",
            "UCSC/{file}_mapped_unique.bed"
    conda:  "../envs/bedtools.yaml"
    shell:  "bedtools bamtobed -i {input[0]} > {output[0]} && bedtools bamtobed -i {input[1]} > {output[1]} "

rule index_fa:
    input:  expand("{ref}/{{org}}/{{gen}}{name}.fa",ref=REFERENCE, name=NAME),
    output: expand("{ref}/{{org}}/{{gen}}{name}.fa.fai",ref=REFERENCE, name=NAME)
    conda:  "../envs/samtools.yaml"
    params: bins = BINS
    threads: 1
    shell:  "for i in {input};do {params.bins}/Preprocessing/indexfa.sh $i 2> {log};done"

rule get_chromsize_genomic:
    input: expand("{ref}/{{org}}/{{gen}}{name}.fa.fai",ref=REFERENCE, name=NAME)
    output: expand("{ref}/{{org}}/{{gen}}{name}.chrom.sizes",ref=REFERENCE, name=NAME)
    conda:  "../envs/samtools.yaml"
    params: bins = BINS
    threads: 1
    shell:  "cut -f1,2 {input} |sed 's/^chr//g' > {output} 2> {log}"

rule BedToBedg:
    input:  rules.bamtobed.output,
            lambda wildcards: "{ref}/{gen}{name}.idx".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=''.join(mapping_params(wildcards.file, None ,config)[1]))
    output: "UCSC/{file}_mapped_sorted.fw.bedg.gz",
            "UCSC/{file}_mapped_sorted.re.bedg.gz",
            "UCSC/{file}_mapped_unique.fw.bedg.gz",
            "UCSC/{file}_mapped_unique.re.bedg.gz"
    conda:  "../envs/perl.yaml"
    params: out=lambda wildcards: expand("QC/{source}",source=source_from_sample(wildcards.file)),
            bins = BINS,
            sizes = lambda wildcards: "{ref}/{gen}{name}.chrom.sizes".format(ref=REFERENCE,gen=genomepath(wildcards.file,config),name=NAME)
    shell:  "perl {params.bins}/Universal/Bed2Bedgraph.pl -f {input[0]} -c {params.sizes} -v on -x {output[0]} -y {output[1]} -a track && perl {params.bins}/Universal/Bed2Bedgraph.pl -f {input[1]} -c {params.sizes} -v on -x {output[2]} -y {output[3]} -a track"

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
    conda:  "../envs/ucsc.yaml"
    params: sizes = lambda wildcards: "{ref}/{gen}{name}.chrom.sizes".format(ref=REFERENCE,gen=genomepath(wildcards.file,config),name=NAME)
    shell:  "export LC_ALL=C; zcat {input[0]} |sort --parallel={threads} -S 25% -T SORTTMP -t$'\t' -k1,1 -k2,2n > {output[4]} && bedGraphToBigWig {output[4]} {params.sizes} {output[0]} && zcat {input[1]} |sort --parallel={threads} -S 25% -T SORTTMP -t$'\t' -k1,1 -k2,2n > {output[5]} && bedGraphToBigWig {output[5]} {params.sizes} {output[1]} && zcat {input[2]} |sort --parallel={threads} -S 25% -T SORTTMP -t$'\t' -k1,1 -k2,2n > {output[6]} && bedGraphToBigWig {output[6]} {params.sizes} {output[2]} && zcat {input[3]} |sort --parallel={threads} -S 25% -T SORTTMP -t$'\t' -k1,1 -k2,2n > {output[7]} && bedGraphToBigWig {output[7]} {params.sizes} {output[3]}"

rule themall:
    input:  rules.BedgToUCSC.output
    output: "DONE/{file}_tracks"
    run:
        for f in output:
            with open(f, "w") as out:
                        out.write("DONE")
