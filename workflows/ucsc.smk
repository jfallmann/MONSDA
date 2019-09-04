include: "header.smk"

wildcard_constraints:
    type="sorted|unique"

rule all:
    input:  expand("DONE/UCSC/{file}_{type}_tracks",file=samplecond(SAMPLES,config),type=['sorted','unique'])

checklist = list()
checklist2 = list()
for file in samplecond(SAMPLES,config):
    for type in ['sorted','unique']:
        checklist.append(os.path.isfile(os.path.abspath('BED/'+file+'_mapped_'+type+'.bed.gz')))
        checklist2.append(os.path.isfile(os.path.abspath('PEAKS/'+file+'_mapped_'+type+'.bed.gz')))

if all(checklist):
    rule BamToBed:
        input:  "BED/{file}_mapped_{type}.bed.gz"
        output: "UCSC/{file}_mapped_{type}.bed.gz"
        log:    "LOGS/UCSC/linkbed{file}_{type}.log"
        conda:  "../envs/base.yaml"
        threads: 1
        params: abs = lambda wildcards: os.path.abspath('BED/'+wildcards.file+'_mapped_'+wildcards.type+'.bed.gz')
        shell:  "ln -s {params.abs} {output}"
elif all(checklist2):
    rule BamToBed:
        input:  "PEAKS/{file}_mapped_{type}.bed.gz"
        output: "UCSC/{file}_mapped_{type}.bed.gz"
        log:    "LOGS/UCSC/linkbed{file}_{type}.log"
        conda:  "../envs/base.yaml"
        threads: 1
        params: abs = lambda wildcards: os.path.abspath('PEAKS/'+wildcards.file+'_mapped_'+wildcards.type+'.bed.gz')
        shell:  "ln -s {params.abs} {output}"
else:
    rule BamToBed:
        input:  "SORTED_MAPPED/{file}_mapped_sorted.bam",
                "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
        output: "UCSC/{file}_mapped_sorted.bed.gz",
                "UCSC/{file}_mapped_unique.bed.gz"
        log:    "LOGS/UCSC/{file}_ucscbamtobed"
        conda:  "../envs/bedtools.yaml"
        threads: 1
        shell:  "bedtools bamtobed -i {input[0]} |gzip > {output[0]} 2> {log} && bedtools bamtobed -i {input[1]} |gzip > {output[1]} 2>> {log}"

rule index_fa:
    input:  expand("{ref}/{{org}}/{{gen}}{{name}}.fa.gz",ref=REFERENCE),
    output: expand("{ref}/{{org}}/{{gen}}{{name}}.fa.fai",ref=REFERENCE)
    log:    "LOGS/UCSC/{org}/{gen}{name}_ucscindexfa"
    conda:  "../envs/samtools.yaml"
    params: bins = BINS
    threads: 1
    shell:  "for i in {input};do {params.bins}/Preprocessing/indexfa.sh $i 2> {log};done"

rule get_chromsize_genomic:
    input: expand("{ref}/{{org}}/{{gen}}{{name}}.fa.fai",ref=REFERENCE)
    output: expand("{ref}/{{org}}/{{gen}}{{name}}.chrom.sizes",ref=REFERENCE)
    log:    "LOGS/UCSC/{org}/{gen}{name}_ucscgetcrom"
    conda:  "../envs/samtools.yaml"
    threads: 1
    params: bins = BINS
    shell:  "cut -f1,2 {input} |sed 's/^chr//g' > {output} 2> {log}"

checklist = list()
checklist2 = list()
for file in samplecond(SAMPLES,config):
    for type in ['sorted','unique']:
        for orient in ['fw','rw']:
            checklist.append(os.path.isfile(os.path.abspath('BED/'+file+'_mapped_'+type+'.'+orient+'.bedg.gz')))
            checklist2.append(os.path.isfile(os.path.abspath('PEAKS/'+file+'_mapped_'+type+'.'+orient+'.bedg.gz')))

if all(checklist):
    rule BedToBedg:
        input:  "BED/{file}_mapped_{type}.{orient}.bedg.gz",
                sizes = lambda wildcards: "{ref}/{gen}{name}.chrom.sizes".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config))
        output: "UCSC/{file}_mapped_{type}.{orient}.bedg.gz"
        log:    "LOGS/UCSC/{file}_{type}_{orient}_ucscbedtobedgraph"
        conda:  "../envs/ucsc.yaml"
        #    conda:  "../envs/perl.yaml"
        threads: 1
        params: abs = lambda wildcards: os.path.abspath('BED/'+wildcards.file+'_mapped_'+wildcards.type+'.'+wildcards.orient+'.bedg.gz')
        shell:  "ln -s {params.abs} {output}"

elif all(checklist2):
    rule BedToBedg:
        input:  "PEAKS/{file}_mapped_{type}.{orient}.bedg.gz",
                sizes = lambda wildcards: "{ref}/{gen}{name}.chrom.sizes".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config))
        output: "UCSC/{file}_mapped_{type}.{orient}.bedg.gz"
        log:    "LOGS/UCSC/{file}_{type}_{orient}_ucscbedtobedgraph"
        conda:  "../envs/ucsc.yaml"
        #    conda:  "../envs/perl.yaml"
        threads: 1
        params: abs = lambda wildcards: os.path.abspath('PEAKS/'+wildcards.file+'_mapped_'+wildcards.type+'.'+wildcards.orient+'.bedg.gz')
        shell:  "ln -s {params.abs} {output}"

else:
    rule BedToBedg:
        input:  bed = "UCSC/{file}_mapped_{type}.bed.gz",
                fai = lambda wildcards: "{ref}/{gen}{name}.fa.fai".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config)),
                sizes = lambda wildcards: "{ref}/{gen}{name}.chrom.sizes".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config))
        output: fw = "UCSC/{file}_mapped_{type}.fw.bedg.gz",
                re = "UCSC/{file}_mapped_{type}.re.bedg.gz"
        log:    "LOGS/UCSC/{file}_{type}_ucscbedtobedgraph"
        conda:  "../envs/bedtools.yaml"
        #    conda:  "../envs/perl.yaml"
        threads: 1
        params: out=lambda wildcards: expand("QC/{source}",source=source_from_sample(wildcards.file)),
                bins = BINS
        shell: "bedtools genomecov -i {input.bed} -bga -split -strand + -g {input.sizes} |perl -wlane 'print join(\"\t\",$F[0],$F[1],$F[1]+1,$F[2])'|sort -V |gzip > {output.fw} 2> {log} && bedtools genomecov -i {input.bed} -bga -split -strand - -g {input.sizes} |perl -wlane 'print join(\"\t\",$F[0],$F[1],$F[1]+1,$F[2])'|sort -V |gzip > {output.re} 2>> {log}"
        #        shell:  "awk '{{if($6==\"+\") print}}' {input[0]} | bedItemOverlapCount {params.genome} -chromSize={params.sizes} stdin |sort -k1,1 -k2,2n|gzip > {output[0]} 2> {log} && awk '{{if($6==\"-\") print}}' {input[0]} | bedItemOverlapCount {params.genome} -chromSize={params.sizes} stdin |sort -k1,1 -k2,2n|gzip > {output[1]} 2>> {log} && awk '{{if($6==\"+\") print}}' {input[1]} | bedItemOverlapCount {params.genome} -chromSize={params.sizes} stdin |sort -k1,1 -k2,2n|gzip > {output[2]} 2>> {log} && awk '{{if($6==\"-\") print}}' {input[1]} | bedItemOverlapCount {params.genome} -chromSize={params.sizes} stdin |sort -k1,1 -k2,2n|gzip > {output[3]} 2>> {log}"
        #    shell:  "perl {params.bins}/Universal/Bed2Bedgraph.pl -f {input[0]} -c {params.sizes} -v on -x {output[0]} -y {output[1]} -a track 2> {log} && perl {params.bins}/Universal/Bed2Bedgraph.pl -f {input[1]} -c {params.sizes} -v on -x {output[2]} -y {output[3]} -a track 2>> {log}"

### This step generates bigwig files for peaks which can then be copied to a web-browsable directory and uploaded to UCSC via the track field
rule BedgToUCSC:
    input:  rules.BedToBedg.output,
    output: "UCSC/{file}_mapped_{type}.fw.bw",
            "UCSC/{file}_mapped_{type}.re.bw",
            temp("UCSC/{file}_{type}_fw_tmp"),
            temp("UCSC/{file}_{type}_re_tmp")
    log:    "LOGS/UCSC/{file}_{type}_bedgtoucsc"
    conda:  "../envs/ucsc.yaml"
    threads: 1
    params: sizes = lambda wildcards: "{ref}/{gen}{name}.chrom.sizes".format(ref=REFERENCE,gen=genomepath(wildcards.file,config),name=namefromfile(wildcards.file, config))
    shell:  "export LC_ALL=C; zcat {input[0]} |sort --parallel={threads} -S 25% -T SORTTMP -t$'\t' -k1,1 -k2,2n > {output[4]} 2> {log} && bedGraphToBigWig {output[4]} {params.sizes} {output[0]} 2>> {log} && zcat {input[1]} |sort --parallel={threads} -S 25% -T SORTTMP -t$'\t' -k1,1 -k2,2n > {output[5]} 2>> {log} && bedGraphToBigWig {output[5]} {params.sizes} {output[1]} 2>> {log} && zcat {input[2]} |sort --parallel={threads} -S 25% -T SORTTMP -t$'\t' -k1,1 -k2,2n > {output[6]} 2>> {log} && bedGraphToBigWig {output[6]} {params.sizes} {output[2]} 2>> {log} && zcat {input[3]} |sort --parallel={threads} -S 25% -T SORTTMP -t$'\t' -k1,1 -k2,2n > {output[7]} 2>> {log} && bedGraphToBigWig {output[7]} {params.sizes} {output[3]} 2>> {log}"

rule themall:
    input:  rules.BedgToUCSC.output
    output: "DONE/UCSC/{file}_{type}_tracks"
    run:
        for f in output:
            with open(f, "w") as out:
                        out.write("DONE")
onsuccess:
    print("Workflow finished, no error")
