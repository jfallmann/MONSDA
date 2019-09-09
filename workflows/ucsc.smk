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
        # Here I use the strand of the first read in pair as the one determining the strand
        shell:  "bedtools bamtobed -split -i {input[0]} |perl -wlane \'if($F[3]=~/\/2$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])\' |gzip > {output[0]} 2> {log} && bedtools bamtobed -split -i {input[1]} |perl -wlane \'if($F[3]=~/\/2$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])\' |gzip > {output[1]} 2>> {log}"
        #        shell:  "bedtools bamtobed -split -i {input[0]} |gzip > {output[0]} && bedtools bamtobed -split -i {input[1]} |gzip > {output[1]}"
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
    shell:  "cut -f1,2 {input} > {output} 2> {log}"

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
        shell: "export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i {input.bed} -bg -split -strand + -g {input.sizes} | sort --parallel={threads} -S 25% -T SORTTMP -t$'\t' -k1,1 -k2,2n |gzip > {output.fw} 2> {log} && bedtools genomecov -i {input.bed} -bg -split -strand - -g {input.sizes} |sort --parallel={threads} -S 25% -T SORTTMP -t$'\t' -k1,1 -k2,2n |gzip > {output.re} 2>> {log}"
        #        shell:  "awk '{{if($6==\"+\") print}}' {input[0]} | bedItemOverlapCount {params.genome} -chromSize={params.sizes} stdin |sort -k1,1 -k2,2n|gzip > {output[0]} 2> {log} && awk '{{if($6==\"-\") print}}' {input[0]} | bedItemOverlapCount {params.genome} -chromSize={params.sizes} stdin |sort -k1,1 -k2,2n|gzip > {output[1]} 2>> {log} && awk '{{if($6==\"+\") print}}' {input[1]} | bedItemOverlapCount {params.genome} -chromSize={params.sizes} stdin |sort -k1,1 -k2,2n|gzip > {output[2]} 2>> {log} && awk '{{if($6==\"-\") print}}' {input[1]} | bedItemOverlapCount {params.genome} -chromSize={params.sizes} stdin |sort -k1,1 -k2,2n|gzip > {output[3]} 2>> {log}"
        #    shell:  "perl {params.bins}/Universal/Bed2Bedgraph.pl -f {input[0]} -c {params.sizes} -v on -x {output[0]} -y {output[1]} -a track 2> {log} && perl {params.bins}/Universal/Bed2Bedgraph.pl -f {input[1]} -c {params.sizes} -v on -x {output[2]} -y {output[3]} -a track 2>> {log}"

### This step generates bigwig files for bedg which can then be copied to a web-browsable directory and uploaded to UCSC via the track field
rule BedgToUCSC:
    input:  rules.BedToBedg.output,
    output: fw = "UCSC/{file}_mapped_{type}.fw.bw",
            re = "UCSC/{file}_mapped_{type}.re.bw",
            t1 = temp("UCSC/{file}_mapped_{type}.fw.tmp"),
            t2 = temp("UCSC/{file}_mapped_{type}.re.tmp")
    log:    "LOGS/UCSC/{file}_{type}_bedgtoucsc"
    conda:  "../envs/ucsc.yaml"
    threads: 1
    priority: 100               # This should be finished before we generate tracks
    params: sizes = lambda wildcards: "{ref}/{gen}{name}.chrom.sizes".format(ref=REFERENCE,gen=genomepath(wildcards.file,config),name=namefromfile(wildcards.file, config))
    shell:  "zcat {input[0]} > {output.t1} && bedGraphToBigWig {output.t1} {params.sizes} {output.fw} 2> {log} && zcat {input[1]} > {output.t2} && bedGraphToBigWig {output.t2} {params.sizes} {output.re} 2>> {log}"

rule GenerateTrack:
    input:  rules.BedgToUCSC.output
    output: "UCSC/{file}_mapped_{type}.fw.bw.trackdone",
            "UCSC/{file}_mapped_{type}.re.bw.trackdone"
    log:    "LOGS/UCSC/{file}_track_{type}.log"
    conda:  "../envs/base.yaml"
    threads: 1
    params: bwdir = lambda wildcards: "UCSC/{src}".format(src=source_from_sample(wildcards.file)),
            bins = os.path.abspath(BINS),
            gen = lambda wildcards: os.path.basename(genomepath(wildcards.file,config))
    shell: "ls {params.bwdir}/*.bw|python3 {params.bins}/Analysis/GenerateTrackDb.py -e 1 -f STDIN -n AutoHub -s AutoHub -l 'UCSC track AutoGen' -u {params.bwdir} -g {params.gen} -b UCSCHub && for i in {params.bwdir}/*.bw; do touch $i\.trackdone;done 2> {log}"

rule themall:
    input:  rules.GenerateTrack.output
    output: "DONE/UCSC/{file}_{type}_tracks"
    run:
        for f in output:
            with open(f, "w") as out:
                out.write("DONE")
onsuccess:
    print("Workflow finished, no error")
