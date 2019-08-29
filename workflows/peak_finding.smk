include: "header.smk"

wildcard_constraints:
    type="sorted|unique"

rule all:
    input:  expand("DONE/{file}_peaks_{type}",file=samplecond(SAMPLES,config), type=['sorted','unique'])

checklist = list()
checklist2 = list()
for file in samplecond(SAMPLES,config):
    for type in ['sorted','unique']:
        checklist.append(os.path.isfile(os.path.abspath('BED/'+file+'_mapped_'+type+'.bed.gz')))
        checklist2.append(os.path.isfile(os.path.abspath('UCSC/'+file+'_mapped_'+type+'.bed.gz')))

if all(checklist):
    rule BamToBed:
        input:  "BED/{file}_mapped_{type}.bed.gz"
        output: "PEAKS/{file}_mapped_{type}.bed.gz"
        log:    "LOGS/Peaks/linkbed{file}_{type}.log"
        conda:  "../envs/base.yaml"
        threads: 1
        params: abs = lambda wildcards: os.path.abspath('BED/'+wildcards.file+'_mapped_'+wildcards.type+'.bed.gz')
        shell:  "ln -s {params.abs} {output}"
elif all(checklist2):
    rule BamToBed:
        input:  "UCSC/{file}_mapped_{type}.bed.gz"
        output: "PEAKS/{file}_mapped_{type}.bed.gz"
        log:    "LOGS/Peaks/linkbed{file}_{type}.log"
        conda:  "../envs/base.yaml"
        threads: 1
        params: abs = lambda wildcards: os.path.abspath('UCSC/'+wildcards.file+'_mapped_'+wildcards.type+'.bed.gz')
        shell:  "ln -s {params.abs} {output}"
else:
    rule bamtobed:
        input:  "SORTED_MAPPED/{file}_mapped_sorted.bam",
                "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
        output: "PEAKS/{file}_mapped_sorted.bed.gz",
                "PEAKS/{file}_mapped_sorted_unique.bed.gz"
        log:    "LOGS/Peaks/bam2bed_{file}.log"
        threads: 1
        conda:  "../envs/bedtools.yaml"
        shell:  "bedtools bamtobed -i {input[0]} |gzip > {output[0]} && bedtools bamtobed -i {input[1]} |gzip > {output[1]} "

rule index_fa:
    input:  expand("{ref}/{{org}}/{{gen}}{{name}}.fa.gz",ref=REFERENCE),
    output: expand("{ref}/{{org}}/{{gen}}{{name}}.fa.fai",ref=REFERENCE)
    log:    "LOGS/Peaks/{org}/{gen}{name}/indexfa.log"
    conda:  "../envs/samtools.yaml"
    threads: 1
    params: bins = BINS
    shell:  "for i in {input};do {params.bins}/Preprocessing/indexfa.sh $i 2> {log};done"

rule get_chromsize_genomic:
    input: expand("{ref}/{{org}}/{{gen}}{{name}}.fa.fai",ref=REFERENCE)
    output: expand("{ref}/{{org}}/{{gen}}{{name}}.chrom.sizes",ref=REFERENCE)
    log:    "LOGS/Peaks/{org}/{gen}{name}/chromsize.log"
    conda:  "../envs/samtools.yaml"
    threads: 1
    params: bins = BINS
    shell:  "cut -f1,2 {input} > {output} 2> {log}"

rule extendbed:
    input:  "PEAKS/{file}_mapped_{type}.bed.gz",
            lambda wildcards: "{ref}/{gen}{name}.idx".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=''.join(tool_params(wildcards.file, None ,config, 'MAPPING')[2]))
    output: "PEAKS/{file}_mapped_extended_{type}.bed.gz"
    log:    "LOGS/Peaks/bam2bed{type}_{file}.log"
    conda:  "../envs/perl.yaml"
    threads: 1
    params: gen=lambda wildcards: "{ref}/{gen}{name}.chrom.sizes".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config)),
            bins=BINS
    shell:  "{params.bins}/Universal/ExtendBed.pl -u 1 -b {input[0]} -o {output[0]} -g {params.gen}"

rule bedtobedgraph:
    input:  beds = "PEAKS/{file}_mapped_extended_{type}.bed.gz" if CLIP == 'iCLIP' else "PEAKS/{file}_mapped_{type}.bed.gz",
            sizes = lambda wildcards: "{ref}/{gen}{name}.chrom.sizes".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config))
    output: "PEAKS/{file}_mapped_{type}.bedg.gz"
#            temp("PEAKS/{file}_unsrted{type}.bedg.gz")
    log:    "LOGS/Peaks/bed2bedgraph{type}_{file}.log"
    #    conda:  "../envs/perl.yaml"
    conda:  "../envs/ucsc.yaml"
    threads: 1
    params: genome = lambda wildcards: "{gen}".format(gen=genome(wildcards.file,config)),
            bins=BINS,
            odir=lambda wildcards,output:(os.path.dirname(output[0]))
    shell: "zcat {input.beds}| bedItemOverlapCount {params.genome} -chromSize={input.sizes} stdin |sort -k1,1 -k2,2n|gzip > {output[0]} 2> {log}"
#    shell: "bedtools genomecov -i {input[0]} -dz -split -strand + |perl -wlane 'print join(\"\t\",$F[0],$F[1],$F[1]+1,$F[2])'|gzip > {output.0} 2> {log} && bedtools genomecov -i {input[1]} -dz -split -strand - |perl -wlane 'print join(\"\t\",$F[0],$F[1],$F[1]+1,$F[2])'|gzip > {output.1} 2>> {log} && bedtools genomecov -i {input[2]} -dz -split -strand + |perl -wlane 'print join(\"\t\",$F[0],$F[1],$F[1]+1,$F[2])'|gzip > {output.2} 2>> {log} && bedtools genomecov -i {input[3]} -dz -split -strand - |perl -wlane 'print join(\"\t\",$F[0],$F[1],$F[1]+1,$F[2])'|gzip > {output.3} 2>> {log}"
#    shell:  "{params.bins}/Universal/Bed2Bedgraph.pl -v on -c {params.sizes} -o {params.odir} -f {input[0]} -x {output[1]} -y {output[1]} && zcat {output[1]} |sort -k1,1 -k2,2n|gzip > {output[0]}"

rule PreprocessPeaks:
    input:  "PEAKS/{file}_mapped_{type}.bedg.gz"
    output: "PEAKS/{file}_prepeak_{type}.bed.gz",
    log:    "LOGS/Peaks/prepeak_{type}_{file}.log"
    conda:  "../envs/perl.yaml"
    threads: 1
    params:  bins=BINS
    shell:  "perl {params.bins}/Analysis/PreprocessPeaks.pl -p {input[0]} |sort -k1,1 -k2,2n | gzip > {output[0]}"

rule Find_Peaks:
    input:  "PEAKS/{file}_prepeak_{type}.bed.gz"
    output: "PEAKS/{file}_peak_{type}.bed.gz"
    log:    "LOGS/Peaks/findpeaks{type}_{file}.log"
    conda:  "../envs/perl.yaml"
    threads: 1
    params: limitratio=config["MINPEAKRATIO"],
            ratio=config["PEAKCUTOFF"],
            distance=config["PEAKDISTANCE"],
            width=config["PEAKWIDTH"],
            cutoff=config["MINPEAKHEIGHT"],
            userlimit=config["USRLIMIT"],
            bins=BINS
    shell:  "perl {params.bins}/Analysis/FindPeaks.pl -p {input[0]} -r {params.ratio} -l {params.limitratio} -t {params.distance} -w {params.width} -c {params.cutoff} -a {params.userlimit} | sort -k1,1 -k2,2n  |gzip > {output[0]}"

#rule QuantPeaks:
#   input:  "PEAKS/{source}/Peak_{file}.bed.gz"
#   output: "PEAKS/{source}/QuantPeak_{file}.bed.gz"
#   params: limit=config["MINPEAKHEIGHT},
#       distance=config["PEAKDISTANCE},
#       width=config["PEAKWIDTH},
#       ratio=config["PEAKCUTOFF}
#   shell:

rule AddSequenceToPeak:
    input:  "PEAKS/{file}_peak_{type}.bed.gz"
    output: "PEAKS/{file}_peak_seq_{type}.bed.gz",
            temp("PEAKS/{file}_peak_seq_{type}.tmp")
    log:    "LOGS/Peaks/seq2peaks{type}_{file}.log"
    conda:  "../envs/bedtools.yaml"
    threads: 1
    params: fasta = lambda wildcards: "{ref}/{gen}{name}.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config)),
    shell:  "fastaFromBed -fi {params.fasta} -bed {input[0]} -name+ -tab -s -fullHeader -fo {output[1]} && cut -d$'\t' -f2 {output[1]}|sed 's/t/u/ig'|paste -d$'\t' <(zcat {input[0]}) -|gzip  > {output[0]}"

rule AnnotatePeak:
    input:  "PEAKS/{file}_peak_seq_{type}.bed.gz"
    output: "PEAKS/{file}_peak_anno_{type}.bed.gz"
    log:    "LOGS/Peaks/seq2peaks{type}_{file}.log"
    conda:  "../envs/perl.yaml"
    threads: 1
    params: fasta = lambda wildcards: "{ref}/{gen}{name}.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config)),
            bins=BINS,
            anno=lambda wildcards: anno_from_file(wildcards.file, config, 'annotation')
    shell:  "perl {params.bins}/Universal/AnnotateBed.pl -b {input} -a {params.anno} |gzip > {output}"

rule PeakToBedg:
    input:  "PEAKS/{file}_peak_{type}.bed.gz"
    output: "UCSC/{file}_peak_{type}.fw.bedg.gz",
            "UCSC/{file}_peak_{type}.re.bedg.gz",
            temp("UCSC/{file}_peak_{type}.fw.tmp.gz"),
            temp("UCSC/{file}_peak_{type}.re.tmp.gz"),
    log:    "LOGS/Peaks/peak2bedg{file}_{type}.log"
    conda:  "../envs/perl.yaml"
    threads: 1
    params: out=expand("UCSC/{source}",source=SOURCE),
            bins=BINS,
            sizes = lambda wildcards: "{ref}/{gen}{name}.chrom.sizes".format(ref=REFERENCE,gen=genomepath(wildcards.file,config),name=namefromfile(wildcards.file, config))
    shell:  "perl {params.bins}/Universal/Bed2Bedgraph.pl -f {input[0]} -c {params.sizes} -v on -p peak -x {output[2]} -y {output[3]} -a track 2>> {log} && zcat {output[2]}|sort -k1,1 -k2,2n |gzip > {output[0]} 2>> {log} &&  zcat {output[2]}|sort -k1,1 -k2,2n |gzip > {output[1]} 2>> {log}"

#rule QuantPeakToBedg:
#   input:  "PEAKS/{source}/QuantPeak_{file}.bed.gz"
#   output: "UCSC/{source}/QuantPeak_{file}.fw.bedg.gz",
#       "UCSC/{source}/QuantPeak_{file}.re.bedg.gz"
#   params: out="UCSC/"{source},
#       source=QuantPeak_{file}
#   shell:  "perl {BINS}/Bed2Bedgraph.pl -f {input[0]} -t {params.source} -c /scratch2/fall/Data/GenomeStuff/ChromSizes/hg38.chrom.size -v on -p peak -x {params.out}"

### This step generates bigwig files for peaks which can then be copied to a web-browsable directory and uploaded to UCSC via the track field
rule PeakToUCSC:
    input:  "UCSC/{file}_peak_{type}.fw.bedg.gz",
            "UCSC/{file}_peak_{type}.re.bedg.gz"
    output: "UCSC/{file}_peak_{type}.fw.bw",
            "UCSC/{file}_peak_{type}.re.bw",
            temp("UCSC/{file}_{type}fw_tmp"),
            temp("UCSC/{file}_{type}re_tmp")
    conda:  "../envs/ucsc.yaml"
    threads: 1
    params: sizes = lambda wildcards: "{ref}/{gen}{name}.chrom.sizes".format(ref=REFERENCE,gen=genomepath(wildcards.file,config),name=namefromfile(wildcards.file, config))
    shell:  "zcat {input[0]} |sort -k1,1 -k2,2n > {output[2]} && bedGraphToBigWig {output[2]} {params.sizes} {output[0]} && zcat {input[1]} |sort -k1,1 -k2,2n > {output[3]} && bedGraphToBigWig {output[3]} {params.sizes} {output[1]}"

#rule QuantPeakToUCSC:
#   input:  "UCSC/{source}/Peak_{file}.fw.bedg.gz",
#       "UCSC/{source}/Peak_{file}.re.bedg.gz"
#   output: "UCSC/{source}/Peak_{file}.fw.bw",
#       "UCSC/{source}/Peak_{file}.re.bw"
#   params: out="UCSC/"{source},
#       source=Peak_{file}
#   shell:  "gunzip -c {input[0]} > tmp && bedGraphToBigWig tmp {params.ref}/{params.gen}/{params.gen}.chrom.sizes {params.source}.fw.bw && gunzip -c {input[1]} > tmp && bedGraphToBigWig tmp {params.ref}/{params.gen}/{params.gen}.chrom.sizes {params.source}.re.bw && rm -f tmp"

rule themall:
    input:  "PEAKS/{file}_mapped_{type}.bedg.gz",
            "UCSC/{file}_peak_{type}.fw.bw",
            "UCSC/{file}_peak_{type}.re.bw",
            "UCSC/{file}_peak_{type}.fw.bedg.gz",
            "UCSC/{file}_peak_{type}.re.bedg.gz",
            "PEAKS/{file}_peak_{type}.bed.gz",
            "PEAKS/{file}_prepeak_{type}.bed.gz",
            "PEAKS/{file}_peak_seq_{type}.bed.gz",
            "PEAKS/{file}_peak_anno_{type}.bed.gz" if config["ANNOTATE"] == "ON" else
            "PEAKS/{file}_mapped_{type}.bedg.gz",
            "UCSC/{file}_peak_{type}.fw.bw",
            "UCSC/{file}_peak_{type}.re.bw",
            "UCSC/{file}_peak_{type}.fw.bedg.gz",
            "UCSC/{file}_peak_{type}.re.bedg.gz",
            "PEAKS/{file}_peak_{type}.bed.gz",
            "PEAKS/{file}_prepeak_{type}.bed.gz",
            "PEAKS/{file}_peak_seq_{type}.bed.gz"
    output: "DONE/{file}_peaks_{type}"
    run:
        for f in output:
            with open(f, "w") as out:
                        out.write("DONE")
onsuccess:
    print("Workflow finished, no error")
