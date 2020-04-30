wildcard_constraints:
    type="sorted|unique"
if ANNOPEAK is not None:
    rule themall:
        input:  expand("PEAKS/{file}_mapped_{type}.bedg.gz",file=samplecond(SAMPLES,config), type=['sorted','unique']),
                expand("UCSC/{file}_peak_{type}.fw.bw",file=samplecond(SAMPLES,config), type=['sorted','unique']),
                expand("UCSC/{file}_peak_{type}.re.bw",file=samplecond(SAMPLES,config), type=['sorted','unique']),
                expand("UCSC/{file}_peak_{type}.fw.bedg.gz",file=samplecond(SAMPLES,config), type=['sorted','unique']),
                expand("UCSC/{file}_peak_{type}.re.bedg.gz",file=samplecond(SAMPLES,config), type=['sorted','unique']),
                expand("PEAKS/{file}_peak_{type}.bed.gz",file=samplecond(SAMPLES,config), type=['sorted','unique']),
                expand("PEAKS/{file}_prepeak_{type}.bed.gz",file=samplecond(SAMPLES,config), type=['sorted','unique']),
                expand("PEAKS/{file}_peak_seq_{type}.bed.gz",file=samplecond(SAMPLES,config), type=['sorted','unique']),
                expand("PEAKS/{file}_peak_anno_{type}.bed.gz",file=samplecond(SAMPLES,config), type=['sorted','unique'])
else:
    rule themall:
        input:  expand("PEAKS/{file}_mapped_{type}.bedg.gz",file=samplecond(SAMPLES,config), type=['sorted','unique']),
                expand("UCSC/{file}_peak_{type}.fw.bw",file=samplecond(SAMPLES,config), type=['sorted','unique']),
                expand("UCSC/{file}_peak_{type}.re.bw",file=samplecond(SAMPLES,config), type=['sorted','unique']),
                expand("UCSC/{file}_peak_{type}.fw.bedg.gz",file=samplecond(SAMPLES,config), type=['sorted','unique']),
                expand("UCSC/{file}_peak_{type}.re.bedg.gz",file=samplecond(SAMPLES,config), type=['sorted','unique']),
                expand("PEAKS/{file}_peak_{type}.bed.gz",file=samplecond(SAMPLES,config), type=['sorted','unique']),
                expand("PEAKS/{file}_prepeak_{type}.bed.gz",file=samplecond(SAMPLES,config), type=['sorted','unique']),
                expand("PEAKS/{file}_peak_seq_{type}.bed.gz",file=samplecond(SAMPLES,config), type=['sorted','unique'])

checklist = list()
checklist2 = list()
for file in samplecond(SAMPLES,config):
    for type in ['sorted','unique']:
        checklist.append(os.path.isfile(os.path.abspath('BED/'+file+'_mapped_'+type+'.bed.gz')) and not os.path.islink(os.path.abspath('BED/'+file+'_mapped_'+type+'.bed.gz')))
        checklist2.append(os.path.isfile(os.path.abspath('UCSC/'+file+'_mapped_'+type+'.bed.gz')) and not os.path.islink(os.path.abspath('UCSC/'+file+'_mapped_'+type+'.bed.gz')))

if all(checklist):
    rule BamToBed:
        input:  "BED/{file}_mapped_{type}.bed.gz"
        output: "PEAKS/{file}_mapped_{type}.bed.gz"
        log:    "LOGS/PEAKS/linkbed{file}_{type}.log"
        conda:  "snakes/envs/base.yaml"
        threads: 1
        params: abs = lambda wildcards: os.path.abspath('BED/'+wildcards.file+'_mapped_'+wildcards.type+'.bed.gz')
        shell:  "ln -s {params.abs} {output}"
elif all(checklist2):
    rule BamToBed:
        input:  "UCSC/{file}_mapped_{type}.bed.gz"
        output: "PEAKS/{file}_mapped_{type}.bed.gz"
        log:    "LOGS/PEAKS/linkbed{file}_{type}.log"
        conda:  "snakes/envs/base.yaml"
        threads: 1
        params: abs = lambda wildcards: os.path.abspath('UCSC/'+wildcards.file+'_mapped_'+wildcards.type+'.bed.gz')
        shell:  "ln -s {params.abs} {output}"
else:
    if not stranded or stranded == 'fr':
        rule BamToBed:
            input:  "MAPPED/{file}_mapped_sorted.bam",
                    "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
            output: "PEAKS/{file}_mapped_sorted.bed.gz",
                    "PEAKS/{file}_mapped_unique.bed.gz"
            log:    "LOGS/PEAKS/bam2bed_{file}.log"
            threads: 1
            conda:  "snakes/envs/bedtools.yaml"
            shell:  "bedtools bamtobed -split -i {input[0]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/2$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |gzip > {output[0]} 2> {log} && bedtools bamtobed -split -i {input[1]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/2$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |gzip > {output[1]} 2>> {log}"
    elif stranded and stranded == 'rf':
        rule BamToBed:
            input:  "MAPPED/{file}_mapped_sorted.bam",
                    "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
            output: "PEAKS/{file}_mapped_sorted.bed.gz",
                    "PEAKS/{file}_mapped_unique.bed.gz"
            log:    "LOGS/PEAKS/bam2bed_{file}.log"
            threads: 1
            conda:  "snakes/envs/bedtools.yaml"
            shell:  "bedtools bamtobed -split -i {input[0]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/1$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |gzip > {output[0]} 2> {log} && bedtools bamtobed -split -i {input[1]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/1$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |gzip > {output[1]} 2>> {log}"

rule index_fa:
    input:  expand("{ref}/{{org}}/{{gen}}{{name}}.fa.gz",ref=REFERENCE),
    output: expand("{ref}/{{org}}/{{gen}}{{name}}.fa.fai",ref=REFERENCE)
    log:    "LOGS/PEAKS/{org}/{gen}{name}/indexfa.log"
    conda:  "snakes/envs/samtools.yaml"
    threads: 1
    params: bins = BINS
    shell:  "for i in {input};do {params.bins}/Preprocessing/indexfa.sh $i 2> {log};done"

rule get_chromsize_genomic:
    input: expand("{ref}/{{org}}/{{gen}}{{name}}.fa.fai",ref=REFERENCE)
    output: expand("{ref}/{{org}}/{{gen}}{{name}}.chrom.sizes",ref=REFERENCE)
    log:    "LOGS/PEAKS/{org}/{gen}{name}/chromsize.log"
    conda:  "snakes/envs/samtools.yaml"
    threads: 1
    params: bins = BINS
    shell:  "cut -f1,2 {input} > {output} 2> {log}"

rule extendbed:
    input:  "PEAKS/{file}_mapped_{type}.bed.gz",
            lambda wildcards: "{ref}/{gen}{name}.fa.fai".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config))#''.join(tool_params(wildcards.file, None ,config, 'PEAKS')[2]))
    output: "PEAKS/{file}_mapped_extended_{type}.bed.gz"
    log:    "LOGS/PEAKS/bam2bed{type}_{file}.log"
    conda:  "snakes/envs/perl.yaml"
    threads: 1
    params: gen=lambda wildcards: "{ref}/{gen}{name}.chrom.sizes".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config)),
            bins=BINS
    shell:  "{params.bins}/Universal/ExtendBed.pl -u 1 -b {input[0]} -o {output[0]} -g {params.gen} 2> {log}"

rule rev_extendbed:
    input:  "PEAKS/{file}_mapped_{type}.bed.gz",
            lambda wildcards: "{ref}/{gen}{name}.fa.fai".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config))
    output: "PEAKS/{file}_mapped_revtrimmed_{type}.bed.gz"
    log:    "LOGS/PEAKS/bam2bed{type}_{file}.log"
    conda:  "snakes/envs/perl.yaml"
    threads: 1
    params: gen=lambda wildcards: "{ref}/{gen}{name}.chrom.sizes".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config)),
            bins=BINS
    shell:  "{params.bins}/Universal/ExtendBed.pl -d 1 -b {input[0]} -o {output[0]} -g {params.gen}  2> {log}"

checklist = list()
checklist2 = list()
for file in samplecond(SAMPLES,config):
    for type in ['sorted','unique']:
        for orient in ['fw','rw']:
            checklist.append(os.path.isfile(os.path.abspath('UCSC/'+file+'_mapped_'+type+'.'+orient+'.bedg.gz')))
            checklist2.append(os.path.isfile(os.path.abspath('BED/'+file+'_mapped_'+type+'.'+orient+'.bedg.gz')))

if all(checklist) and CLIP not in ['iCLIP', 'revCLIP']:
    rule BedToBedg:
        input:  expand("UCSC/{{file}}_mapped_{{type}}.{orient}.bedg.gz",orient=['fw','rw']),
                lambda wildcards: "{ref}/{gen}{name}.fa.fai".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config))
        output: fwd = "PEAKS/{file}_mapped_{type}.fw.bedg.gz",
                rev = "PEAKS/{file}_mapped_{type}.re.bedg.gz",
                concat = "PEAKS/{file}_mapped_{type}.bedg.gz"
        log:    "LOGS/PEAKS/{file}_{type}.ucscbedtobedgraph"
        conda:  "snakes/envs/base.yaml"
        threads: 1
        params: absf = lambda wildcards: os.path.abspath('UCSC/'+wildcards.file+'_mapped_'+wildcards.type+'.fw.bedg.gz'),
                absr = lambda wildcards: os.path.abspath('UCSC/'+wildcards.file+'_mapped_'+wildcards.type+'.re.bedg.gz')
        shell:  "export LC_ALL=C; export LC_COLLATE=C; ln -s {params.absf} {output.fwd} && ln -s {params.absr} {output.rev} && zcat {params.absf} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"+\")'| sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat} && zcat {params.absr} | perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"-\")'| sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip >> {output.concat}  2> {log}"

elif all(checklist2) and CLIP not in ['iCLIP', 'revCLIP']:
    rule BedToBedg:
        input:  expand("BED/{{file}}_mapped_{{type}}.{orient}.bedg.gz",orient=['fw','rw']),
                lambda wildcards: "{ref}/{gen}{name}.fa.fai".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config))
        output: fwd = "PEAKS/{file}_mapped_{type}.fw.bedg.gz",
                rev = "PEAKS/{file}_mapped_{type}.re.bedg.gz",
                concat = "PEAKS/{file}_mapped_{type}.bedg.gz"
        log:    "LOGS/PEAKS/{file}_{type}.ucscbedtobedgraph"
        conda:  "snakes/envs/base.yaml"
        threads: 1
        params: absf = lambda wildcards: os.path.abspath('BED/'+wildcards.file+'_mapped_'+wildcards.type+'.fw.bedg.gz'),
                absr = lambda wildcards: os.path.abspath('BED/'+wildcards.file+'_mapped_'+wildcards.type+'.re.bedg.gz')
        shell:  "export LC_ALL=C; export LC_COLLATE=C; ln -s {params.absf} {output.fwd} && ln -s {params.absr} {output.rev} && zcat {params.absf} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"+\")'| sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat} && zcat {params.absr} | perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"-\")'| sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip >> {output.concat}  2> {log}"

elif CLIP == 'iCLIP':
     rule BedToBedg:
        input:  bed = "PEAKS/{file}_mapped_extended_{type}.bed.gz",
                fai = lambda wildcards: "{ref}/{gen}{name}.fa.fai".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config)),
                sizes = lambda wildcards: "{ref}/{gen}{name}.chrom.sizes".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config))
        output: concat = "PEAKS/{file}_mapped_{type}.bedg.gz"
        log:    "LOGS/PEAKS/bed2bedgraph{type}_{file}.log"
        conda:  "snakes/envs/bedtools.yaml"
        threads: 1
        params: genome = lambda wildcards: "{gen}".format(gen=genome(wildcards.file,config)),
                bins=BINS,
                odir=lambda wildcards,output:(os.path.dirname(output[0]))
        shell: "export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i {input.bed} -bg -split -strand + -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"+\")'| sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat} 2> {log} && bedtools genomecov -i {input.bed} -bg -split -strand - -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"-\")'|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip >> {output.concat} 2>> {log}"

elif CLIP == 'revCLIP':
    rule BedToBedg:
        input:  bed = "PEAKS/{file}_mapped_revtrimmed_{type}.bed.gz",
                fai = lambda wildcards: "{ref}/{gen}{name}.fa.fai".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config)),
                sizes = lambda wildcards: "{ref}/{gen}{name}.chrom.sizes".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config))
        output: concat = "PEAKS/{file}_mapped_{type}.bedg.gz"
        log:    "LOGS/PEAKS/bed2bedgraph{type}_{file}.log"
        conda:  "snakes/envs/bedtools.yaml"
        threads: 1
        params: genome = lambda wildcards: "{gen}".format(gen=genome(wildcards.file,config)),
                bins=BINS,
                odir=lambda wildcards,output:(os.path.dirname(output[0]))
        shell: "export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i {input.bed} -bg -split -strand + -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"+\")'| sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat} 2> {log} && bedtools genomecov -i {input.bed} -bg -split -strand - -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"-\")'|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip >> {output.concat} 2>> {log}"

else:
    rule BedToBedg:
        input:  bed = "PEAKS/{file}_mapped_{type}.bed.gz",
                fai = lambda wildcards: "{ref}/{gen}{name}.fa.fai".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config)),
                sizes = lambda wildcards: "{ref}/{gen}{name}.chrom.sizes".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config))
        output: concat = "PEAKS/{file}_mapped_{type}.bedg.gz"
        log:    "LOGS/PEAKS/bed2bedgraph{type}_{file}.log"
        conda:  "snakes/envs/bedtools.yaml"
        threads: 1
        params: genome = lambda wildcards: "{gen}".format(gen=genome(wildcards.file,config)),
                bins=BINS,
                odir=lambda wildcards,output:(os.path.dirname(output[0]))
        shell: "export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i {input.bed} -bg -split -strand + -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"+\")'| sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat} 2> {log} && bedtools genomecov -i {input.bed} -bg -split -strand - -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"-\")'|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip >> {output.concat} 2>> {log}"
        #shell: "export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i {input.bed} -bg -split -strand -g {input.sizes} | sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat} 2> {log}"
        #shell: "export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i {input.bed} -bg -split -strand + -g {input.sizes} | sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |perl -wlne 'print join(\"\t\",$_,\".\",\"+\")'|gzip > {output.concat} 2> {log} && bedtools genomecov -i {input.bed} -bg -split -strand - -g {input.sizes} |sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |perl -wlne 'print join(\"\t\",$_,\".\",\"+\")' | gzip >> {output.concat} 2>> {log}"

rule PreprocessPeaks:
    input:  "PEAKS/{file}_mapped_{type}.bedg.gz"
    output: "PEAKS/{file}_prepeak_{type}.bed.gz",
    log:    "LOGS/PEAKS/prepeak_{type}_{file}.log"
    conda:  "snakes/envs/perl.yaml"
    threads: 1
    params:  bins=BINS,
             opts=PREPROCESS
    shell:  "perl {params.bins}/Analysis/PreprocessPeaks.pl -p {input[0]} {params.opts} |sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n | gzip > {output[0]} 2> {log}"

rule Find_Peaks:
    input:  "PEAKS/{file}_prepeak_{type}.bed.gz"
    output: "PEAKS/{file}_peak_{type}.bed.gz"
    log:    "LOGS/PEAKS/findpeaks{type}_{file}.log"
    conda:  "snakes/envs/perl.yaml"
    threads: 1
    params: limitratio=MINPEAKRATIO,
            ratio=PEAKCUTOFF,
            distance=PEAKDISTANCE,
            width=PEAKWIDTH,
            cutoff=MINPEAKHEIGHT,
            userlimit=USRLIMIT,
            bins=BINS
    shell:  "perl {params.bins}/Analysis/FindPeaks.pl -p {input[0]} -r {params.ratio} -l {params.limitratio} -t {params.distance} -w {params.width} -c {params.cutoff} -a {params.userlimit} | sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output[0]} 2> {log}"

#rule QuantPeaks:
#   input:  "PEAKS/{source}/Peak_{file}.bed.gz"
#   output: "PEAKS/{source}/QuantPeak_{file}.bed.gz"
#   params: limit=config["MINPEAKHEIGHT},
#       distance=config["PEAKDISTANCE},
#       width=config["PEAKWIDTH},
#       ratio=config["PEAKCUTOFF}
#   shell:

rule UnzipGenome:
    input:  expand("{ref}/{{org}}/{{gen}}{{name}}.fa.gz",ref=REFERENCE),
    output: fa = expand("{ref}/{{org}}/{{gen}}{{name}}_fastafrombed.fa",ref=REFERENCE)
    log:    "LOGS/PEAKS/{org}/{gen}{name}/indexfa.log"
    conda:  "snakes/envs/samtools.yaml"
    threads: 1
    params: bins = BINS
    shell:  "zcat {input[0]} |perl -F\\\\040 -wlane 'if($_ =~ /^>/){{($F[0] = $F[0] =~ /^>chr/ ? $F[0] : \">chr\".substr($F[0],1))=~ s/\_/\./g;print $F[0]}}else{{print}}' > {output.fa} && {params.bins}/Preprocessing/indexfa.sh {output.fa} 2> {log}"

rule AddSequenceToPeak:
    input:  pk = "PEAKS/{file}_peak_{type}.bed.gz",
            fa = lambda wildcards: "{ref}/{gen}{name}_fastafrombed.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config))
    output: peak = "PEAKS/{file}_peak_seq_{type}.bed.gz",
            pt = temp("PEAKS/{file}_peak_chr_{type}.tmp"),
            ps = temp("PEAKS/{file}_peak_seq_{type}.tmp")
    log:    "LOGS/PEAKS/seq2peaks{type}_{file}.log"
    conda:  "snakes/envs/bedtools.yaml"
    threads: 1
    params: bins=BINS
    shell:  "export LC_ALL=C; zcat {input.pk} | perl -wlane '$F[0] = $F[0] =~ /^chr/ ? $F[0] : \"chr\".$F[0]; print join(\"\\t\",@F[0..5])' > {output.pt} && fastaFromBed -fi {input.fa} -bed {output.pt} -name -tab -s -fullHeader -fo {output.ps} && cut -d$'\t' -f2 {output.ps}|sed 's/t/u/ig'|paste -d$'\t' <(zcat {input.pk}) - |sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip  > {output.peak} 2> {log}"  # NEED TO GET RID OF SPACES AND WHATEVER IN HEADER

if ANNOPEAK is not None:
    rule AnnotatePeak:
        input:  "PEAKS/{file}_peak_seq_{type}.bed.gz"
        output: "PEAKS/{file}_peak_anno_{type}.bed.gz"
        log:    "LOGS/PEAKS/annotatepeaks{type}_{file}.log"
        conda:  "snakes/envs/perl.yaml"
        threads: 1
        params: bins=BINS,
                anno = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'PEAKS')['ANNOTATION']])
        shell:  "perl {params.bins}/Universal/AnnotateBed.pl -b {input} -a {params.anno} |gzip > {output} 2> {log}"

    rule PeakToBedg:
        input:  pk = "PEAKS/{file}_peak_{type}.bed.gz",
                pa = rules.AnnotatePeak.output
        output: "UCSC/{file}_peak_{type}.fw.bedg.gz",
                "UCSC/{file}_peak_{type}.re.bedg.gz",
                temp("UCSC/{file}_peak_{type}.fw.tmp.gz"),
                temp("UCSC/{file}_peak_{type}.re.tmp.gz"),
        log:    "LOGS/PEAKS/peak2bedg{file}_{type}.log"
        conda:  "snakes/envs/perl.yaml"
        threads: 1
        params: out=expand("UCSC/{source}",source=SOURCE),
                bins=BINS,
                sizes = lambda wildcards: "{ref}/{gen}{name}.chrom.sizes".format(ref=REFERENCE,gen=genomepath(wildcards.file,config),name=namefromfile(wildcards.file, config))
        shell:  "perl {params.bins}/Universal/Bed2Bedgraph.pl -f {input.pk} -c {params.sizes} -p peak -x {output[2]} -y {output[3]} -a track 2>> {log} && zcat {output[2]}|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n  |gzip > {output[0]} 2>> {log} &&  zcat {output[2]}|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output[1]} 2>> {log}"

else:
    rule PeakToBedg:
        input:  pk = "PEAKS/{file}_peak_{type}.bed.gz"
        output: "UCSC/{file}_peak_{type}.fw.bedg.gz",
                "UCSC/{file}_peak_{type}.re.bedg.gz",
                temp("UCSC/{file}_peak_{type}.fw.tmp.gz"),
                temp("UCSC/{file}_peak_{type}.re.tmp.gz"),
        log:    "LOGS/PEAKS/peak2bedg{file}_{type}.log"
        conda:  "snakes/envs/perl.yaml"
        threads: 1
        params: out=expand("UCSC/{source}",source=SOURCE),
                bins=BINS,
                sizes = lambda wildcards: "{ref}/{gen}{name}.chrom.sizes".format(ref=REFERENCE,gen=genomepath(wildcards.file,config),name=namefromfile(wildcards.file, config))
        shell:  "perl {params.bins}/Universal/Bed2Bedgraph.pl -f {input.pk} -c {params.sizes} -p peak -x {output[2]} -y {output[3]} -a track 2>> {log} && zcat {output[2]}|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n  |gzip > {output[0]} 2>> {log} &&  zcat {output[2]}|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output[1]} 2>> {log}"

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
    log:    "LOGS/PEAKS/peak2ucsc{file}_{type}.log"
    conda:  "snakes/envs/ucsc.yaml"
    threads: 1
    params: sizes = lambda wildcards: "{ref}/{gen}{name}.chrom.sizes".format(ref=REFERENCE,gen=genomepath(wildcards.file,config),name=namefromfile(wildcards.file, config))
    shell:  "zcat {input[0]} > {output[2]} 2>> {log} && bedGraphToBigWig {output[2]} {params.sizes} {output[0]} 2>> {log} && zcat {input[1]} > {output[3]} 2>> {log} && bedGraphToBigWig {output[3]} {params.sizes} {output[1]} 2>> {log}"

#rule QuantPeakToUCSC:
#   input:  "UCSC/{source}/Peak_{file}.fw.bedg.gz",
#       "UCSC/{source}/Peak_{file}.re.bedg.gz"
#   output: "UCSC/{source}/Peak_{file}.fw.bw",
#       "UCSC/{source}/Peak_{file}.re.bw"
#   params: out="UCSC/"{source},
#       source=Peak_{file}
#   shell:  "gunzip -c {input[0]} > tmp && bedGraphToBigWig tmp {params.ref}/{params.gen}/{params.gen}.chrom.sizes {params.source}.fw.bw && gunzip -c {input[1]} > tmp && bedGraphToBigWig tmp {params.ref}/{params.gen}/{params.gen}.chrom.sizes {params.source}.re.bw && rm -f tmp"

onsuccess:
    print("Workflow finished, no error")
onerror:
	print("ERROR: "+str({log}))
