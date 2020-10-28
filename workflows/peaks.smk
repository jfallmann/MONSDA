PEAKBIN, PEAKENV = env_bin_from_config2(SAMPLES,config,'PEAKS')

wildcard_constraints:
    type = "sorted|sorted_unique" if not rundedup else "sorted|unique|sorted_dedup|sorted_unique_dedup"

if ANNOPEAK is not None:
    if not rundedup:
        rule themall:
            input:  expand("PEAKS/{file}_mapped_{type}.bedg.gz",file=samplecond(SAMPLES,config), type=['sorted','sorted_unique']),
                    expand("UCSC/{file}_peak_{type}.fw.bw.trackdone",file=samplecond(SAMPLES,config), type=['sorted','sorted_unique']),
                    expand("UCSC/{file}_peak_{type}.re.bw.trackdone",file=samplecond(SAMPLES,config), type=['sorted','sorted_unique']),
                    expand("UCSC/{file}_peak_{type}.fw.bedg.gz",file=samplecond(SAMPLES,config), type=['sorted','sorted_unique']),
                    expand("UCSC/{file}_peak_{type}.re.bedg.gz",file=samplecond(SAMPLES,config), type=['sorted','sorted_unique']),
                    expand("PEAKS/{file}_peak_{type}.bed.gz",file=samplecond(SAMPLES,config), type=['sorted','sorted_unique']),
                    expand("PEAKS/{file}_prepeak_{type}.bed.gz",file=samplecond(SAMPLES,config), type=['sorted','sorted_unique']),
                    expand("PEAKS/{file}_peak_seq_{type}.bed.gz",file=samplecond(SAMPLES,config), type=['sorted','sorted_unique']),
                    expand("PEAKS/{file}_peak_anno_{type}.bed.gz",file=samplecond(SAMPLES,config), type=['sorted','sorted_unique'])
    else:
        rule themall:
            input:  expand("PEAKS/{file}_mapped_{type}.bedg.gz",file=samplecond(SAMPLES,config), type=['sorted_dedup','sorted_unique_dedup']),
                    expand("UCSC/{file}_peak_{type}.fw.bw.trackdone",file=samplecond(SAMPLES,config), type=['sorted_dedup','sorted_unique_dedup']),
                    expand("UCSC/{file}_peak_{type}.re.bw.trackdone",file=samplecond(SAMPLES,config), type=['sorted_dedup','sorted_unique_dedup']),
                    expand("UCSC/{file}_peak_{type}.fw.bedg.gz",file=samplecond(SAMPLES,config), type=['sorted_dedup','sorted_unique_dedup']),
                    expand("UCSC/{file}_peak_{type}.re.bedg.gz",file=samplecond(SAMPLES,config), type=['sorted_dedup','sorted_unique_dedup']),
                    expand("PEAKS/{file}_peak_{type}.bed.gz",file=samplecond(SAMPLES,config), type=['sorted_dedup','sorted_unique_dedup']),
                    expand("PEAKS/{file}_prepeak_{type}.bed.gz",file=samplecond(SAMPLES,config), type=['sorted_dedup','sorted_unique_dedup']),
                    expand("PEAKS/{file}_peak_seq_{type}.bed.gz",file=samplecond(SAMPLES,config), type=['sorted_dedup','sorted_unique_dedup']),
                    expand("PEAKS/{file}_peak_anno_{type}.bed.gz",file=samplecond(SAMPLES,config), type=['sorted_dedup','sorted_unique_dedup'])

else:
    if not rundedup:
        rule themall:
            input:  expand("PEAKS/{file}_mapped_{type}.bedg.gz",file=samplecond(SAMPLES,config), type=['sorted','sorted_unique']),
                    expand("UCSC/{file}_peak_{type}.fw.bw.trackdone",file=samplecond(SAMPLES,config), type=['sorted','sorted_unique']),
                    expand("UCSC/{file}_peak_{type}.re.bw.trackdone",file=samplecond(SAMPLES,config), type=['sorted','sorted_unique']),
                    expand("UCSC/{file}_peak_{type}.fw.bedg.gz",file=samplecond(SAMPLES,config), type=['sorted','sorted_unique']),
                    expand("UCSC/{file}_peak_{type}.re.bedg.gz",file=samplecond(SAMPLES,config), type=['sorted','sorted_unique']),
                    expand("PEAKS/{file}_peak_{type}.bed.gz",file=samplecond(SAMPLES,config), type=['sorted','sorted_unique']),
                    expand("PEAKS/{file}_prepeak_{type}.bed.gz",file=samplecond(SAMPLES,config), type=['sorted','sorted_unique']),
                    expand("PEAKS/{file}_peak_seq_{type}.bed.gz",file=samplecond(SAMPLES,config), type=['sorted','sorted_unique'])
    else:
        rule themall:
            input:  expand("PEAKS/{file}_mapped_{type}.bedg.gz",file=samplecond(SAMPLES,config), type=['sorted_dedup','sorted_unique_dedup']),
                    expand("UCSC/{file}_peak_{type}.fw.bw.trackdone",file=samplecond(SAMPLES,config), type=['sorted_dedup','sorted_unique_dedup']),
                    expand("UCSC/{file}_peak_{type}.re.bw.trackdone",file=samplecond(SAMPLES,config), type=['sorted_dedup','sorted_unique_dedup']),
                    expand("UCSC/{file}_peak_{type}.fw.bedg.gz",file=samplecond(SAMPLES,config), type=['sorted_dedup','sorted_unique_dedup']),
                    expand("UCSC/{file}_peak_{type}.re.bedg.gz",file=samplecond(SAMPLES,config), type=['sorted_dedup','sorted_unique_dedup']),
                    expand("PEAKS/{file}_peak_{type}.bed.gz",file=samplecond(SAMPLES,config), type=['sorted_dedup','sorted_unique_dedup']),
                    expand("PEAKS/{file}_prepeak_{type}.bed.gz",file=samplecond(SAMPLES,config), type=['sorted_dedup','sorted_unique_dedup']),
                    expand("PEAKS/{file}_peak_seq_{type}.bed.gz",file=samplecond(SAMPLES,config), type=['sorted_dedup','sorted_unique_dedup'])


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
        conda:  "nextsnakes/envs/base.yaml"
        threads: 1
        params: abs = lambda wildcards: os.path.abspath('BED/'+wildcards.file+'_mapped_'+wildcards.type+'.bed.gz')
        shell:  "ln -s {params.abs} {output}"
elif all(checklist2):
    rule BamToBed:
        input:  "UCSC/{file}_mapped_{type}.bed.gz"
        output: "PEAKS/{file}_mapped_{type}.bed.gz"
        log:    "LOGS/PEAKS/linkbed{file}_{type}.log"
        conda:  "nextsnakes/envs/base.yaml"
        threads: 1
        params: abs = lambda wildcards: os.path.abspath('UCSC/'+wildcards.file+'_mapped_'+wildcards.type+'.bed.gz')
        shell:  "ln -s {params.abs} {output}"
else:
    if not stranded or stranded == 'fr':
        rule BamToBed:
            input:  "MAPPED/{file}_mapped_{type}.bam"
            output: "PEAKS/{file}_mapped_{type}.bed.gz"
            log:    "LOGS/PEAKS/bam2bed_{file}_{type}.log"
            threads: 1
            conda:  "nextsnakes/envs/bedtools.yaml"
            shell:  "bedtools bamtobed -split -i {input[0]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/2$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |gzip > {output[0]} 2> {log}"
    elif stranded and stranded == 'rf':
        rule BamToBed:
            input:  "MAPPED/{file}_mapped_{type}.bam"
            output: "PEAKS/{file}_mapped_{type}.bed.gz"
            log:    "LOGS/PEAKS/bam2bed_{file}_{type}.log"
            threads: 1
            conda:  "nextsnakes/envs/bedtools.yaml"
            shell:  "bedtools bamtobed -split -i {input[0]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/1$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |gzip > {output[0]} 2> {log}"

rule index_fa:
    input:  REFERENCE
    output: expand("{ref}.fa.fai",ref=REFERENCE.replace('.fa.gz',''))
    log:    expand("LOGS/PEAKS/{ref}/indexfa.log", ref=REFERENCE.replace('.fa.gz',''))
    conda:  "nextsnakes/envs/samtools.yaml"
    threads: 1
    params: bins = BINS
    shell:  "for i in {input};do {params.bins}/Preprocessing/indexfa.sh $i 2> {log};done"

rule get_chromsize_genomic:
    input:  expand("{ref}.fa.fai", ref=REFERENCE.replace('.fa.gz',''))
    output: expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz',''))
    log:    expand("LOGS/PEAKS/{ref}/chromsize.log", ref=REFERENCE.replace('.fa.gz',''))
    conda:  "nextsnakes/envs/samtools.yaml"
    threads: 1
    params: bins = BINS
    shell:  "cut -f1,2 {input} > {output} 2> {log}"

rule extendbed:
    input:  pks = "PEAKS/{file}_mapped_{type}.bed.gz",
            ref = expand("{ref}.chrom.sizes",ref=REFERENCE.replace('.fa.gz',''))
    output: ext = "PEAKS/{file}_mapped_extended_{type}.bed.gz"
    log:    "LOGS/PEAKS/bam2bed{type}_{file}.log"
    conda:  "nextsnakes/envs/perl.yaml"
    threads: 1
    params: bins = BINS
    shell:  "{params.bins}/Universal/ExtendBed.pl -u 1 -b {input.pks} -o {output.ext} -g {input.ref} 2> {log}"

rule rev_extendbed:
    input:  pks = "PEAKS/{file}_mapped_{type}.bed.gz",
            ref = expand("{ref}.chrom.sizes",ref=REFERENCE.replace('.fa.gz',''))
    output: "PEAKS/{file}_mapped_revtrimmed_{type}.bed.gz"
    log:    "LOGS/PEAKS/bam2bed{type}_{file}.log"
    conda:  "nextsnakes/envs/perl.yaml"
    threads: 1
    params: bins=BINS
    shell:  "{params.bins}/Universal/ExtendBed.pl -d 1 -b {input.pks} -o {output[0]} -g {input.ref}  2> {log}"

checklist = list()
checklist2 = list()
for file in samplecond(SAMPLES,config):
    for type in ['sorted','unique']:
        for orient in ['fw','rw']:
            checklist.append(os.path.isfile(os.path.abspath('UCSC/'+file+'_mapped_'+type+'.'+orient+'.bedg.gz')))
            checklist2.append(os.path.isfile(os.path.abspath('BED/'+file+'_mapped_'+type+'.'+orient+'.bedg.gz')))

if all(checklist) and CLIP not in ['iCLIP', 'revCLIP']:
    rule BedToBedg:
        input:  bed = expand("UCSC/{{file}}_mapped_{{type}}.{orient}.bedg.gz",orient=['fw','rw']),
                ref = REFERENCE
        output: fwd = "PEAKS/{file}_mapped_{type}.fw.bedg.gz",
                rev = "PEAKS/{file}_mapped_{type}.re.bedg.gz",
                concat = "PEAKS/{file}_mapped_{type}.bedg.gz"
        log:    "LOGS/PEAKS/{file}_{type}.ucscbedtobedgraph"
        conda:  "nextsnakes/envs/base.yaml"
        threads: 1
        params: absf = lambda wildcards: os.path.abspath('UCSC/'+wildcards.file+'_mapped_'+wildcards.type+'.fw.bedg.gz'),
                absr = lambda wildcards: os.path.abspath('UCSC/'+wildcards.file+'_mapped_'+wildcards.type+'.re.bedg.gz')
        shell:  "export LC_ALL=C; export LC_COLLATE=C; ln -s {params.absf} {output.fwd} && ln -s {params.absr} {output.rev} && zcat {params.absf} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"+\")'| sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat} && zcat {params.absr} | perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"-\")'| sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip >> {output.concat}  2> {log}"

elif all(checklist2) and CLIP not in ['iCLIP', 'revCLIP']:
    rule BedToBedg:
        input:  bed = expand("BED/{{file}}_mapped_{{type}}.{orient}.bedg.gz",orient=['fw','rw']),
                ref = REFERENCE
        output: fwd = "PEAKS/{file}_mapped_{type}.fw.bedg.gz",
                rev = "PEAKS/{file}_mapped_{type}.re.bedg.gz",
                concat = "PEAKS/{file}_mapped_{type}.bedg.gz"
        log:    "LOGS/PEAKS/{file}_{type}.ucscbedtobedgraph"
        conda:  "nextsnakes/envs/base.yaml"
        threads: 1
        params: absf = lambda wildcards: os.path.abspath('BED/'+wildcards.file+'_mapped_'+wildcards.type+'.fw.bedg.gz'),
                absr = lambda wildcards: os.path.abspath('BED/'+wildcards.file+'_mapped_'+wildcards.type+'.re.bedg.gz')
        shell:  "export LC_ALL=C; export LC_COLLATE=C; ln -s {params.absf} {output.fwd} && ln -s {params.absr} {output.rev} && zcat {params.absf} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"+\")'| sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat} && zcat {params.absr} | perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"-\")'| sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip >> {output.concat}  2> {log}"

elif CLIP == 'iCLIP':
     rule BedToBedg:
        input:  bed = "PEAKS/{file}_mapped_extended_{type}.bed.gz",
                fai = expand("{ref}.fa.fai",ref=REFERENCE.replace('.fa.gz','')),
                sizes = expand("{ref}.chrom.sizes",ref=REFERENCE.replace('.fa.gz',''))
        output: concat = "PEAKS/{file}_mapped_{type}.bedg.gz"
        log:    "LOGS/PEAKS/bed2bedgraph{type}_{file}.log"
        conda:  "nextsnakes/envs/bedtools.yaml"
        threads: 1
        params: bins=BINS,
                odir=lambda wildcards,output:(os.path.dirname(output[0]))
        shell: "export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i {input.bed} -bg -split -strand + -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"+\")'| sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat} 2> {log} && bedtools genomecov -i {input.bed} -bg -split -strand - -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"-\")'|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip >> {output.concat} 2>> {log}"

elif CLIP == 'revCLIP':
    rule BedToBedg:
        input:  bed = "PEAKS/{file}_mapped_revtrimmed_{type}.bed.gz",
                fai = expand("{ref}.fa.fai",ref=REFERENCE.replace('.fa.gz','')),
                sizes = expand("{ref}.chrom.sizes",ref=REFERENCE.replace('.fa.gz',''))
        output: concat = "PEAKS/{file}_mapped_{type}.bedg.gz"
        log:    "LOGS/PEAKS/bed2bedgraph{type}_{file}.log"
        conda:  "nextsnakes/envs/bedtools.yaml"
        threads: 1
        params: bins=BINS,
                odir=lambda wildcards,output:(os.path.dirname(output[0]))
        shell: "export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i {input.bed} -bg -split -strand + -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"+\")'| sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat} 2> {log} && bedtools genomecov -i {input.bed} -bg -split -strand - -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"-\")'|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip >> {output.concat} 2>> {log}"

else:
    rule BedToBedg:
        input:  bed = "PEAKS/{file}_mapped_{type}.bed.gz",
                fai = expand("{ref}.fa.fai",ref=REFERENCE.replace('.fa.gz','')),
                sizes = expand("{ref}.chrom.sizes",ref=REFERENCE.replace('.fa.gz',''))
        output: concat = "PEAKS/{file}_mapped_{type}.bedg.gz"
        log:    "LOGS/PEAKS/bed2bedgraph{type}_{file}.log"
        conda:  "nextsnakes/envs/bedtools.yaml"
        threads: 1
        params: bins=BINS,
                odir=lambda wildcards,output:(os.path.dirname(output[0]))
        shell: "export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i {input.bed} -bg -split -strand + -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"+\")'| sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat} 2> {log} && bedtools genomecov -i {input.bed} -bg -split -strand - -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"-\")'|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip >> {output.concat} 2>> {log}"

rule PreprocessPeaks:
    input:  "PEAKS/{file}_mapped_{type}.bedg.gz"
    output: "PEAKS/{file}_prepeak_{type}.bed.gz",
    log:    "LOGS/PEAKS/prepeak_{type}_{file}.log"
    conda:  "nextsnakes/envs/perl.yaml"
    threads: 1
    params:  bins=BINS,
             opts=lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "PEAKS")['OPTIONS'][0].items()),
    shell:  "perl {params.bins}/Analysis/PreprocessPeaks.pl -p {input[0]} {params.opts} |sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n | gzip > {output[0]} 2> {log}"

rule Find_Peaks:
    input:  "PEAKS/{file}_prepeak_{type}.bed.gz"
    output: "PEAKS/{file}_peak_{type}.bed.gz"
    log:    "LOGS/PEAKS/findpeaks{type}_{file}.log"
    conda:  "nextsnakes/envs/perl.yaml"
    threads: 1
    params: opts=lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "PEAKS")['OPTIONS'][1].items()),
            bins=BINS
    shell:  "perl {params.bins}/Analysis/FindPeaks.pl {params.opts} | sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output[0]} 2> {log}"

#rule QuantPeaks:
#   input:  "PEAKS/{source}/Peak_{file}.bed.gz"
#   output: "PEAKS/{source}/QuantPeak_{file}.bed.gz"
#   params: limit=config["MINPEAKHEIGHT},
#       distance=config["PEAKDISTANCE},
#       width=config["PEAKWIDTH},
#       ratio=config["PEAKCUTOFF}
#   shell:

rule UnzipGenome:
    input:  ref = REFERENCE,
    output: fa = expand("{ref}_fastafrombed.fa",ref=REFERENCE.replace('.fa.gz',''))
    log:    "LOGS/PEAKS/indexfa.log"
    conda:  "nextsnakes/envs/samtools.yaml"
    threads: 1
    params: bins = BINS
    shell:  "zcat {input[0]} |perl -F\\\\040 -wlane 'if($_ =~ /^>/){{($F[0] = $F[0] =~ /^>chr/ ? $F[0] : \">chr\".substr($F[0],1))=~ s/\_/\./g;print $F[0]}}else{{print}}' > {output.fa} && {params.bins}/Preprocessing/indexfa.sh {output.fa} 2> {log}"

rule AddSequenceToPeak:
    input:  pk = "PEAKS/{file}_peak_{type}.bed.gz",
            fa = expand("{ref}_fastafrombed.fa",ref=REFERENCE.replace('.fa.gz',''))
    output: peak = "PEAKS/{file}_peak_seq_{type}.bed.gz",
            pt = temp("PEAKS/{file}_peak_chr_{type}.tmp"),
            ps = temp("PEAKS/{file}_peak_seq_{type}.tmp")
    log:    "LOGS/PEAKS/seq2peaks{type}_{file}.log"
    conda:  "nextsnakes/envs/bedtools.yaml"
    threads: 1
    params: bins=BINS
    shell:  "export LC_ALL=C; zcat {input.pk} | perl -wlane '$F[0] = $F[0] =~ /^chr/ ? $F[0] : \"chr\".$F[0]; print join(\"\\t\",@F[0..5])' > {output.pt} && fastaFromBed -fi {input.fa} -bed {output.pt} -name -tab -s -fullHeader -fo {output.ps} && cut -d$'\t' -f2 {output.ps}|sed 's/t/u/ig'|paste -d$'\t' <(zcat {input.pk}) - |sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip  > {output.peak} 2> {log}"  # NEED TO GET RID OF SPACES AND WHATEVER IN HEADER

if ANNOPEAK is not None:
    rule AnnotatePeak:
        input:  "PEAKS/{file}_peak_seq_{type}.bed.gz"
        output: "PEAKS/{file}_peak_anno_{type}.bed.gz"
        log:    "LOGS/PEAKS/annotatepeaks{type}_{file}.log"
        conda:  "nextsnakes/envs/perl.yaml"
        threads: 1
        params: bins=BINS,
                anno = ANNOTATION
        shell:  "perl {params.bins}/Universal/AnnotateBed.pl -b {input} -a {params.anno} |gzip > {output} 2> {log}"

    rule PeakToBedg:
        input:  pk = "PEAKS/{file}_peak_{type}.bed.gz",
                pa = rules.AnnotatePeak.output
        output: fw = "UCSC/{file}_peak_{type}.fw.bedg.gz",
                re = "UCSC/{file}_peak_{type}.re.bedg.gz",
                tfw = temp("UCSC/{file}_peak_{type}.fw.tmp.gz"),
                tre = temp("UCSC/{file}_peak_{type}.re.tmp.gz"),
        log:    "LOGS/PEAKS/peak2bedg{file}_{type}.log"
        conda:  "nextsnakes/envs/perl.yaml"
        threads: 1
        params: bins=BINS,
                sizes = expand("{ref}.chrom.sizes",ref=REFERENCE.replace('.fa.gz',''))
        shell:  "perl {params.bins}/Universal/Bed2Bedgraph.pl -f {input.pk} -c {params.sizes} -p peak -x {output.tfw} -y {output.tre} -a track 2>> {log} && zcat {output.tfw}|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n  |gzip > {output.fw} 2>> {log} &&  zcat {output.tre}|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.re} 2>> {log}"

else:
    rule PeakToBedg:
        input:  pk = "PEAKS/{file}_peak_{type}.bed.gz"
        output: fw = "UCSC/{file}_peak_{type}.fw.bedg.gz",
                re = "UCSC/{file}_peak_{type}.re.bedg.gz",
                tfw = temp("UCSC/{file}_peak_{type}.fw.tmp.gz"),
                tre = temp("UCSC/{file}_peak_{type}.re.tmp.gz"),
        log:    "LOGS/PEAKS/peak2bedg{file}_{type}.log"
        conda:  "nextsnakes/envs/perl.yaml"
        threads: 1
        params: bins=BINS,
                sizes = expand("{ref}.chrom.sizes",ref=REFERENCE.replace('.fa.gz',''))
        shell:  "perl {params.bins}/Universal/Bed2Bedgraph.pl -f {input.pk} -c {params.sizes} -p peak -x {output.tfw} -y {output.tre} -a track 2>> {log} && zcat {output.tfw}|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n  |gzip > {output.fw} 2>> {log} &&  zcat {output.tre}|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.re} 2>> {log}"

### This step normalized the bedg files for comparison in the browser
rule NormalizeBedg:
    input:  fw = rules.PeakToBedg.output.fw,
            re = rules.PeakToBedg.output.re
    output: fw = "UCSC/{file}_peak_{type}.fw.norm.bedg.gz",
            re = "UCSC/{file}_peak_{type}.re.norm.bedg.gz"
    log:    "LOGS/UCSC/{file}_{type}_ucscpeaknormalizebedgraph.log"
    conda:  "nextsnakes/envs/perl.yaml"
    threads: 1
    shell: "export LC_ALL=C; if [[ -n \"$(zcat {input.fw} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then scale=$(bc <<< \"scale=6;1000000/$(zcat {input.fw}|cut -f4|sort -u|wc -l)\") perl -wlane '$sc=$ENV{{scale}};print join(\"\t\",@F[0..$#F-1]),\"\t\",$F[-1]/$sc' <(zcat {input.fw}) |gzip > {output.fw} 2> {log}; else gzip < /dev/null > {output.fw}; echo \"File {input.fw} empty\" >> {log}; fi && if [[ -n \"$(zcat {input.re} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then scale=$(bc <<< \"scale=6;1000000/$(zcat {input.re}|cut -f4|sort -u|wc -l)\") perl -wlane '$sc=$ENV{{scale}};print join(\"\t\",@F[0..$#F-1]),\"\t\",$F[-1]/$sc' <(zcat {input.re})|gzip > {output.re} 2> {log}; else gzip < /dev/null > {output.re}; echo \"File {input.re} empty\" >> {log}; fi"


### This step generates bigwig files for peaks which can then be copied to a web-browsable directory and uploaded to UCSC via the track field
rule PeakToUCSC:
    input:  fw = rules.NormalizeBedg.output.fw,
            re = rules.NormalizeBedg.output.re
    output: fw = "UCSC/{file}_peak_{type}.fw.bw",
            re = "UCSC/{file}_peak_{type}.re.bw",
            tfw = temp("UCSC/{file}_{type}fw_tmp"),
            tre = temp("UCSC/{file}_{type}re_tmp")
    log:    "LOGS/PEAKS/peak2ucsc{file}_{type}.log"
    conda:  "nextsnakes/envs/ucsc.yaml"
    threads: 1
    params: sizes = expand("{ref}.chrom.sizes",ref=REFERENCE.replace('.fa.gz',''))
    shell:  "zcat {input.fw} > {output.tfw} 2>> {log} && bedGraphToBigWig {output.tfw} {params.sizes} {output.fw} 2>> {log} && zcat {input.re} > {output.tre} 2>> {log} && bedGraphToBigWig {output.tre} {params.sizes} {output.re} 2>> {log}"

rule GenerateTrack:
    input:  fw = rules.PeakToUCSC.output.fw,
            re = rules.PeakToUCSC.output.re
    output: "UCSC/{file}_peak_{type}.fw.bw.trackdone",
            "UCSC/{file}_peak_{type}.re.bw.trackdone"
    log:    "LOGS/UCSC/{file}_peak_{type}.log"
    conda:  "nextsnakes/envs/base.yaml"
    threads: MAXTHREAD
    params: bwdir = lambda wildcards: "UCSC/{src}".format(src=SETS),
            bins = os.path.abspath(BINS),
            gen = REFDIR,#lambda wildcards: os.path.basename(genomepath(wildcards.file,config)),
            options = '-n Peaks_'+str(PEAKENV)+' -s peaks -l UCSC_peaks_'+str(PEAKENV)+' -b UCSC_'+str(PEAKENV),
            uid = lambda wildcards: "{src}".format(src='UCSC'+os.sep+"PEAKS_"+SETS.replace(os.sep,'_'))
    shell: "echo -e \"{input.fw}\\n{input.re}\"|python3 {params.bins}/Analysis/GenerateTrackDb.py -i {params.uid} -e 1 -f STDIN -u '' -g {params.gen} {params.options} && touch {input.fw}\.trackdone && touch {input.re}.trackdone 2> {log}"
