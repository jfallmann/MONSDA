PEAKBIN, PEAKENV = env_bin_from_config3(config, 'PEAKS')

if ANNOPEAK is not None:
    if not rundedup:
        rule themall:
            input:  expand("UCSC/PEAKS/{combo}{file}_peak_{type}.fw.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("UCSC/PEAKS/{combo}{file}_peak_{type}.re.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}{file}_peak_{type}.fw.bedg.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}{file}_peak_{type}.re.bedg.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}{file}_peak_{type}.bed.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}{file}_peak_seq_{type}.bed.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}{file}_peak_anno_{type}.bed.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique'])
    else:
        rule themall:
            input:  expand("UCSC/PEAKS/{combo}{file}_peak_{type}.fw.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("UCSC/PEAKS/{combo}{file}_peak_{type}.re.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}{file}_peak_{type}.fw.bedg.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}{file}_peak_{type}.re.bedg.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}{file}_peak_{type}.bed.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}{file}_peak_seq_{type}.bed.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}{file}_peak_anno_{type}.bed.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup'])

else:
    if not rundedup:
        rule themall:
            input:  expand("UCSC/PEAKS/{combo}{file}_peak_{type}.fw.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("UCSC/PEAKS/{combo}{file}_peak_{type}.re.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}{file}_peak_{type}.fw.bedg.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}{file}_peak_{type}.re.bedg.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}{file}_peak_{type}.bed.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}{file}_peak_seq_{type}.bed.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique'])
    else:
        rule themall:
            input:  expand("UCSC/PEAKS/{combo}{file}_peak_{type}.fw.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("UCSC/PEAKS/{combo}{file}_peak_{type}.re.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}{file}_peak_{type}.fw.bedg.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}{file}_peak_{type}.re.bedg.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}{file}_peak_{type}.bed.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}{file}_peak_seq_{type}.bed.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup'])


checklist = list()
for file in samplecond(SAMPLES, config):
    checktype = ['sorted', 'unique'] if not rundedup else ['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']
    for type in checktype:
        checklist.append(os.path.isfile(os.path.abspath('BED/'+scombo+file+'_mapped_'+type+'.bed.gz')) and not os.path.islink(os.path.abspath('BED/'+scombo+file+'_mapped_'+type+'.bed.gz')))

if not all(checklist):
    if not stranded or stranded == 'fr':
        rule BamToBed:
            input:  "MAPPED/{scombo}{file}_mapped_{type}.bam"
            output: "BED/{scombo}{file}_mapped_{type}.bed.gz"
            log:    "LOGS/PEAKS/{scombo}{file}bam2bed_{type}.log"
            threads: 1
            conda:  "nextsnakes/envs/bedtools.yaml"
            shell:  "bedtools bamtobed -split -i {input[0]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/2$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |gzip > {output[0]} 2> {log}"

    elif stranded and stranded == 'rf':
        rule BamToBed:
            input:  "MAPPED/{scombo}{file}_mapped_{type}.bam"
            output: "BED/{scombo}{file}_mapped_{type}.bed.gz"
            log:    "LOGS/PEAKS/{scombo}{file}bam2bed_{type}.log"
            threads: 1
            conda:  "nextsnakes/envs/bedtools.yaml"
            shell:  "bedtools bamtobed -split -i {input[0]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/1$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |gzip > {output[0]} 2> {log}"

rule index_fa:
    input:  REFERENCE
    output: expand("{ref}.fa.fai", ref=REFERENCE.replace('.fa.gz', ''))
    log:    expand("LOGS/PEAKS/{combo}{ref}/indexfa.log", ref=REFERENCE.replace('.fa.gz', ''), combo=combo)
    conda:  "nextsnakes/envs/samtools.yaml"
    threads: 1
    params: bins = BINS
    shell:  "for i in {input};do {params.bins}/Preprocessing/indexfa.sh $i 2> {log};done"

rule get_chromsize_genomic:
    input:  expand("{ref}.fa.fai", ref=REFERENCE.replace('.fa.gz', ''))
    output: expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
    log:    expand("LOGS/PEAKS/{combo}{ref}/chromsize.log", ref=REFERENCE.replace('.fa.gz', ''), combo=combo)
    conda:  "nextsnakes/envs/samtools.yaml"
    threads: 1
    params: bins = BINS
    shell:  "cut -f1,2 {input} > {output} 2> {log}"

rule extendbed:
    input:  pks = "BED/{scombo}{file}_mapped_{type}.bed.gz",
            ref = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
    output: ext = "BED/{scombo}{file}_mapped_extended_{type}.bed.gz"
    log:    "LOGS/PEAKS/{scombo}{file}bam2bed_{type}.log"
    conda:  "nextsnakes/envs/perl.yaml"
    threads: 1
    params: bins = BINS
    shell:  "{params.bins}/Universal/ExtendBed.pl -u 1 -b {input.pks} -o {output.ext} -g {input.ref} 2> {log}"

rule rev_extendbed:
    input:  pks = "BED/{scombo}{file}_mapped_{type}.bed.gz",
            ref = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
    output: ext = "BED/{scombo}{file}_mapped_revtrimmed_{type}.bed.gz"
    log:    "LOGS/PEAKS/{scombo}{file}bam2bed_{type}.log"
    conda:  "nextsnakes/envs/perl.yaml"
    threads: 1
    params: bins=BINS
    shell:  "{params.bins}/Universal/ExtendBed.pl -d 1 -b {input.pks} -o {output.ext} -g {input.ref}  2> {log}"

if IP == 'iCLIP':
     rule BedToBedg:
        input:  bed = expand("BED/{combo}{{file}}_mapped_extended_{{type}}.bed.gz", combo=scombo),
                fai = expand("{ref}.fa.fai", ref=REFERENCE.replace('.fa.gz', '')),
                sizes = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
        output: concat = "PEAKS/{combo}{file}_mapped_{type}.bedg.gz"
        log:    "LOGS/PEAKS/{combo}{file}bed2bedgraph_{type}.log"
        conda:  "nextsnakes/envs/bedtools.yaml"
        threads: 1
        params: bins=BINS,
                odir=lambda wildcards, output:(os.path.dirname(output[0]))
        shell: "export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i {input.bed} -bg -split -strand + -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"+\")'| sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat} 2> {log} && bedtools genomecov -i {input.bed} -bg -split -strand - -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"-\")'|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip >> {output.concat} 2>> {log}"

elif IP == 'revCLIP':
    rule BedToBedg:
        input:  bed = expand("BED/{combo}{{file}}_mapped_revtrimmed_{{type}}.bed.gz", combo=scombo),
                fai = expand("{ref}.fa.fai", ref=REFERENCE.replace('.fa.gz', '')),
                sizes = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
        output: concat = "PEAKS/{combo}{file}_mapped_{type}.bedg.gz"
        log:    "LOGS/PEAKS/{combo}bed2bedgraph_{type}_{file}.log"
        conda:  "nextsnakes/envs/bedtools.yaml"
        threads: 1
        params: bins=BINS,
                odir=lambda wildcards, output:(os.path.dirname(output[0]))
        shell: "export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i {input.bed} -bg -split -strand + -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"+\")'| sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat} 2> {log} && bedtools genomecov -i {input.bed} -bg -split -strand - -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"-\")'|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip >> {output.concat} 2>> {log}"

else:
    rule BedToBedg:
        input:  bed = expand("BED/{combo}{{file}}_mapped_{{type}}.bed.gz", combo=scombo),
                fai = expand("{ref}.fa.fai", ref=REFERENCE.replace('.fa.gz', '')),
                sizes = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
        output: concat = "PEAKS/{combo}{file}_mapped_{type}.bedg.gz"
        log:    "LOGS/PEAKS/{combo}bed2bedgraph_{type}_{file}.log"
        conda:  "nextsnakes/envs/bedtools.yaml"
        threads: 1
        params: bins=BINS,
                odir=lambda wildcards, output:(os.path.dirname(output[0]))
        shell: "export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i {input.bed} -bg -split -strand + -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"+\")'| sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat} 2> {log} && bedtools genomecov -i {input.bed} -bg -split -strand - -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"-\")'|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip >> {output.concat} 2>> {log}"

rule PreprocessPeaks:
    input:  bedg = rules.BedToBedg.output.concat
    output: pre = "PEAKS/{combo}{file}_prepeak_{type}.bed.gz",
    log:    "LOGS/PEAKS/{combo}prepeak_{type}_{file}.log"
    conda:  "nextsnakes/envs/perl.yaml"
    threads: 1
    params:  bins=BINS,
             opts=lambda wildcards: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None, config, "PEAKS", PEAKENV)['OPTIONS'][0].items()),
    shell:  "perl {params.bins}/Analysis/PreprocessPeaks.pl -p {input.bedg} {params.opts} |sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n | gzip > {output.pre} 2> {log}"

rule FindPeaks:
    input:  pre = rules.PreprocessPeaks.output.pre
    output: peak = "PEAKS/{combo}{file}_peak_{type}.bed.gz"
    log:    "LOGS/PEAKS/{combo}findpeaks_{type}_{file}.log"
    conda:  "nextsnakes/envs/"+PEAKENV+".yaml"
    threads: 1
    params: ppara = lambda wildcards: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None, config, "PEAKS", PEAKENV)['OPTIONS'][1].items()),
            peak = PEAKBIN
    shell:  "set +o pipefail; export LC_ALL=C; if [[ -n \"$(zcat {input.pre} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then {params.peak} {params.ppara} <(zcat {input.pre}) |tail -n+2| sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.peak} 2> {log}; else gzip < /dev/null > {output.peak}; echo \"File {input.pre} empty\" >> {log}; fi"

rule UnzipGenome:
    input:  ref = REFERENCE,
    output: fa = expand("{ref}_fastafrombed.fa", ref=REFERENCE.replace('.fa.gz', ''))
    log:    expand("LOGS/{combo}indexfa.log", combo=combo)
    conda:  "nextsnakes/envs/samtools.yaml"
    threads: 1
    params: bins = BINS
    shell:  "set +o pipefail; zcat {input[0]} |perl -F\\\\040 -wlane 'if($_ =~ /^>/){{($F[0] = $F[0] =~ /^>chr/ ? $F[0] : \">chr\".substr($F[0],1))=~ s/\_/\./g;print $F[0]}}else{{print}}' > {output.fa} && {params.bins}/Preprocessing/indexfa.sh {output.fa} 2> {log}"

rule AddSequenceToPeak:
    input:  pk = rules.FindPeaks.output.peak,
            fa = expand("{ref}_fastafrombed.fa", ref=REFERENCE.replace('.fa.gz', ''))
    output: peak = "PEAKS/{combo}{file}_peak_seq_{type}.bed.gz",
            pt = temp("PEAKS/{combo}{file}_peak_chr_{type}.tmp"),
            ps = temp("PEAKS/{combo}{file}_peak_seq_{type}.tmp")
    log:    "LOGS/PEAKS/{combo}seq2peaks_{type}_{file}.log"
    conda:  "nextsnakes/envs/bedtools.yaml"
    threads: 1
    params: bins=BINS
    shell:  "set +o pipefail; export LC_ALL=C; if [[ -n \"$(zcat {input.pk} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then export LC_ALL=C; zcat {input.pk} | perl -wlane '$F[0] = $F[0] =~ /^chr/ ? $F[0] : \"chr\".$F[0]; print join(\"\\t\",@F[0..5])' > {output.pt} && fastaFromBed -fi {input.fa} -bed {output.pt} -name -tab -s -fullHeader -fo {output.ps} && cut -d$'\t' -f2 {output.ps}|sed 's/t/u/ig'|paste -d$'\t' <(zcat {input.pk}) - |sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip  > {output.peak} 2> {log}; else gzip < /dev/null > {output.peak} && touch {output.pt} {output.ps}; fi"  # NEED TO GET RID OF SPACES AND WHATEVER IN HEADER

if ANNOPEAK is not None:
    rule AnnotatePeak:
        input:  "PEAKS/{combo}{file}_peak_seq_{type}.bed.gz"
        output: "PEAKS/{combo}{file}_peak_anno_{type}.bed.gz"
        log:    "LOGS/PEAKS/{combo}annotatepeaks_{type}_{file}.log"
        conda:  "nextsnakes/envs/perl.yaml"
        threads: 1
        params: bins=BINS,
                anno = ANNOTATION
        shell:  "set +o pipefail; perl {params.bins}/Universal/AnnotateBed.pl -b {input} -a {params.anno} |gzip > {output} 2> {log}"

    rule PeakToBedg:
        input:  pk = "PEAKS/{combo}{file}_peak_{type}.bed.gz",
                pa = rules.AnnotatePeak.output
        output: fw = "PEAKS/{combo}{file}_peak_{type}.fw.bedg.gz",
                re = "PEAKS/{combo}{file}_peak_{type}.re.bedg.gz",
                tfw = temp("PEAKS/{combo}{file}_peak_{type}.fw.tmp.gz"),
                trw = temp("PEAKS/{combo}{file}_peak_{type}.re.tmp.gz"),
        log:    "LOGS/PEAKS/{combo}peak2bedg_{type}_{file}.log"
        conda:  "nextsnakes/envs/perl.yaml"
        threads: 1
        params: bins=BINS,
                sizes = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
        shell:  "set +o pipefail; perl {params.bins}/Universal/Bed2Bedgraph.pl -f {input.pk} -c {params.sizes} -p peak -x {output.tfw} -y {output.tre} -a track 2>> {log} && zcat {output.tfw}|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n  |gzip > {output.fw} 2>> {log} &&  zcat {output.tre}|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.re} 2>> {log}"

else:
    rule PeakToBedg:
        input:  pk = "PEAKS/{combo}{file}_peak_{type}.bed.gz"
        output: fw = "PEAKS/{combo}{file}_peak_{type}.fw.bedg.gz",
                re = "PEAKS/{combo}{file}_peak_{type}.re.bedg.gz",
                tfw = temp("PEAKS/{combo}{file}_peak_{type}.fw.tmp.gz"),
                tre = temp("PEAKS/{combo}{file}_peak_{type}.re.tmp.gz"),
        log:    "LOGS/PEAKS/{combo}peak2bedg_{type}_{file}.log"
        conda:  "nextsnakes/envs/perl.yaml"
        threads: 1
        params: bins=BINS,
                sizes = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
        shell:  "set +o pipefail; perl {params.bins}/Universal/Bed2Bedgraph.pl -f {input.pk} -c {params.sizes} -p peak -x {output.tfw} -y {output.tre} -a track 2>> {log} && zcat {output.tfw}|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n  |gzip > {output.fw} 2>> {log} &&  zcat {output.tre}|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.re} 2>> {log}"

### This step normalized the bedg files for comparison in the browser
rule NormalizeBedg:
    input:  fw = rules.PeakToBedg.output.fw,
            re = rules.PeakToBedg.output.re
    output: fw = "PEAKS/{combo}{file}_peak_{type}.fw.norm.bedg.gz",
            re = "PEAKS/{combo}{file}_peak_{type}.re.norm.bedg.gz"
    log:    "LOGS/PEAKS/{combo}ucscpeaknormalizebedgraph_{type}_{file}.log"
    conda:  "nextsnakes/envs/perl.yaml"
    threads: 1
    shell: "set +o pipefail; export LC_ALL=C; if [[ -n \"$(zcat {input.fw} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then scale=$(bc <<< \"scale=6;1000000/$(zcat {input.fw}|cut -f4|sort -u|wc -l)\") perl -wlane '$sc=$ENV{{scale}};print join(\"\t\",@F[0..$#F-1]),\"\t\",$F[-1]/$sc' <(zcat {input.fw}) |gzip > {output.fw} 2> {log}; else gzip < /dev/null > {output.fw}; echo \"File {input.fw} empty\" >> {log}; fi && if [[ -n \"$(zcat {input.re} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then scale=$(bc <<< \"scale=6;1000000/$(zcat {input.re}|cut -f4|sort -u|wc -l)\") perl -wlane '$sc=$ENV{{scale}};print join(\"\t\",@F[0..$#F-1]),\"\t\",$F[-1]/$sc' <(zcat {input.re})|gzip > {output.re} 2> {log}; else gzip < /dev/null > {output.re}; echo \"File {input.re} empty\" >> {log}; fi"


### This step generates bigwig files for peaks which can then be copied to a web-browsable directory and uploaded to UCSC via the track field
rule PeakToUCSC:
    input:  fw = rules.NormalizeBedg.output.fw,
            re = rules.NormalizeBedg.output.re
    output: fw = "UCSC/PEAKS/{combo}{file}_peak_{type}.fw.bw",
            re = "UCSC/PEAKS/{combo}{file}_peak_{type}.re.bw",
            tfw = temp("UCSC/PEAKS/{combo}{file}_{type}fw_tmp"),
            tre = temp("UCSC/PEAKS/{combo}{file}_{type}re_tmp")
    log:    "LOGS/PEAKS/{combo}peak2ucsc_{type}_{file}.log"
    conda:  "nextsnakes/envs/ucsc.yaml"
    threads: 1
    params: sizes = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
    shell:  "set +o pipefail; zcat {input.fw} > {output.tfw} 2>> {log} && bedGraphToBigWig {output.tfw} {params.sizes} {output.fw} 2>> {log} && zcat {input.re} > {output.tre} 2>> {log} && bedGraphToBigWig {output.tre} {params.sizes} {output.re} 2>> {log}"

rule GenerateTrack:
    input:  fw = rules.PeakToUCSC.output.fw,
            re = rules.PeakToUCSC.output.re
    output: "UCSC/PEAKS/{combo}{file}_peak_{type}.fw.bw.trackdone",
            "UCSC/PEAKS/{combo}{file}_peak_{type}.re.bw.trackdone"
    log:    "LOGS/PEAKS/{combo}generatetrack_{type}_{file}.log"
    conda:  "nextsnakes/envs/base.yaml"
    threads: MAXTHREAD
    params: bwdir = lambda wildcards: "UCSC/PEAKS/{combo}{src}".format(combo=combo, src=SETS),
            bins = os.path.abspath(BINS),
            gen = REFDIR,#lambda wildcards: os.path.basename(genomepath(wildcards.file, config)),
            options = '-n Peaks_'+str(PEAKENV)+' -s peaks -l UCSC_peaks_'+str(PEAKENV)+' -b UCSC_'+str(PEAKENV),
            uid = lambda wildcards: "{src}".format(src='UCSC'+os.sep+'PEAKS'+os.sep+combo+SETS.replace(os.sep, '_'))
    shell: "set +o pipefail; echo -e \"{input.fw}\\n{input.re}\"|python3 {params.bins}/Analysis/GenerateTrackDb.py -i {params.uid} -e 1 -f STDIN -u '' -g {params.gen} {params.options} && touch {input.fw}\.trackdone && touch {input.re}.trackdone 2> {log}"
