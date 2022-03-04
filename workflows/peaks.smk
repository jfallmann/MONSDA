PEAKBIN, PEAKENV = env_bin_from_config3(config, 'PEAKS')

if ANNOPEAK is not None:
    if not rundedup:
        rule themall:
            input:  expand("TRACKS/PEAKS/{combo}/{file}_peak_{type}.fw.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("TRACKS/PEAKS/{combo}/{file}_peak_{type}.re.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.fw.bedg.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.re.bedg.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.bed.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_seq_{type}.bed.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_anno_{type}.bed.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique'])
    else:
        rule themall:
            input:  expand("TRACKS/PEAKS/{combo}/{file}_peak_{type}.fw.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("TRACKS/PEAKS/{combo}/{file}_peak_{type}.re.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.fw.bedg.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.re.bedg.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.bed.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}/{file}_peak_seq_{type}.bed.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}/{file}_peak_anno_{type}.bed.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup'])

else:
    if not rundedup:
        rule themall:
            input:  expand("TRACKS/PEAKS/{combo}/{file}_peak_{type}.fw.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("TRACKS/PEAKS/{combo}/{file}_peak_{type}.re.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.fw.bedg.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.re.bedg.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.bed.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_seq_{type}.bed.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique'])
    else:
        rule themall:
            input:  expand("TRACKS/PEAKS/{combo}/{file}_peak_{type}.fw.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("TRACKS/PEAKS/{combo}/{file}_peak_{type}.re.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.fw.bedg.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.re.bedg.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.bed.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}/{file}_peak_seq_{type}.bed.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup'])


checklist = list()
for file in samplecond(SAMPLES, config):
    checktype = ['sorted', 'unique'] if not rundedup else ['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']
    for type in checktype:
        checklist.append(os.path.isfile(os.path.abspath('BED/'+scombo+file+'_mapped_'+type+'.bed.gz')) and not os.path.islink(os.path.abspath('BED/'+scombo+file+'_mapped_'+type+'.bed.gz')))

if not all(checklist):
    if not stranded or (stranded == 'fr' or stranded == 'ISF'):
        rule BamToBed:
            input:  expand("MAPPED/{scombo}/{{file}}_mapped_{{type}}_nosoftclip.bam", scombo=scombo)
            output: "BED/{scombo}/{file}_mapped_{type}_nosoftclip.bed.gz"
            log:    "LOGS/PEAKS/{scombo}/{file}bam2bed_{type}.log"
            conda:  "bedtools.yaml"
            threads: 1
            params: sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
            shell:  "bedtools bamtobed -split -i {input[0]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/2$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |sort -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output[0]} 2> {log}"

    elif stranded and (stranded == 'rf' or stranded == 'ISR'):
        rule BamToBed:
            input:  expand("MAPPED/{scombo}/{{file}}_mapped_{{type}}_nosoftclip.bam", scombo=scombo)
            output: "BED/{scombo}/{file}_mapped_{type}_nosoftclip.bed.gz"
            log:    "LOGS/PEAKS/{scombo}/{file}bam2bed_{type}.log"
            conda:  "bedtools.yaml"
            threads: 1
            params: sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
            shell:  "bedtools bamtobed -split -i {input[0]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/1$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |sort -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output[0]} 2> {log}"

include: "manipulate_genome.smk"

rule extendbed:
    input:  pks = "BED/{scombo}/{file}_mapped_{type}.bed.gz",
            ref = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
    output: ext = "BED/{scombo}/{file}_mapped_extended_{type}.bed.gz"
    log:    "LOGS/PEAKS/{scombo}/{file}bam2bed_{type}.log"
    conda:  "perl.yaml"
    threads: 1
    params: bins = BINS
    shell:  "{params.bins}/Universal/ExtendBed.pl -u 1 -b {input.pks} -o {output.ext} -g {input.ref} 2> {log}"

rule rev_extendbed:
    input:  pks = "BED/{scombo}/{file}_mapped_{type}.bed.gz",
            ref = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
    output: ext = "BED/{scombo}/{file}_mapped_revtrimmed_{type}.bed.gz"
    log:    "LOGS/PEAKS/{scombo}/{file}bam2bed_{type}.log"
    conda:  "perl.yaml"
    threads: 1
    params: bins = BINS
    shell:  "{params.bins}/Universal/ExtendBed.pl -d 1 -b {input.pks} -o {output.ext} -g {input.ref}  2> {log}"

if IP == 'iCLIP':
     rule BedToBedg:
        input:  bed = expand("BED/{scombo}/{{file}}_mapped_extended_{{type}}.bed.gz", scombo=scombo),
                fai = expand("{ref}.fa.fai", ref=REFERENCE.replace('.fa.gz', '')),
                sizes = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
        output: concat = "PEAKS/{combo}/{file}_mapped_{type}.bedg.gz",
                tosrt = temp("PEAKS/{combo}/{file}_mapped_{type}.unsrt")
        log:    "LOGS/PEAKS/{combo}/{file}bed2bedgraph_{type}.log"
        conda:  "bedtools.yaml"
        threads: 1
        params: bins = BINS,
                odir = lambda wildcards, output:(os.path.dirname(output[0])),
                sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
        shell: "export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i {input.bed} -bg -split -strand + -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"+\")' > {output.tosrt} 2> {log} && bedtools genomecov -i {input.bed} -bg -split -strand - -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"-\")' >> {output.tosrt} && cat {output.tosrt}| sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat} 2>> {log}"
elif IP == 'revCLIP':
    rule BedToBedg:
        input:  bed = expand("BED/{scombo}/{{file}}_mapped_revtrimmed_{{type}}.bed.gz", scombo=scombo),
                fai = expand("{ref}.fa.fai", ref=REFERENCE.replace('.fa.gz', '')),
                sizes = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
        output: concat = "PEAKS/{combo}/{file}_mapped_{type}.bedg.gz",
                tosrt = temp("PEAKS/{combo}/{file}_mapped_{type}.unsrt")
        log:    "LOGS/PEAKS/{combo}/bed2bedgraph_{type}_{file}.log"
        conda:  "bedtools.yaml"
        threads: 1
        params: bins = BINS,
                odir = lambda wildcards, output:(os.path.dirname(output[0])),
                sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
        shell: "export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i {input.bed} -bg -split -strand + -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"+\")' > {output.tosrt} 2> {log} && bedtools genomecov -i {input.bed} -bg -split -strand - -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"-\")' >> {output.tosrt} && cat {output.tosrt}| sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat} 2>> {log}"
else:
    rule BedToBedg:
        input:  bed = expand("BED/{scombo}/{{file}}_mapped_{{type}}.bed.gz", scombo=scombo),
                fai = expand("{ref}.fa.fai", ref=REFERENCE.replace('.fa.gz', '')),
                sizes = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
        output: concat = "PEAKS/{combo}/{file}_mapped_{type}.bedg.gz",
                tosrt = temp("PEAKS/{combo}/{file}_mapped_{type}.unsrt")
        log:    "LOGS/PEAKS/{combo}/bed2bedgraph_{type}_{file}.log"
        conda:  "bedtools.yaml"
        threads: 1
        params: bins = BINS,
                odir = lambda wildcards, output:(os.path.dirname(output[0])),
                sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
        shell: "export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i {input.bed} -bg -split -strand + -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"+\")' > {output.tosrt} 2> {log} && bedtools genomecov -i {input.bed} -bg -split -strand - -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"-\")' >> {output.tosrt} && cat {output.tosrt}| sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat} 2>> {log}"

rule PreprocessPeaks:
    input:  bedg = rules.BedToBedg.output.concat
    output: pre = "PEAKS/{combo}/{file}_prepeak_{type}.bed.gz"
    log:    "LOGS/PEAKS/{combo}/prepeak_{type}_{file}.log"
    conda:  "perl.yaml"
    threads: 1
    params:  bins = BINS,
             opts = lambda wildcards: tool_params(wildcards.file, None, config, "PEAKS", PEAKENV)['OPTIONS'].get('PREPROCESS', ""),
             sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
    shell:  "perl {params.bins}/Analysis/PreprocessPeaks.pl -p <(zcat {input.bedg}) {params.opts} |sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k3,3n -k2,2n -k6,6 | gzip > {output.pre} 2> {log}"

rule FindPeaks:
    input:  pre = "PEAKS/{combo}/{file}_prepeak_{type}.bed.gz"
    output: peak = "PEAKS/{combo}/{file}_peak_{type}.bed.gz"
    log:    "LOGS/PEAKS/{combo}/{file}findpeaks_{type}.log"
    conda:  "perl.yaml"
    threads: 1
    params: opts = lambda wildcards: tool_params(wildcards.file, None, config, "PEAKS", PEAKENV)['OPTIONS'].get('FINDPEAKS', ""),
            bins = BINS,
            sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
    shell:  "perl {params.bins}/Analysis/FindPeaks.pl {params.opts} -p <(zcat {input.pre}) 2> {log}| sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.peak} 2>> {log}"

#rule QuantPeaks:
#   input:  "{outdir}{source}/Peak_{file}.bed.gz"
#   output: "{outdir}{source}/QuantPeak_{file}.bed.gz"
#   params: limit=config["MINPEAKHEIGHT},
#       distance=config["PEAKDISTANCE},
#       width=config["PEAKWIDTH},
#       ratio=config["PEAKCUTOFF}
#   shell:

rule AddSequenceToPeak:
    input:  pk = rules.FindPeaks.output.peak,
            fa = expand("{ref}.fa", ref=REFERENCE.replace('.fa.gz', '')),
    output: peak = "PEAKS/{combo}/{file}_peak_seq_{type}.bed.gz",
            pt = temp("PEAKS/{combo}/{file}_peak_chr_{type}.tmp"),
            ps = temp("PEAKS/{combo}/{file}_peak_seq_{type}.tmp")
    log:    "LOGS/PEAKS/{combo}/seq2peaks_{type}_{file}.log"
    conda:  "bedtools.yaml"
    threads: 1
    params: bins = BINS,
            sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
    shell:  "set +o pipefail; export LC_ALL=C; if [[ -n \"$(zcat {input.pk} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then export LC_ALL=C; zcat {input.pk} | perl -wlane '$F[0] = $F[0] =~ /^chr/ ? $F[0] : \"chr\".$F[0]; print join(\"\\t\",@F[0..5])' > {output.pt} && bedtools getfasta -fi {input.fa} -bed {output.pt} -name -tab -s -fullHeader -fo {output.ps} && cut -d$'\t' -f2 {output.ps}|paste -d$'\t' <(zcat {input.pk}) - |sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip  > {output.peak} 2> {log}; else gzip < /dev/null > {output.peak} && touch {output.pt} {output.ps}; fi"  # NEED TO GET RID OF SPACES AND WHATEVER IN HEADER

if ANNOPEAK is not None:
    rule AnnotatePeak:
        input:  "PEAKS/{combo}/{file}_peak_seq_{type}.bed.gz"
        output: "PEAKS/{combo}/{file}_peak_anno_{type}.bed.gz"
        log:    "LOGS/PEAKS/{combo}/{file}annotatepeaks_{type}.log"
        conda:  "perl.yaml"
        threads: 1
        params: bins = BINS,
                anno = ANNOTATION
        shell:  "perl {params.bins}/Universal/AnnotateBed.pl -b <(zcat {input}) -a {params.anno} |gzip > {output} 2> {log}"

    rule PeakToBedg:
        input:  pk = "PEAKS/{combo}/{file}_peak_{type}.bed.gz",
                pa = rules.AnnotatePeak.output,
                sizes = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
        output: fw = "PEAKS/{combo}/{file}_peak_{type}.fw.bedg.gz",
                re = "PEAKS/{combo}/{file}_peak_{type}.re.bedg.gz",
                tfw = temp("PEAKS/{combo}/{file}_peak_{type}.fw.tmp.gz"),
                tre = temp("PEAKS/{combo}/{file}_peak_{type}.re.tmp.gz"),
        log:    "LOGS/PEAKS/{combo}/{file}peak2bedg_{type}.log"
        conda:  "perl.yaml"
        threads: 1
        params: bins = BINS,
                sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
        shell:  "perl {params.bins}/Universal/Bed2Bedgraph.pl -f <(zcat {input.pk}) -c {input.sizes} -p peak -x {output.tfw} -y {output.tre} -a track 2>> {log} && zcat {output.tfw}|sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n  |gzip > {output.fw} 2>> {log} &&  zcat {output.tre}|sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.re} 2>> {log}"

else:
    rule PeakToBedg:
        input:  pk = "PEAKS/{combo}/{file}_peak_{type}.bed.gz",
                sizes = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
        output: fw = "PEAKS/{combo}/{file}_peak_{type}.fw.bedg.gz",
                re = "PEAKS/{combo}/{file}_peak_{type}.re.bedg.gz",
                tfw = temp("PEAKS/{combo}/{file}_peak_{type}.fw.tmp.gz"),
                tre = temp("PEAKS/{combo}/{file}_peak_{type}.re.tmp.gz"),
        log:    "LOGS/PEAKS/{combo}/{file}peak2bedg_{type}.log"
        conda:  "perl.yaml"
        threads: 1
        params: bins = BINS,
                sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
        shell:  "perl {params.bins}/Universal/Bed2Bedgraph.pl -f <(zcat {input.pk}) -c {input.sizes} -p peak -x {output.tfw} -y {output.tre} -a track 2>> {log} && zcat {output.tfw}|sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n  |gzip > {output.fw} 2>> {log} &&  zcat {output.tre}|sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.re} 2>> {log}"

### This step normalized the bedg files for comparison in the browser
rule NormalizeBedg:
    input:  fw = rules.PeakToBedg.output.fw,
            re = rules.PeakToBedg.output.re
    output: fw = "PEAKS/{combo}/{file}_peak_{type}.fw.norm.bedg.gz",
            re = "PEAKS/{combo}/{file}_peak_{type}.re.norm.bedg.gz"
    log:    "LOGS/PEAKS/{combo}/{file}ucscpeaknormalizebedgraph_{type}.log"
    conda:  "perl.yaml"
    threads: 1
    params: sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
    shell: "export LC_ALL=C; if [[ -n \"$(zcat {input.fw} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then scale=$(bc <<< \"scale=6;$(zcat {input.fw}|cut -f4|perl -wne '{{$x+=$_;}}END{{if ($x == 0){{$x=1}} print $x}}')/1000000\") perl -wlane '$sc=$ENV{{scale}};print join(\"\t\",@F[0..$#F-1]),\"\t\",$F[-1]/$sc' <(zcat {input.fw})| sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.fw} 2> {log}; else gzip < /dev/null > {output.fw}; echo \"File {input.fw} empty\" >> {log}; fi && if [[ -n \"$(zcat {input.re} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then scale=$(bc <<< \"scale=6;$(zcat {input.re}|cut -f4|perl -wne '{{$x+=$_;}}END{{if ($x == 0){{$x=1}} print $x}}')/1000000\") perl -wlane '$sc=$ENV{{scale}};print join(\"\t\",@F[0..$#F-1]),\"\t\",$F[-1]/$sc' <(zcat {input.re})| sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.re} 2> {log}; else gzip < /dev/null > {output.re}; echo \"File {input.re} empty\" >> {log}; fi"


### This step generates bigwig files for peaks which can then be copied to a web-browsable directory and uploaded to TRACKS via the track field
rule PeakToTRACKS:
    input:  fw = rules.NormalizeBedg.output.fw,
            re = rules.NormalizeBedg.output.re,
            sizes = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
    output: fw = "TRACKS/PEAKS/{combo}/{file}_peak_{type}.fw.bw",
            re = "TRACKS/PEAKS/{combo}/{file}_peak_{type}.re.bw",
            tfw = temp("TRACKS/PEAKS/{combo}/{file}_{type}fw_tmp"),
            tre = temp("TRACKS/PEAKS/{combo}/{file}_{type}re_tmp")
    log:    "LOGS/PEAKS/{combo}/{file}peak2ucsc_{type}.log"
    conda:  "ucsc.yaml"
    threads: 1
    shell:  "zcat {input.fw} > {output.tfw} 2>> {log} && bedGraphToBigWig {output.tfw} {input.sizes} {output.fw} 2>> {log} && zcat {input.re} > {output.tre} 2>> {log} && bedGraphToBigWig {output.tre} {input.sizes} {output.re} 2>> {log}"

rule GenerateTrack:
    input:  fw = rules.PeakToTRACKS.output.fw,
            re = rules.PeakToTRACKS.output.re
    output: "TRACKS/PEAKS/{combo}/{file}_peak_{type}.fw.bw.trackdone",
            "TRACKS/PEAKS/{combo}/{file}_peak_{type}.re.bw.trackdone"
    log:    "LOGS/PEAKS/{combo}/{file}generatetrack_{type}_peak.log"
    conda:  "base.yaml"
    threads: MAXTHREAD
    params: bwdir = lambda wildcards: "TRACKS/PEAKS/{combo}/{src}".format(combo=combo, src=SETS),
            bins = os.path.abspath(BINS),
            gen = REFDIR,#lambda wildcards: os.path.basename(genomepath(wildcards.file, config)),
            options = '-n Peaks_'+str(PEAKENV)+' -s peaks -l TRACKS_peaks_'+str(PEAKENV)+' -b TRACKS_'+str(PEAKENV),
            uid = lambda wildcards: "{src}".format(src='TRACKS'+os.sep+'PEAKS'+os.sep+combo+SETS.replace(os.sep, '_'))
    shell: "echo -e \"{input.fw}\\n{input.re}\"|python3 {params.bins}/Analysis/GenerateTrackDb.py -i {params.uid} -e 1 -f STDIN -u '' -g {params.gen} {params.options} && touch {input.fw}\.trackdone && touch {input.re}.trackdone 2> {log}"
