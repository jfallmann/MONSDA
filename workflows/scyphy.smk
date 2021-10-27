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


rule UnzipGenome:
    input:  ref = REFERENCE,
    output: fa = expand("{ref}.fa", ref=REFERENCE.replace('.fa.gz', '')),
            fai = expand("{ref}.fa.fai", ref=REFERENCE.replace('.fa.gz', '')),
            fas = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
    log:    expand("LOGS/PEAKS/{combo}/indexfa.log", combo=combo)
    conda:  "samtools.yaml"
    threads: 1
    params: bins = BINS
    shell:  "set +o pipefail; zcat {input[0]} |perl -F\\\\040 -wane 'if($_ =~ /^>/){{$F[0] = $F[0] =~ /^>chr/ ? $F[0] : \">chr\".substr($F[0],1);chomp($F[0]);print \"\\n\".$F[0].\"\\n\"}} else{{($line=$_)=~s/\\r[\\n]*/\\n/gm; chomp($line=$_); print $line}}' |tail -n+2 > {output.fa} && {params.bins}/Preprocessing/indexfa.sh {output.fa} 2> {log} && cut -f1,2 {output.fai} > {output.fas}"

rule UnzipGenome_no_us:
    input:  ref = REFERENCE,
    output: fa = expand("{ref}_us.fa", ref=REFERENCE.replace('.fa.gz', '')),
            fai = expand("{ref}_us.fa.fai", ref=REFERENCE.replace('.fa.gz', '')),
            fas = expand("{ref}_us.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
    log:    expand("LOGS/PEAKS/{combo}/indexfa.log", combo=combo)
    conda:  "samtools.yaml"
    threads: 1
    params: bins = BINS
    shell:  "set +o pipefail; zcat {input[0]} |perl -F\\\\040 -wane 'if($_ =~ /^>/){{$F[0] = $F[0] =~ /^>chr/ ? $F[0] : \">chr\".substr($F[0],1))=~ s/\_/\./g;chomp($F[0]);print \"\\n\".$F[0].\"\\n\"}} else{{($line=$_)=~s/\\r[\\n]*/\\n/gm; chomp($line=$_); print $line}}' |tail -n+2 > {output.fa} && {params.bins}/Preprocessing/indexfa.sh {output.fa} 2> {log} && cut -f1,2 {output.fai} > {output.fas}"

rule remove_softclip:
    input:  bam = "MAPPED/{scombo}/{file}_mapped_{type}.bam",
            fa = REFERENCE,
            refi = expand("{ref}.fai", ref=REFERENCE.replace('.fa.gz', '')),
    output: bam = "MAPPED/{scombo}/{file}_mapped_{type}_nosoftclip.bam",
            bai = "MAPPED/{scombo}/{file}_mapped_{type}_nosoftclip.bam.bai"
    log:    "LOGS/PEAKS/{scombo}/{file}_removesoftclip_{type}.log"
    conda:  "scyphy.yaml"
    threads: 1
    params: bins = BINS,
            tmpidx = lambda x: tempfile.mkdtemp(dir='TMP')
    shell: "python {params.bins}/Analysis/RemoveSoftClip.py -f {input.fa} -b {input.bam} -c -o \'-\' | samtools sort -T {params.tmpidx}/SAMSORT -o {output.bam} --threads {threads} \'-\' 2>> {log} && samtools index {output.bam} 2>> {log} && rm -rf TMP/{params.tmpidx}"


if not all(checklist):
    if not stranded or stranded == 'fr':
        rule BamToBed:
            input:  expand("MAPPED/{scombo}/{file}_mapped_{type}_nosoftclip.bam", scombo=scombo)
            output: "BED/{combo}/{file}_mapped_{type}.bed.gz"
            log:    "LOGS/PEAKS/{combo}/{file}bam2bed_{type}.log"
            conda:  "bedtools.yaml"
            threads: 1
            shell:  "bedtools bamtobed -split -i {input[0]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/2$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |sort -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output[0]} 2> {log}"


    elif stranded and stranded == 'rf':
        rule BamToBed:
            input:  expand("MAPPED/{scombo}/{file}_mapped_{type}_nosoftclip.bam", scombo=scombo))
            output: "BED/{combo}/{file}_mapped_{type}.bed.gz"
            log:    "LOGS/PEAKS/{combo}/{file}bam2bed_{type}.log"
            conda:  "bedtools.yaml"
            threads: 1
            shell:  "bedtools bamtobed -split -i {input[0]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/1$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |sort -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output[0]} 2> {log}"


rule extendbed:
    input:  pks = "BED/{combo}/{file}_mapped_{type}.bed.gz",
            ref = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
    output: ext = "BED/{combo}/{file}_mapped_extended_{type}.bed.gz"
    log:    "LOGS/PEAKS/{combo}/{file}extendbed_{type}.log"
    conda:  "perl.yaml"
    threads: 1
    params: bins = BINS
    shell:  "{params.bins}/Universal/ExtendBed.pl -u 0 -b {input.pks} -o {output.ext} -g {input.ref} 2> {log}"

rule rev_extendbed:
    input:  pks = "BED/{combo}/{file}_mapped_{type}.bed.gz",
            ref = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
    output: ext = "BED/{combo}/{file}_mapped_revtrimmed_{type}.bed.gz"
    log:    "LOGS/PEAKS/{combo}/{file}extendbed_{type}.log"
    conda:  "perl.yaml"
    threads: 1
    params: bins = BINS
    shell:  "{params.bins}/Universal/ExtendBed.pl -d 0 -b {input.pks} -o {output.ext} -g {input.ref}  2> {log}"

if IP == 'iCLIP':
     rule BedToBedg:
        input:  bed = "BED/{combo}/{{file}}_mapped_extended_{{type}}.bed.gz",
                fai = expand("{ref}.fa.fai", ref=REFERENCE.replace('.fa.gz', '')),
                sizes = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
        output: concat = "PEAKS/{combo}/{file}_mapped_{type}.bedg.gz"
        log:    "LOGS/PEAKS/{combo}/{file}bed2bedgraph_{type}.log"
        conda:  "bedtools.yaml"
        threads: 1
        params: bins = BINS,
                odir = lambda wildcards, output:(os.path.dirname(output[0]))
        shell: "export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i {input.bed} -bg -split -strand + -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"+\")'| sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat} 2> {log} && bedtools genomecov -i {input.bed} -bg -split -strand - -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"-\")'|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip >> {output.concat} 2>> {log}"

elif IP == 'revCLIP':
    rule BedToBedg:
        input:  bed = "BED/{combo}/{{file}}_mapped_revtrimmed_{{type}}.bed.gz",
                fai = expand("{ref}.fa.fai", ref=REFERENCE.replace('.fa.gz', '')),
                sizes = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
        output: concat = "PEAKS/{combo}/{file}_mapped_{type}.bedg.gz"
        log:    "LOGS/PEAKS/{combo}/bed2bedgraph_{type}_{file}.log"
        conda:  "bedtools.yaml"
        threads: 1
        params: bins = BINS,
                odir = lambda wildcards, output:(os.path.dirname(output[0]))
        shell: "export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i {input.bed} -bg -split -strand + -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"+\")'| sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat} 2> {log} && bedtools genomecov -i {input.bed} -bg -split -strand - -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"-\")'|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip >> {output.concat} 2>> {log}"

else:
    rule BedToBedg:
        input:  bed = "BED/{combo}/{{file}}_mapped_{{type}}.bed.gz",
                fai = expand("{ref}.fa.fai", ref=REFERENCE.replace('.fa.gz', '')),
                sizes = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
        output: concat = "PEAKS/{combo}/{file}_mapped_{type}.bedg.gz"
        log:    "LOGS/PEAKS/{combo}/bed2bedgraph_{type}_{file}.log"
        conda:  "bedtools.yaml"
        threads: 1
        params: bins = BINS,
                odir = lambda wildcards, output:(os.path.dirname(output[0]))
        shell: "export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i {input.bed} -bg -split -strand + -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"+\")'| sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat} 2> {log} && bedtools genomecov -i {input.bed} -bg -split -strand - -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"-\")'|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip >> {output.concat} 2>> {log}"

rule PreprocessPeaks:
    input:  bedg = rules.BedToBedg.output.concat
    output: pre = "PEAKS/{combo}/{file}_prepeak_{type}.bed.gz",
    log:    "LOGS/PEAKS/{combo}/prepeak_{type}_{file}.log"
    conda:  "perl.yaml"
    threads: 1
    params:  bins = BINS,
             opts = lambda wildcards: tool_params(wildcards.file, None, config, "PEAKS", PEAKENV)['OPTIONS'].get('PREPROCESS', ""),
    shell:  "perl {params.bins}/Analysis/PreprocessPeaks.pl -p <(zcat {input.bedg}) {params.opts} |sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k3,3n -k2,2n -k6,6 | gzip > {output.pre} 2> {log}"

rule FindPeaks:
    input:  pre = rules.PreprocessPeaks.output.pre
    output: peak = "PEAKS/{combo}/{file}_peak_{type}.bed.gz"
    log:    "LOGS/PEAKS/{combo}/findpeaks_{type}_{file}.log"
    conda:  ""+PEAKENV+".yaml"
    threads: 1
    params: ppara = lambda wildcards: tool_params(wildcards.file, None, config, "PEAKS", PEAKENV)['OPTIONS'].get('FINDPEAKS', ""),
            peak = PEAKBIN
    shell:  "set +o pipefail; export LC_ALL=C; if [[ -n \"$(zcat {input.pre} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then {params.peak} {params.ppara} <(zcat {input.pre}) 2> {log}|tail -n+2| sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |grep -v 'nan'| gzip > {output.peak} 2>> {log}; else gzip < /dev/null > {output.peak}; echo \"File {input.pre} empty\" >> {log}; fi"

rule AddSequenceToPeak:
    input:  pk = rules.FindPeaks.output.peak,
            fa = rules.UnzipGenome.output.fa
    output: peak = "PEAKS/{combo}/{file}_peak_seq_{type}.bed.gz",
            pt = temp("PEAKS/{combo}/{file}_peak_chr_{type}.tmp"),
            ps = temp("PEAKS/{combo}/{file}_peak_seq_{type}.tmp")
    log:    "LOGS/PEAKS/{combo}/seq2peaks_{type}_{file}.log"
    conda:  "bedtools.yaml"
    threads: 1
    params: bins=BINS
    shell:  "set +o pipefail; export LC_ALL=C; if [[ -n \"$(zcat {input.pk} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then export LC_ALL=C; zcat {input.pk} | perl -wlane '$F[0] = $F[0] =~ /^chr/ ? $F[0] : \"chr\".$F[0]; print join(\"\\t\",@F[0..5])' > {output.pt} && bedtools getfasta -fi {input.fa} -bed {output.pt} -name -tab -s -fullHeader -fo {output.ps} && cut -d$'\t' -f2 {output.ps}|sed 's/t/u/ig'|paste -d$'\t' <(zcat {input.pk}) - |sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip  > {output.peak} 2> {log}; else gzip < /dev/null > {output.peak} && touch {output.pt} {output.ps}; fi"  # NEED TO GET RID OF SPACES AND WHATEVER IN HEADER

if ANNOPEAK is not None:
    rule AnnotatePeak:
        input:  "PEAKS/{combo}/{file}_peak_seq_{type}.bed.gz"
        output: "PEAKS/{combo}/{file}_peak_anno_{type}.bed.gz"
        log:    "LOGS/PEAKS/{combo}/annotatepeaks_{type}_{file}.log"
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
                trw = temp("PEAKS/{combo}/{file}_peak_{type}.re.tmp.gz"),
        log:    "LOGS/PEAKS/{combo}/{file}_peak2bedg_{type}.log"
        conda:  "perl.yaml"
        threads: 1
        params: bins=BINS
        shell:  "perl {params.bins}/Universal/Bed2Bedgraph.pl -f {input.pk} -c {input.sizes} -p score -x {output.tfw} -y {output.tre} -a track 2>> {log} && zcat {output.tfw}|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n  |gzip > {output.fw} 2>> {log} &&  zcat {output.tre}|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.re} 2>> {log}"

else:
    rule PeakToBedg:
        input:  pk = "PEAKS/{combo}/{file}_peak_{type}.bed.gz",
                sizes = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
        output: fw = "PEAKS/{combo}/{file}_peak_{type}.fw.bedg.gz",
                re = "PEAKS/{combo}/{file}_peak_{type}.re.bedg.gz",
                tfw = temp("PEAKS/{combo}/{file}_peak_{type}.fw.tmp.gz"),
                tre = temp("PEAKS/{combo}/{file}_peak_{type}.re.tmp.gz"),
        log:    "LOGS/PEAKS/{combo}/{file}_peak2bedg_{type}.log"
        conda:  "perl.yaml"
        threads: 1
        params: bins=BINS
        shell:  "perl {params.bins}/Universal/Bed2Bedgraph.pl -f {input.pk} -c {input.sizes} -p score -x {output.tfw} -y {output.tre} -a track 2>> {log} && zcat {output.tfw}|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n  |gzip > {output.fw} 2>> {log} &&  zcat {output.tre}|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.re} 2>> {log}"

### This step normalized the bedg files for comparison in the browser
rule NormalizeBedg:
    input:  fw = rules.PeakToBedg.output.fw,
            re = rules.PeakToBedg.output.re
    output: fw = "PEAKS/{combo}/{file}_peak_{type}.fw.norm.bedg.gz",
            re = "PEAKS/{combo}/{file}_peak_{type}.re.norm.bedg.gz"
    log:    "LOGS/PEAKS/{combo}/ucscpeaknormalizebedgraph_{type}_{file}.log"
    conda:  "perl.yaml"
    threads: 1
    shell: "set +o pipefail; export LC_ALL=C; if [[ -n \"$(zcat {input.fw} | head -c 1 | tr \'\\0\\n\' __)\" ]]; then scale=$(bc <<< \"scale=6;$(zcat {input.fw}|cut -f4|perl -wne '{{$x+=$_;}}END{{if ($x == 0){{$x=1}} print $x}}')/1000000\") perl -wlane '$sc=$ENV{{scale}};print join(\"\t\",@F[0..$#F-1]),\"\t\",$F[-1]/$sc' <(zcat {input.fw})| sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.fw} 2> {log}; else gzip < /dev/null > {output.fw}; echo \"File {input.fw} empty\" >> {log}; fi && if [[ -n \"$(zcat {input.re} | head -c 1 | tr \'\\0\\n\' __)\" ]]; then scale=$(bc <<< \"scale=6;$(zcat {input.re}|cut -f4|perl -wne '{{$x+=$_;}}END{{if ($x == 0){{$x=1}} print $x}}')/1000000\") perl -wlane '$sc=$ENV{{scale}};print join(\"\t\",@F[0..$#F-1]),\"\t\",$F[-1]/$sc' <(zcat {input.re})| sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.re} 2> {log}; else gzip < /dev/null > {output.re}; echo \"File {input.re} empty\" >> {log}; fi"


### This step generates bigwig files for peaks which can then be copied to a web-browsable directory and uploaded to UCSC via the track field
rule PeakToTRACKS:
    input:  fw = rules.NormalizeBedg.output.fw,
            re = rules.NormalizeBedg.output.re,
            fas = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
    output: fw = "TRACKS/PEAKS/{combo}/{file}_peak_{type}.fw.bw",
            re = "TRACKS/PEAKS/{combo}/{file}_peak_{type}.re.bw",
            tfw = temp("TRACKS/PEAKS/{combo}/{file}_{type}fw_tmp"),
            tre = temp("TRACKS/PEAKS/{combo}/{file}_{type}re_tmp")
    log:    "LOGS/PEAKS/{combo}/{file}_peak2ucsc_{type}.log"
    conda:  "ucsc.yaml"
    threads: 1
    shell:  "zcat {input.fw} > {output.tfw} 2>> {log} && bedGraphToBigWig {output.tfw} {input.fas} {output.fw} 2>> {log} && zcat {input.re} > {output.tre} 2>> {log} && bedGraphToBigWig {output.tre} {input.fas} {output.re} 2>> {log}"


rule GenerateTrack:
    input:  fw = rules.PeakToTRACKS.output.fw,
            re = rules.PeakToTRACKS.output.re
    output: "TRACKS/PEAKS/{combo}/{file}_peak_{type}.fw.bw.trackdone",
            "TRACKS/PEAKS/{combo}/{file}_peak_{type}.re.bw.trackdone"
    log:    "LOGS/PEAKS/{combo}/generatetrack_{type}_{file}.log"
    conda:  "base.yaml"
    threads: MAXTHREAD
    params: bwdir = lambda wildcards: "TRACKS/PEAKS/{combo}/{src}".format(combo=combo, src=SETS),
            bins = os.path.abspath(BINS),
            gen = REFDIR,#lambda wildcards: os.path.basename(genomepath(wildcards.file, config)),
            options = '-n Peaks_'+str(PEAKENV)+' -s peaks -l TRACKS_peaks_'+str(PEAKENV)+' -b TRACKS_'+str(PEAKENV),
            uid = lambda wildcards: "{src}".format(src='TRACKS'+os.sep+'PEAKS'+os.sep+combo+SETS.replace(os.sep, '_'))
    shell: "echo -e \"{input.fw}\\n{input.re}\"|python3 {params.bins}/Analysis/GenerateTrackDb.py -i {params.uid} -e 1 -f STDIN -u '' -g {params.gen} {params.options} && touch {input.fw}.trackdone && touch {input.re}.trackdone 2> {log}"