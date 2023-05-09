PEAKBIN, PEAKENV = env_bin_from_config(config, 'PEAKS')

if ANNOPEAK is not None:
    if not rundedup:
        rule themall:
            input:  expand("TRACKS/PEAKS/{combo}/{file}_peak_{type}.fw.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("TRACKS/PEAKS/{combo}/{file}_peak_{type}.re.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("TRACKS/PEAKS/{combo}/{file}_mapped_{type}.fw.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("TRACKS/PEAKS/{combo}/{file}_mapped_{type}.re.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.fw.bedg.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.re.bedg.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.fw.norm.bedg.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.re.norm.bedg.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_seq_{type}.bed.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_anno_{type}.bed.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique'])
    else:
        rule themall:
            input:  expand("TRACKS/PEAKS/{combo}/{file}_peak_{type}.fw.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("TRACKS/PEAKS/{combo}/{file}_peak_{type}.re.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("TRACKS/PEAKS/{combo}/{file}_mapped_{type}.fw.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("TRACKS/PEAKS/{combo}/{file}_mapped_{type}.re.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
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
                    expand("TRACKS/PEAKS/{combo}/{file}_mapped_{type}.fw.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("TRACKS/PEAKS/{combo}/{file}_mapped_{type}.re.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.fw.bedg.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.re.bedg.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.bed.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_seq_{type}.bed.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique'])
    else:
        rule themall:
            input:  expand("TRACKS/PEAKS/{combo}/{file}_peak_{type}.fw.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("TRACKS/PEAKS/{combo}/{file}_peak_{type}.re.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("TRACKS/PEAKS/{combo}/{file}_mapped_{type}.fw.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("TRACKS/PEAKS/{combo}/{file}_mapped_{type}.re.bw.trackdone", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.fw.bedg.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.re.bedg.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.bed.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}/{file}_peak_seq_{type}.bed.gz", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup'])

checklist = list()
for file in samplecond(SAMPLES, config):
    checktype = ['sorted', 'unique'] if not rundedup else ['sorted', 'unique', 'sorted_dedup', 'sorted_unique_dedup']
    for type in checktype:
        checklist.append(os.path.isfile(os.path.abspath('BED/'+scombo+file+'_mapped_'+type+'_nosoftclip.bed.gz')) and not os.path.islink(os.path.abspath('BED/'+scombo+file+'_mapped_'+type+'_nosoftclip.bed.gz')) and os.path.isfile(os.path.abspath('MAPPED/'+scombo+file+'_mapped_'+type+'_nosoftclip.bam')) and not os.path.islink(os.path.abspath('MAPPED/'+scombo+file+'_mapped_'+type+'_nosoftclip.bam')))


include: "manipulate_genome.smk"

if not all(checklist):
    rule remove_softclip:
        input:  bam = "MAPPED/{scombo}/{file}_mapped_{type}.bam",
                fa = REFERENCE,
                refi = expand("{ref}.fa.fai", ref=REFERENCE.replace('.fa.gz', '')),
        output: bam = "MAPPED/{scombo}/{file}_mapped_{type}_nosoftclip.bam",
                bai = "MAPPED/{scombo}/{file}_mapped_{type}_nosoftclip.bam.bai"
        log:    "LOGS/PEAKS/{scombo}/{file}_removesoftclip_{type}.log"
        conda:  "scyphy.yaml"
        threads: 1
        params: bins = BINS,
                tmpidx = lambda x: tempfile.mkdtemp(dir='TMP'),
                opts = lambda wildcards: tool_params(wildcards.file, None, config, "PEAKS", PEAKENV)['OPTIONS'].get('SOFTCLIP', "")
        shell: "python {params.bins}/Analysis/RemoveSoftClip.py -f {input.fa} -b {input.bam} {params.opts} -o \'-\' | samtools sort -T {params.tmpidx}/SAMSORT -o {output.bam} --threads {threads} \'-\' 2>> {log} && samtools index {output.bam} 2>> {log} && rm -rf TMP/{params.tmpidx}"

    if not stranded or (stranded == 'fr' or stranded == 'ISF'):
        rule BamToBed:
            input:  "MAPPED/{scombo}/{file}_mapped_{type}_nosoftclip.bam",
            output: "BED/{scombo}/{file}_mapped_{type}_nosoftclip.bed.gz"
            log:    "LOGS/PEAKS/{scombo}/{file}bam2bed_{type}.log"
            conda:  "bedtools.yaml"
            threads: 1
            params: sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
            shell:  "bedtools bamtobed -split -i {input[0]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/2$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |sort -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output[0]} 2> {log}"

    elif stranded and (stranded == 'rf' or stranded == 'ISR'):
        rule BamToBed:
            input:  "MAPPED/{scombo}/{file}_mapped_{type}_nosoftclip.bam",
            output: "BED/{scombo}/{file}_mapped_{type}_nosoftclip.bed.gz"
            log:    "LOGS/PEAKS/{scombo}/{file}bam2bed_{type}.log"
            conda:  "bedtools.yaml"
            threads: 1
            params: sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
            shell:  "bedtools bamtobed -split -i {input[0]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/1$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |sort -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output[0]} 2> {log}"


rule extendbed:
    input:  pks = "BED/{scombo}/{file}_mapped_{type}_nosoftclip.bed.gz",
            ref = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
    output: ext = "BED/{scombo}/{file}_mapped_extended_{type}_nosoftclip.bed.gz"
    log:    "LOGS/PEAKS/{scombo}/{file}extendbed_{type}.log"
    conda:  "perl.yaml"
    threads: 1
    params: bins = BINS
    shell:  "{params.bins}/Universal/ExtendBed.pl -u 0 -b {input.pks} -o {output.ext} -g {input.ref} 2> {log}"

rule rev_extendbed:
    input:  pks = "BED/{scombo}/{file}_mapped_{type}_nosoftclip.bed.gz",
            ref = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
    output: ext = "BED/{scombo}/{file}_mapped_revtrimmed_{type}_nosoftclip.bed.gz"
    log:    "LOGS/PEAKS/{scombo}/{file}extendbed_{type}.log"
    conda:  "perl.yaml"
    threads: 1
    params: bins = BINS
    shell:  "{params.bins}/Universal/ExtendBed.pl -d 0 -b {input.pks} -o {output.ext} -g {input.ref}  2> {log}"

rule BedToBedg_ext:
    input:  bed = expand("BED/{scombo}/{{file}}_mapped_extended_{{type}}_nosoftclip.bed.gz", scombo=scombo),
            fai = expand("{ref}.fa.fai", ref=REFERENCE.replace('.fa.gz', '')),
            sizes = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
    output: concat_fw = "PEAKS/{combo}/{file}_mapped_{type}_nosoftclip_ext.fw.bedg.gz",
            concat_re = "PEAKS/{combo}/{file}_mapped_{type}_nosoftclip_ext.re.bedg.gz",
            concat = "PEAKS/{combo}/{file}_mapped_{type}_nosoftclip_ext.bedg.gz",
            tosrt = temp("PEAKS/{combo}/{file}_mapped_{type}_ext.unsrt")
    log:    "LOGS/PEAKS/{combo}/{file}bed2bedgraph_ext_{type}.log"
    conda:  "bedtools.yaml"
    threads: 1
    params: bins = BINS,
            odir = lambda wildcards, output:(os.path.dirname(output[0])),
            sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
    shell: "export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i {input.bed} -bg -split -strand + -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"+\")' > {output.tosrt} 2> {log} && cat {output.tosrt}| sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat_fw} 2>> {log} && bedtools genomecov -i {input.bed} -bg -split -strand - -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"-\")' > {output.tosrt} && cat {output.tosrt}| sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat_re} 2>> {log} && zcat {output.concat_fw} {output.concat_re} | sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat}"

rule BedToBedg_rev:
    input:  bed = expand("BED/{scombo}/{{file}}_mapped_revtrimmed_{{type}}_nosoftclip.bed.gz", scombo=scombo),
            fai = expand("{ref}.fa.fai", ref=REFERENCE.replace('.fa.gz', '')),
            sizes = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
    output: concat_fw = "PEAKS/{combo}/{file}_mapped_{type}_nosoftclip_rev.fw.bedg.gz",
            concat_re = "PEAKS/{combo}/{file}_mapped_{type}_nosoftclip_rev.re.bedg.gz",
            concat = "PEAKS/{combo}/{file}_mapped_{type}_nosoftclip_rev.bedg.gz",
            tosrt = temp("PEAKS/{combo}/{file}_mapped_{type}_rev.unsrt")
    log:    "LOGS/PEAKS/{combo}/bed2bedgraph_rev_{type}_{file}.log"
    conda:  "bedtools.yaml"
    threads: 1
    params: bins = BINS,
            odir = lambda wildcards, output:(os.path.dirname(output[0])),
            sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
    shell: "export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i {input.bed} -bg -split -strand + -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"+\")' > {output.tosrt} 2> {log} && cat {output.tosrt}| sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat_fw} 2>> {log} && bedtools genomecov -i {input.bed} -bg -split -strand - -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"-\")' > {output.tosrt} && cat {output.tosrt}| sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat_re} 2>> {log} && zcat {output.concat_fw} {output.concat_re} | sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat}"

rule BedToBedg:
    input:  bed = expand("BED/{scombo}/{{file}}_mapped_{{type}}_nosoftclip.bed.gz", scombo=scombo),
            fai = expand("{ref}.fa.fai", ref=REFERENCE.replace('.fa.gz', '')),
            sizes = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
    output: concat_fw = "PEAKS/{combo}/{file}_mapped_{type}_nosoftclip.fw.bedg.gz",
            concat_re = "PEAKS/{combo}/{file}_mapped_{type}_nosoftclip.re.bedg.gz",
            concat = "PEAKS/{combo}/{file}_mapped_{type}_nosoftclip.bedg.gz",
            tosrt = temp("PEAKS/{combo}/{file}_mapped_{type}.unsrt")
    log:    "LOGS/PEAKS/{combo}/bed2bedgraph_{type}_{file}.log"
    conda:  "bedtools.yaml"
    threads: 1
    params: bins = BINS,
            odir = lambda wildcards, output:(os.path.dirname(output[0])),
            sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
    shell: "export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i {input.bed} -bg -split -strand + -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"+\")' > {output.tosrt} 2> {log} && cat {output.tosrt}| sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat_fw} 2>> {log} && bedtools genomecov -i {input.bed} -bg -split -strand - -g {input.sizes} |perl -wlane 'print join(\"\t\",@F[0..2],\".\",$F[3],\"-\")' > {output.tosrt} && cat {output.tosrt}| sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat_re} 2>> {log} && zcat {output.concat_fw} {output.concat_re} | sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.concat}"

if any(x == IP for x in ['iCLIP', 'CLIP']):
    rule PreprocessPeaks:
        input:  bedg = rules.BedToBedg_ext.output.concat
        output: pre = "PEAKS/{combo}/{file}_prepeak_{type}_nosoftclip.bed.gz"
        log:    "LOGS/PEAKS/{combo}/prepeak_{type}_{file}.log"
        conda:  "perl.yaml"
        threads: 1
        params: bins = BINS,
                opts = lambda wildcards: tool_params(wildcards.file, None, config, "PEAKS", PEAKENV)['OPTIONS'].get('PREPROCESS', "")
        shell:  "perl {params.bins}/Analysis/PreprocessPeaks.pl -p <(zcat {input.bedg}) {params.opts} | sort -t$'\t' -k1,1 -k3,3n -k2,2n -k6,6 | gzip > {output.pre} 2> {log}"

elif IP == 'revCLIP':
    rule PreprocessPeaks:
        input:  bedg = rules.BedToBedg_rev.output.concat
        output: pre = "PEAKS/{combo}/{file}_prepeak_{type}_nosoftclip.bed.gz"
        log:    "LOGS/PEAKS/{combo}/prepeak_{type}_{file}.log"
        conda:  "perl.yaml"
        threads: 1
        params:  bins = BINS,
                opts = lambda wildcards: tool_params(wildcards.file, None, config, "PEAKS", PEAKENV)['OPTIONS'].get('PREPROCESS', "")
        shell:  "perl {params.bins}/Analysis/PreprocessPeaks.pl -p <(zcat {input.bedg}) {params.opts} | sort -t$'\t' -k1,1 -k3,3n -k2,2n -k6,6 | gzip > {output.pre} 2> {log}"

rule FindPeaks:
    input:  pre = rules.PreprocessPeaks.output.pre
    output: peak = "PEAKS/{combo}/{file}_peak_{type}.bed.gz"
    log:    "LOGS/PEAKS/{combo}/findpeaks_{type}_{file}.log"
    conda:  ""+PEAKENV+".yaml"
    threads: 1
    params: ppara = lambda wildcards: tool_params(wildcards.file, None, config, "PEAKS", PEAKENV)['OPTIONS'].get('FINDPEAKS', ""),
            peak = PEAKBIN,
            sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
    shell:  "set +o pipefail; export LC_ALL=C; if [[ -n \"$(zcat {input.pre} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then {params.peak} {params.ppara} <(zcat {input.pre}|sort -t$'\t' -k1,1 -k3,3n -k2,2n -k6,6) 2> {log}|tail -n+2| sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |grep -v 'nan'| gzip > {output.peak} 2>> {log}; else gzip < /dev/null > {output.peak}; echo \"File {input.pre} empty\" >> {log}; fi"

rule AddSequenceToPeak:
    input:  pk = rules.FindPeaks.output.peak,
            fa = rules.UnzipGenome.output.fa
    output: peak = "PEAKS/{combo}/{file}_peak_seq_{type}.bed.gz",
            pt = temp("PEAKS/{combo}/{file}_peak_chr_{type}.tmp"),
            ps = temp("PEAKS/{combo}/{file}_peak_seq_{type}.tmp")
    log:    "LOGS/PEAKS/{combo}/seq2peaks_{type}_{file}.log"
    conda:  "bedtools.yaml"
    threads: 1
    params: bins=BINS,
            sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
    shell:  "set +o pipefail; export LC_ALL=C; if [[ -n \"$(zcat {input.pk} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then zcat {input.pk} | perl -wlane '$F[0] = $F[0] =~ /^chr/ ? $F[0] : \"chr\".$F[0]; print join(\"\\t\",@F[0..5])' > {output.pt} && bedtools getfasta -fi {input.fa} -bed {output.pt} -name -tab -s -fullHeader -fo {output.ps} && cut -d$'\t' -f2 {output.ps}|sed 's/t/u/ig'|paste -d$'\t' <(zcat {input.pk}) - |sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip  > {output.peak} 2> {log}; else gzip < /dev/null > {output.peak} && touch {output.pt} {output.ps}; fi" 

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
        params: bins=BINS,
                sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
        shell:  "perl {params.bins}/Universal/Bed2Bedgraph.pl -f {input.pk} -c {input.sizes} -p score -x {output.tfw} -y {output.tre} -a track 2>> {log} && zcat {output.tfw}|sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n  |gzip > {output.fw} 2>> {log} &&  zcat {output.tre}|sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.re} 2>> {log}"

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
        params: bins=BINS,
                sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
        shell:  "perl {params.bins}/Universal/Bed2Bedgraph.pl -f {input.pk} -c {input.sizes} -p score -x {output.tfw} -y {output.tre} -a track 2>> {log} && zcat {output.tfw}|sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n  |gzip > {output.fw} 2>> {log} &&  zcat {output.tre}|sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.re} 2>> {log}"

### This step normalized the bedg files for comparison in the browser
rule NormalizeBedg:
    input:  fw = rules.PeakToBedg.output.fw,
            re = rules.PeakToBedg.output.re,
            map_fw = rules.BedToBedg.output.concat_fw,
            map_re = rules.BedToBedg.output.concat_re,
            map_fw_ext = rules.BedToBedg_ext.output.concat_fw,
            map_re_ext = rules.BedToBedg_ext.output.concat_re,
            map_fw_rev = rules.BedToBedg_rev.output.concat_fw,
            map_re_rev = rules.BedToBedg_rev.output.concat_re
    output: fw = "PEAKS/{combo}/{file}_peak_{type}.fw.norm.bedg.gz",
            re = "PEAKS/{combo}/{file}_peak_{type}.re.norm.bedg.gz",
            map_fw = "PEAKS/{combo}/{file}_mapped_{type}.fw.norm.bedg.gz",
            map_re = "PEAKS/{combo}/{file}_mapped_{type}.re.norm.bedg.gz",
            map_fw_ext = "PEAKS/{combo}/{file}_mapped_{type}_ext.fw.norm.bedg.gz",
            map_re_ext = "PEAKS/{combo}/{file}_mapped_{type}_ext.re.norm.bedg.gz",
            map_fw_rev = "PEAKS/{combo}/{file}_mapped_{type}_rev.fw.norm.bedg.gz",
            map_re_rev = "PEAKS/{combo}/{file}_mapped_{type}_rev.re.norm.bedg.gz"
    log:    "LOGS/PEAKS/{combo}/ucscpeaknormalizebedgraph_{type}_{file}.log"
    conda:  "perl.yaml"
    threads: 1
    params: sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
    shell: "set +o pipefail; export LC_ALL=C; if [[ -n \"$(zcat {input.fw} | head -c 1 | tr \'\\0\\n\' __)\" ]]; then scale=$(bc <<< \"scale=6;$(zcat {input.fw}|cut -f4|perl -wne '{{$x+=$_;}}END{{if ($x == 0){{$x=1}} print $x}}')/1000000\") perl -wlane '$sc=$ENV{{scale}};print join(\"\t\",@F[0..$#F-1]),\"\t\",$F[-1]/$sc' <(zcat {input.fw})| sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.fw} 2> {log}; else gzip < /dev/null > {output.fw}; echo \"File {input.fw} empty\" >> {log}; fi && if [[ -n \"$(zcat {input.re} | head -c 1 | tr \'\\0\\n\' __)\" ]]; then scale=$(bc <<< \"scale=6;$(zcat {input.re}|cut -f4|perl -wne '{{$x+=$_;}}END{{if ($x == 0){{$x=1}} print $x}}')/1000000\") perl -wlane '$sc=$ENV{{scale}};print join(\"\t\",@F[0..$#F-1]),\"\t\",$F[-1]/$sc' <(zcat {input.re})| sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.re} 2> {log}; else gzip < /dev/null > {output.re}; echo \"File {input.re} empty\" >> {log}; fi && if [[ -n \"$(zcat {input.map_fw} | head -c 1 | tr \'\\0\\n\' __)\" ]]; then scale=$(bc <<< \"scale=6;$(zcat {input.map_fw}|cut -f5|perl -wne '{{$x+=$_;}}END{{if ($x == 0){{$x=1}} print $x}}')/1000000\") perl -wlane '$sc=$ENV{{scale}};print join(\"\t\",@F[0..2]),\"\t\",$F[4]/$sc' <(zcat {input.map_fw})| sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.map_fw} 2> {log}; else gzip < /dev/null > {output.map_fw}; echo \"File {input.map_fw} empty\" >> {log}; fi && if [[ -n \"$(zcat {input.map_re} | head -c 1 | tr \'\\0\\n\' __)\" ]]; then scale=$(bc <<< \"scale=6;$(zcat {input.map_re}|cut -f5|perl -wne '{{$x+=$_;}}END{{if ($x == 0){{$x=1}} print $x}}')/1000000\") perl -wlane '$sc=$ENV{{scale}};print join(\"\t\",@F[0..2]),\"\t\",$F[4]/$sc' <(zcat {input.map_re})| sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.map_re} 2> {log}; else gzip < /dev/null > {output.map_re}; echo \"File {input.map_re} empty\" >> {log}; fi && if [[ -n \"$(zcat {input.map_fw_ext} | head -c 1 | tr \'\\0\\n\' __)\" ]]; then scale=$(bc <<< \"scale=6;$(zcat {input.map_fw_ext}|cut -f5|perl -wne '{{$x+=$_;}}END{{if ($x == 0){{$x=1}} print $x}}')/1000000\") perl -wlane '$sc=$ENV{{scale}};print join(\"\t\",@F[0..2]),\"\t\",$F[4]/$sc' <(zcat {input.map_fw_ext})| sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.map_fw_ext} 2> {log}; else gzip < /dev/null > {output.map_fw_ext}; echo \"File {input.map_fw_ext} empty\" >> {log}; fi && if [[ -n \"$(zcat {input.map_re_ext} | head -c 1 | tr \'\\0\\n\' __)\" ]]; then scale=$(bc <<< \"scale=6;$(zcat {input.map_re_ext}|cut -f5|perl -wne '{{$x+=$_;}}END{{if ($x == 0){{$x=1}} print $x}}')/1000000\") perl -wlane '$sc=$ENV{{scale}};print join(\"\t\",@F[0..2]),\"\t\",$F[4]/$sc' <(zcat {input.map_re_ext})| sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.map_re_ext} 2> {log}; else gzip < /dev/null > {output.map_re_ext}; echo \"File {input.map_re_ext} empty\" >> {log}; fi && if [[ -n \"$(zcat {input.map_fw_rev} | head -c 1 | tr \'\\0\\n\' __)\" ]]; then scale=$(bc <<< \"scale=6;$(zcat {input.map_fw_rev}|cut -f5|perl -wne '{{$x+=$_;}}END{{if ($x == 0){{$x=1}} print $x}}')/1000000\") perl -wlane '$sc=$ENV{{scale}};print join(\"\t\",@F[0..2]),\"\t\",$F[4]/$sc' <(zcat {input.map_fw_rev})| sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.map_fw_rev} 2> {log}; else gzip < /dev/null > {output.map_fw_rev}; echo \"File {input.map_fw_rev} empty\" >> {log}; fi && if [[ -n \"$(zcat {input.map_re_rev} | head -c 1 | tr \'\\0\\n\' __)\" ]]; then scale=$(bc <<< \"scale=6;$(zcat {input.map_re_rev}|cut -f5|perl -wne '{{$x+=$_;}}END{{if ($x == 0){{$x=1}} print $x}}')/1000000\") perl -wlane '$sc=$ENV{{scale}};print join(\"\t\",@F[0..2]),\"\t\",$F[4]/$sc' <(zcat {input.map_re_rev})| sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.map_re_rev} 2> {log}; else gzip < /dev/null > {output.map_re_rev}; echo \"File {input.map_re_rev} empty\" >> {log}; fi"


### This step generates bigwig files for peaks which can then be copied to a web-browsable directory and uploaded to a genome browser via the track field
rule PeakToTRACKS:
    input:  fw = rules.NormalizeBedg.output.fw,
            re = rules.NormalizeBedg.output.re,
            map_fw = rules.NormalizeBedg.output.map_fw,
            map_re = rules.NormalizeBedg.output.map_re,
            map_fw_ext = rules.NormalizeBedg.output.map_fw_ext,
            map_re_ext = rules.NormalizeBedg.output.map_re_ext,
            map_fw_rev = rules.NormalizeBedg.output.map_fw_rev,
            map_re_rev = rules.NormalizeBedg.output.map_re_rev,
            fas = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
    output: fw = "TRACKS/PEAKS/{combo}/{file}_peak_{type}.fw.bw",
            re = "TRACKS/PEAKS/{combo}/{file}_peak_{type}.re.bw",
            map_fw = "TRACKS/PEAKS/{combo}/{file}_mapped_{type}.fw.bw",
            map_re = "TRACKS/PEAKS/{combo}/{file}_mapped_{type}.re.bw",
            map_fw_ext = "TRACKS/PEAKS/{combo}/{file}_mapped_{type}_ext.fw.bw",
            map_re_ext = "TRACKS/PEAKS/{combo}/{file}_mapped_{type}_ext.re.bw",
            map_fw_rev = "TRACKS/PEAKS/{combo}/{file}_mapped_{type}_rev.fw.bw",
            map_re_rev = "TRACKS/PEAKS/{combo}/{file}_mapped_{type}_rev.re.bw",
            tfw = temp("TRACKS/PEAKS/{combo}/{file}_{type}fw_tmp"),
            tre = temp("TRACKS/PEAKS/{combo}/{file}_{type}re_tmp"),
            tmapfw = temp("TRACKS/PEAKS/{combo}/{file}_{type}mapfw_tmp"),
            tmapre = temp("TRACKS/PEAKS/{combo}/{file}_{type}mapre_tmp"),
            tmapfw_ext = temp("TRACKS/PEAKS/{combo}/{file}_{type}mapfw_ext_tmp"),
            tmapre_ext = temp("TRACKS/PEAKS/{combo}/{file}_{type}mapre_ext_tmp"),
            tmapfw_rev = temp("TRACKS/PEAKS/{combo}/{file}_{type}mapfw_rev_tmp"),
            tmapre_rev = temp("TRACKS/PEAKS/{combo}/{file}_{type}mapre_rev_tmp")
    log:    "LOGS/PEAKS/{combo}/{file}_peak2ucsc_{type}.log"
    conda:  "ucsc.yaml"
    threads: 1
    shell:  "export LC_ALL=C; if [[ -n \"$(zcat {input.fw} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then zcat {input.fw} > {output.tfw} 2>> {log} && bedGraphToBigWig {output.tfw} {input.fas} {output.fw} 2>> {log}; else touch {output.tfw}; gzip < /dev/null > {output.fw}; echo \"File {input.fw} empty\" >> {log}; fi && if [[ -n \"$(zcat {input.re} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then zcat {input.re} > {output.tre} 2>> {log} && bedGraphToBigWig {output.tre} {input.fas} {output.re} 2>> {log}; else touch {output.tre}; gzip < /dev/null > {output.re}; echo \"File {input.re} empty\" >> {log}; fi && zcat {input.map_fw} > {output.tmapfw} 2>> {log} && bedGraphToBigWig {output.tmapfw} {input.fas} {output.map_fw} 2>> {log} && zcat {input.map_re} > {output.tmapre} 2>> {log} && bedGraphToBigWig {output.tmapre} {input.fas} {output.map_re} 2>> {log} && zcat {input.map_fw_ext} > {output.tmapfw_ext} 2>> {log} && bedGraphToBigWig {output.tmapfw_ext} {input.fas} {output.map_fw_ext} 2>> {log} && zcat {input.map_re_ext} > {output.tmapre_ext} 2>> {log} && bedGraphToBigWig {output.tmapre_ext} {input.fas} {output.map_re_ext} 2>> {log} && zcat {input.map_fw_rev} > {output.tmapfw_rev} 2>> {log} && bedGraphToBigWig {output.tmapfw_rev} {input.fas} {output.map_fw_rev} 2>> {log} && zcat {input.map_re_rev} > {output.tmapre_rev} 2>> {log} && bedGraphToBigWig {output.tmapre_rev} {input.fas} {output.map_re_rev} 2>> {log}"

rule GenerateTrack:
    input:  fw = rules.PeakToTRACKS.output.fw,
            re = rules.PeakToTRACKS.output.re,
            map_fw = rules.PeakToTRACKS.output.map_fw,
            map_re = rules.PeakToTRACKS.output.map_re,
            map_fw_ext = rules.PeakToTRACKS.output.map_fw_ext,
            map_re_ext = rules.PeakToTRACKS.output.map_re_ext,
            map_fw_rev = rules.PeakToTRACKS.output.map_fw_rev,
            map_re_rev = rules.PeakToTRACKS.output.map_re_rev
    output: "TRACKS/PEAKS/{combo}/{file}_peak_{type}.fw.bw.trackdone",
            "TRACKS/PEAKS/{combo}/{file}_peak_{type}.re.bw.trackdone",
            "TRACKS/PEAKS/{combo}/{file}_mapped_{type}.fw.bw.trackdone",
            "TRACKS/PEAKS/{combo}/{file}_mapped_{type}.re.bw.trackdone",
            "TRACKS/PEAKS/{combo}/{file}_mapped_{type}_ext.fw.bw.trackdone",
            "TRACKS/PEAKS/{combo}/{file}_mapped_{type}_ext.re.bw.trackdone",
            "TRACKS/PEAKS/{combo}/{file}_mapped_{type}_rev.fw.bw.trackdone",
            "TRACKS/PEAKS/{combo}/{file}_mapped_{type}_rev.re.bw.trackdone"
    log:    "LOGS/PEAKS/{combo}/generatetrack_{type}_{file}.log"
    conda:  "base.yaml"
    threads: MAXTHREAD
    params: bwdir = lambda wildcards: "TRACKS/PEAKS/{combo}/{src}".format(combo=combo, src=SETS),
            bins = os.path.abspath(BINS),
            gen = REFDIR,#lambda wildcards: os.path.basename(genomepath(wildcards.file, config)),
            options = '-n Peaks_'+str(PEAKENV)+' -s peaks -l TRACKS_peaks_'+str(PEAKENV)+' -b TRACKS_'+str(PEAKENV),
            uid = lambda wildcards: "{src}".format(src='TRACKS'+os.sep+'PEAKS'+os.sep+combo+SETS.replace(os.sep, '_'))
    shell: "echo -e \"{input.fw}\\n{input.re}\\n{input.map_fw}\\n{input.map_re}\"|python3 {params.bins}/Analysis/GenerateTrackDb.py -i {params.uid} -e 1 -f STDIN -u '' -g {params.gen} {params.options} && touch {input.fw}.trackdone && touch {input.re}.trackdone && touch {input.map_fw}.trackdone && touch {input.map_re}.trackdone 2> {log} && touch {input.map_fw_ext}.trackdone && touch {input.map_re_ext}.trackdone 2> {log} && touch {input.map_fw_rev}.trackdone && touch {input.map_re_rev}.trackdone 2> {log}"
