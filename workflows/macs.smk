PEAKBIN, PEAKENV = env_bin_from_config2(SAMPLES,config,'PEAKS')

SAMPLES = set_pairings(SAMPLES, config)

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
                    expand("PEAKS/{file}_peak_seq_{type}.bed.gz",file=samplecond(SAMPLES,config), type=['sorted','sorted_unique'])
    else:
        rule themall:
            input:  expand("PEAKS/{file}_mapped_{type}.bedg.gz",file=samplecond(SAMPLES,config), type=['sorted_dedup','sorted_unique_dedup']),
                    expand("UCSC/{file}_peak_{type}.fw.bw.trackdone",file=samplecond(SAMPLES,config), type=['sorted_dedup','sorted_unique_dedup']),
                    expand("UCSC/{file}_peak_{type}.re.bw.trackdone",file=samplecond(SAMPLES,config), type=['sorted_dedup','sorted_unique_dedup']),
                    expand("UCSC/{file}_peak_{type}.fw.bedg.gz",file=samplecond(SAMPLES,config), type=['sorted_dedup','sorted_unique_dedup']),
                    expand("UCSC/{file}_peak_{type}.re.bedg.gz",file=samplecond(SAMPLES,config), type=['sorted_dedup','sorted_unique_dedup']),
                    expand("PEAKS/{file}_peak_{type}.bed.gz",file=samplecond(SAMPLES,config), type=['sorted_dedup','sorted_unique_dedup']),
                    expand("PEAKS/{file}_peak_seq_{type}.bed.gz",file=samplecond(SAMPLES,config), type=['sorted_dedup','sorted_unique_dedup'])


checklist = list()
checklist2 = list()
for file in samplecond(SAMPLES,config):
    for type in ['sorted','piunique']:
        checklist.append(os.path.isfile(os.path.abspath('BED/'+file+'_mapped_'+type+'.bed.gz')) and not os.path.islink(os.path.abspath('BED/'+file+'_mapped_'+type+'.bed.gz')))
        checklist2.append(os.path.isfile(os.path.abspath('UCSC/'+file+'_mapped_'+type+'.bed.gz')) and not os.path.islink(os.path.abspath('UCSC/'+file+'_mapped_'+type+'.bed.gz')))

rule RunMacs:
    input:  bam = "MAPPED/{file}_mapped_{type}.bam"
    output: peak = "PEAKS/{file}_peak_{type}_peaks.xls"
    log:    "LOGS/PEAKS/findpeaks_macs_{type}_{file}.log"
    conda:  "nextsnakes/envs/"+PEAKENV+".yaml"
    threads: 1
    params: pairing = lambda wildcards: get_pairing(wildcards.file, wildcards.type, config),
            ppara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "PEAKS")['OPTIONS'][1].items()),
            peak = PEAKBIN,
            outdir = lambda wildcards, output: os.path.dirname(output.peak),
            outname = lambda wildcards: os.path.filename(wildcards.file),
    shell:  "{params.peak} callpeak -t {input.bam} {params.pairing} --outdir {params.outdir} -n {params.outname} -f BAM {params.ppara} 2> {log}"

rule UnzipGenome:
    input:  ref = REFERENCE,
    output: fa = expand("{ref}_fastafrombed.fa",ref=REFERENCE.replace('.fa.gz',''))
    log:    "LOGS/PEAKS/indexfa.log"
    conda:  "nextsnakes/envs/samtools.yaml"
    threads: 1
    params: bins = BINS
    shell:  "zcat {input[0]} |perl -F\\\\040 -wlane 'if($_ =~ /^>/){{($F[0] = $F[0] =~ /^>chr/ ? $F[0] : \">chr\".substr($F[0],1))=~ s/\_/\./g;print $F[0]}}else{{print}}' > {output.fa} && {params.bins}/Preprocessing/indexfa.sh {output.fa} 2> {log}"

rule AddSequenceToPeak:
    input:  pk = rules.FindPeaks.output.peak,
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
                trw = temp("UCSC/{file}_peak_{type}.re.tmp.gz"),
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
            uid = lambda wildcards: "{src}".format(src='UCSC'+os.sep+'PEAKS_'+SETS.replace(os.sep,'_'))
    shell: "echo -e \"{input.fw}\\n{input.re}\"|python3 {params.bins}/Analysis/GenerateTrackDb.py -i {params.uid} -e 1 -f STDIN -u '' -g {params.gen} {params.options} && touch {input.fw}\.trackdone && touch {input.re}.trackdone 2> {log}"
