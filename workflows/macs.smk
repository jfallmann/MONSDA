PEAKBIN, PEAKENV = env_bin_from_config3(config, 'PEAKS')
PEAKSAMPLES = set_pairing(SAMPLES, config)

log.info('PEAKSAMPLES: '+str(PEAKSAMPLES))

if ANNOPEAK is not None:
    if not rundedup:
        rule themall:
            input:  expand("TRACKS/PEAKS/{combo}/{file}_peak_{type}.fw.bw.trackdone", combo=combo, file=samplecond(PEAKSAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("TRACKS/PEAKS/{combo}/{file}_peak_{type}.re.bw.trackdone", combo=combo, file=samplecond(PEAKSAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.fw.bedg.gz", combo=combo, file=samplecond(PEAKSAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.re.bedg.gz", combo=combo, file=samplecond(PEAKSAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.bed.gz", combo=combo, file=samplecond(PEAKSAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_seq_{type}.bed.gz", combo=combo, file=samplecond(PEAKSAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_anno_{type}.bed.gz", combo=combo, file=samplecond(PEAKSAMPLES, config), type=['sorted', 'sorted_unique'])
    else:
        rule themall:
            input:  expand("TRACKS/PEAKS/{combo}/{file}_peak_{type}.fw.bw.trackdone", combo=combo, file=samplecond(PEAKSAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("TRACKS/PEAKS/{combo}/{file}_peak_{type}.re.bw.trackdone", combo=combo, file=samplecond(PEAKSAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.fw.bedg.gz", combo=combo, file=samplecond(PEAKSAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.re.bedg.gz", combo=combo, file=samplecond(PEAKSAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.bed.gz", combo=combo, file=samplecond(PEAKSAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}/{file}_peak_seq_{type}.bed.gz", combo=combo, file=samplecond(PEAKSAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}/{file}_peak_anno_{type}.bed.gz", combo=combo, file=samplecond(PEAKSAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup'])

else:
    if not rundedup:
        rule themall:
            input:  expand("TRACKS/PEAKS/{combo}/{file}_peak_{type}.fw.bw.trackdone", combo=combo, file=samplecond(PEAKSAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("TRACKS/PEAKS/{combo}/{file}_peak_{type}.re.bw.trackdone", combo=combo, file=samplecond(PEAKSAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.fw.bedg.gz", combo=combo, file=samplecond(PEAKSAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.re.bedg.gz", combo=combo, file=samplecond(PEAKSAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.bed.gz", combo=combo, file=samplecond(PEAKSAMPLES, config), type=['sorted', 'sorted_unique']),
                    expand("PEAKS/{combo}/{file}_peak_seq_{type}.bed.gz", combo=combo, file=samplecond(PEAKSAMPLES, config), type=['sorted', 'sorted_unique'])
    else:
        rule themall:
            input:  expand("TRACKS/PEAKS/{combo}/{file}_peak_{type}.fw.bw.trackdone", combo=combo, file=samplecond(PEAKSAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("TRACKS/PEAKS/{combo}/{file}_peak_{type}.re.bw.trackdone", combo=combo, file=samplecond(PEAKSAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.fw.bedg.gz", combo=combo, file=samplecond(PEAKSAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.re.bedg.gz", combo=combo, file=samplecond(PEAKSAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}/{file}_peak_{type}.bed.gz", combo=combo, file=samplecond(PEAKSAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup']),
                    expand("PEAKS/{combo}/{file}_peak_seq_{type}.bed.gz", combo=combo, file=samplecond(PEAKSAMPLES, config), type=['sorted', 'sorted_unique', 'sorted_dedup', 'sorted_unique_dedup'])

include: "manipulate_genome.smk"

rule FindPeaks:
    input:  bam = expand("MAPPED/{scombo}/{{file}}_mapped_{{type}}.bam", scombo=scombo)
    output: peak = "PEAKS/{combo}/{file}_peak_{type}.bed.gz"
    log:    "LOGS/PEAKS/{combo}/{file}_findpeaks_macs_{type}.log"
    conda:  ""+PEAKENV+".yaml"
    threads: 1
    params: pairing = lambda wildcards: get_pairing(wildcards.file, wildcards.type, config, config.get('SAMPLES', get_samples_postprocess(config, 'PEAKS')), scombo),
            ppara = lambda wildcards: tool_params(wildcards.file, None, config, "PEAKS", PEAKENV)['OPTIONS'].get('FINDPEAKS', ""),
            peak = PEAKBIN,
            outdir = lambda wildcards, output: os.path.dirname(output.peak),
            outname = lambda wildcards: os.path.basename(wildcards.file)+'_peak_'+wildcards.type,
            mapmode = lambda wildcards: 'BAMPE' if paired == 'paired' and 'unique' not in wildcards.type else 'BAM'
    shell:  "set +o pipefail; export LC_ALL=C; if [[ -n \"$(samtools view {input.bam} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then {params.peak} callpeak -t {input.bam} {params.pairing} --outdir {params.outdir} -n {params.outname} -f {params.mapmode} {params.ppara} 2> {log} && gzip {params.outdir}/{params.outname}_peaks.narrowPeak 2>> {log} && ln -s {params.outname}_peaks.narrowPeak.gz {output.peak} 2>> {log}; else gzip < /dev/null > {output.peak}; echo \"File {input.bam} empty\" >> {log}; fi"


rule AddSequenceToPeak:
    input:  pk = rules.FindPeaks.output.peak,
            fa = rules.UnzipGenome.output.fa
    output: peak = "PEAKS/{combo}/{file}_peak_seq_{type}.bed.gz",
            pt = temp("PEAKS/{combo}/{file}_peak_chr_{type}.tmp"),
            ps = temp("PEAKS/{combo}/{file}_peak_seq_{type}.tmp")
    log:    "LOGS/PEAKS/{combo}/{file}_seq2peaks_{type}.log"
    conda:  "bedtools.yaml"
    threads: 1
    params: bins=BINS,
            sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
    shell:  "set +o pipefail; export LC_ALL=C; if [[ -n \"$(zcat {input.pk} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then zcat {input.pk} | perl -wlane '$F[0] = $F[0] =~ /^chr/ ? $F[0] : \"chr\".$F[0]; print join(\"\\t\",@F[0..5])' > {output.pt} && bedtools getfasta -fi {input.fa} -bed {output.pt} -name -tab -s -fullHeader -fo {output.ps} && cut -d$'\t' -f2 {output.ps}|paste -d$'\t' <(zcat {input.pk}) - |sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip  > {output.peak} 2> {log}; else gzip < /dev/null > {output.peak} && touch {output.pt} {output.ps}; fi"  # NEED TO GET RID OF SPACES AND WHATEVER IN HEADER

if ANNOPEAK is not None:
    rule AnnotatePeak:
        input:  "PEAKS/{combo}/{file}_peak_seq_{type}.bed.gz"
        output: "PEAKS/{combo}/{file}_peak_anno_{type}.bed.gz"
        log:    "LOGS/PEAKS/{combo}/{file}_annotatepeaks_{type}.log"
        conda:  "perl.yaml"
        threads: 1
        params: bins=BINS,
                anno = ANNOTATION
        shell:  "perl {params.bins}/Universal/AnnotateBed.pl -b {input} -a {params.anno} |gzip > {output} 2> {log}"

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
            re = rules.PeakToBedg.output.re
    output: fw = "PEAKS/{combo}/{file}_peak_{type}.fw.norm.bedg.gz",
            re = "PEAKS/{combo}/{file}_peak_{type}.re.norm.bedg.gz"
    log:    "LOGS/PEAKS/{combo}/{file}_ucscpeaknormalizebedgraph_{type}.log"
    conda:  "perl.yaml"
    threads: 1
    params: sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
    shell: "set +o pipefail; export LC_ALL=C; if [[ -n \"$(zcat {input.fw} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then scale=$(bc <<< \"scale=6;$(zcat {input.fw}|cut -f4|perl -wne '{{$x+=$_;}}END{{if ($x == 0){{$x=1}} print $x}}')/1000000\") perl -wlane '$sc=$ENV{{scale}};print join(\"\t\",@F[0..$#F-1]),\"\t\",$F[-1]/$sc' <(zcat {input.fw}) | sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.fw} 2> {log}; else gzip < /dev/null > {output.fw}; echo \"File {input.fw} empty\" >> {log}; fi && if [[ -n \"$(zcat {input.re} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then scale=$(bc <<< \"scale=6;$(zcat {input.re}|cut -f4|perl -wne '{{$x+=$_;}}END{{if ($x == 0){{$x=1}} print $x}}')/1000000\") perl -wlane '$sc=$ENV{{scale}};print join(\"\t\",@F[0..$#F-1]),\"\t\",$F[-1]/$sc' <(zcat {input.re}) | sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.re} 2> {log}; else gzip < /dev/null > {output.re}; echo \"File {input.re} empty\" >> {log}; fi"


### This step generates bigwig files for peaks which can then be copied to a web-browsable directory and uploaded to TRACKS via the track field
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
    log:    "LOGS/PEAKS/{combo}/{file}_peaktrack_{type}.log"
    conda:  "base.yaml"
    threads: MAXTHREAD
    params: bwdir = lambda wildcards: "TRACKS/{src}".format(src=SETS),
            bins = os.path.abspath(BINS),
            gen = REFDIR,  #lambda wildcards: os.path.basename(genomepath(wildcards.file, config)),
            options = '-n Peaks_'+str(PEAKENV)+' -s peaks -l TRACKS_peaks_'+str(PEAKENV)+' -b TRACKS_'+str(PEAKENV),
            uid = lambda wildcards: "{src}".format(src='TRACKS'+os.sep+'PEAKS'+os.sep+combo+SETS.replace(os.sep, '_'))
    shell: "echo -e \"{input.fw}\\n{input.re}\"|python3 {params.bins}/Analysis/GenerateTrackDb.py -i {params.uid} -e 1 -f STDIN -u '' -g {params.gen} {params.options} && touch {input.fw}.trackdone && touch {input.re}.trackdone 2> {log}"
