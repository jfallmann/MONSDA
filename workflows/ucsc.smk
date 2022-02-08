#combo = combo.split(os.sep)[0]+os.sep  # Unlike countreads we only have one env/bin to run ucsc tracks so we do not need scombo but can directly overwrite combo

wildcard_constraints:
    combo=scombo                # Only needed for ucsc.smk

rule themall:
    input: expand("TRACKS/{combo}/{file}_mapped_{type}.{orient}.bw.trackdone", file=samplecond(SAMPLES, config), type=["sorted", "sorted_unique"], orient=['fw','re'], combo=scombo) if not rundedup else expand("TRACKS/{combo}/{file}_mapped_{type}.{orient}.bw.trackdone", file=samplecond(SAMPLES, config), type=["sorted", "sorted_unique", "sorted_dedup", "sorted_unique_dedup"], orient=['fw', 're'], combo=scombo)

checklist = list()
for file in samplecond(SAMPLES, config):
    checktype = ['sorted', 'unique'] if not rundedup else ['sorted', 'unique', 'sorted_dedup', 'sorted_unique_dedup']
    for type in checktype:
        checklist.append(os.path.isfile(os.path.abspath('BED/'+scombo+file+'_mapped_'+type+'.bed.gz')) and not os.path.islink(os.path.abspath('BED/'+scombo+file+'_mapped_'+type+'.bed.gz')))

if not all(checklist):
    if not stranded or (stranded == 'fr' or stranded == 'ISF'):
        rule BamToBed:
            input:  expand("MAPPED/{scombo}/{{file}}_mapped_{{type}}.bam", scombo=scombo)
            output: "BED/{scombo}/{file}_mapped_{type}.bed.gz"
            log:    "LOGS/PEAKS/{scombo}/{file}bam2bed_{type}.log"
            conda:  "bedtools.yaml"
            threads: 1
            params: sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
            shell:  "bedtools bamtobed -split -i {input[0]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/2$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |sort -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output[0]} 2> {log}"

    elif stranded and (stranded == 'rf' or stranded == 'ISR'):
        rule BamToBed:
            input:  expand("MAPPED/{scombo}/{{file}}_mapped_{{type}}.bam", scombo=scombo)
            output: "BED/{scombo}/{file}_mapped_{type}.bed.gz"
            log:    "LOGS/PEAKS/{scombo}/{file}bam2bed_{type}.log"
            conda:  "bedtools.yaml"
            threads: 1
            params: sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
            shell:  "bedtools bamtobed -split -i {input[0]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/1$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |sort -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output[0]} 2> {log}"

include: "manipulate_genome.smk"

checklist = list()
for file in samplecond(SAMPLES, config):
    for type in ['sorted','unique'] if not rundedup else ['sorted', 'unique', 'sorted_dedup', 'sorted_unique_dedup']:
        for orient in ['fw','rw']:
            checklist.append(os.path.isfile(os.path.abspath('BED/'+scombo+file+'_mapped_'+type+'.'+orient+'.bedg.gz')))


if not all(checklist):
    rule BedToBedg:
        input:  bed = "BED/{combo}/{file}_mapped_{type}.bed.gz",
                fai = expand("{ref}.fa.fai", ref=REFERENCE.replace('.fa.gz', '')),
                sizes = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
        output: fw = "TRACKS/{combo}/{file}_mapped_{type}.fw.bedg.gz",
                re = "TRACKS/{combo}/{file}_mapped_{type}.re.bedg.gz"
        log:    "LOGS/TRACKS/{combo}/{file}_{type}_ucscbedtobedgraph.log"
        conda:  "bedtools.yaml"
        threads: 1
        params: sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
        shell: "export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i {input.bed} -bg -split -strand + -g {input.sizes} | sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.fw} 2> {log} && bedtools genomecov -i {input.bed} -bg -split -strand - -g {input.sizes} |sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.re} 2>> {log}"

### This step normalizes the bedg files for comparison in the browser
rule NormalizeBedg:
    input:  fw = rules.BedToBedg.output.fw,
            re = rules.BedToBedg.output.re
    output: fw = "TRACKS/{combo}/{file}_mapped_{type}.fw.norm.bedg.gz",
            re = "TRACKS/{combo}/{file}_mapped_{type}.re.norm.bedg.gz"
    log:    "LOGS/TRACKS/{combo}/{file}_{type}_ucscnormalizebedgraph.log"
    conda:  "perl.yaml"
    threads: 1
    params: sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
    shell: "export LC_ALL=C; if [[ -n \"$(zcat {input.fw} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then scale=$(bc <<< \"scale=6;$(zcat {input.fw}|cut -f4|perl -wne '{{$x+=$_;}}END{{if ($x == 0){{$x=1}} print $x}}')/1000000\") perl -wlane '$sc=$ENV{{scale}};print join(\"\t\",@F[0..$#F-1]),\"\t\",$F[-1]/$sc' <(zcat {input.fw})| sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.fw} 2> {log}; else gzip < /dev/null > {output.fw}; echo \"File {input.fw} empty\" >> {log}; fi && if [[ -n \"$(zcat {input.re} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then scale=$(bc <<< \"scale=6;$(zcat {input.re}|cut -f4|perl -wne '{{$x+=$_;}}END{{if ($x == 0){{$x=1}} print $x}}')/1000000\") perl -wlane '$sc=$ENV{{scale}};print join(\"\t\",@F[0..$#F-1]),\"\t\",$F[-1]/$sc' <(zcat {input.re})| sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n|gzip > {output.re} 2> {log}; else gzip < /dev/null > {output.re}; echo \"File {input.re} empty\" >> {log}; fi"

### This step generates bigwig files for bedg which can then be copied to a web-browsable directory and uploaded to TRACKS via the track field
rule BedgToTRACKS:
    input:  fw = rules.NormalizeBedg.output.fw,
            re = rules.NormalizeBedg.output.re,
            sizes = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
    output: fw = "TRACKS/{combo}/{file}_mapped_{type}.fw.bw",
            re = "TRACKS/{combo}/{file}_mapped_{type}.re.bw",
            t1 = temp("TRACKS/{combo}/{file}_mapped_{type}.fw.tmp"),
            t2 = temp("TRACKS/{combo}/{file}_mapped_{type}.re.tmp")
    log:    "LOGS/TRACKS/{combo}/{file}_{type}_bedgtoucsc.log"
    conda:  "ucsc.yaml"
    threads: 1
    priority: 10               # This should be finished before we generate tracks
    shell:  "export LC_ALL=C; if [[ -n \"$(zcat {input.fw} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then zcat {input.fw} > {output.t1} && bedGraphToBigWig {output.t1} {input.sizes} {output.fw} 2> {log}; else touch {output.t1}; gzip < /dev/null > {output.fw}; echo \"File {input.fw} empty\" >> {log}; fi && if [[ -n \"$(zcat {input.re} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then zcat {input.re} > {output.t2} && bedGraphToBigWig {output.t2} {input.sizes} {output.re} 2>> {log}; else touch {output.t2}; gzip < /dev/null > {output.re}; echo \"File {input.re} empty\" >> {log}; fi"

rule GenerateTrack:
    input:  fw = rules.BedgToTRACKS.output.fw,
            re = rules.BedgToTRACKS.output.re
    output: "TRACKS/{combo}/{file}_mapped_{type}.fw.bw.trackdone",
            "TRACKS/{combo}/{file}_mapped_{type}.re.bw.trackdone"
    log:    "LOGS/TRACKS/{combo}/{file}_track_{type}.log"
    conda:  "base.yaml"
    threads: MAXTHREAD
    params: bwdir = lambda wildcards: "TRACKS/{src}".format(src=SETS),
            bins = os.path.abspath(BINS),
            gen = REFDIR,#lambda wildcards: os.path.basename(genomepath(wildcards.file, config)),
            options = lambda wildcards: tool_params(wildcards.file, None, config, 'TRACKS', 'ucsc')['OPTIONS'].get('TRACKS', "") if os.sep in wildcards.file else tool_params(os.sep.join([wildcards.combo[1:], wildcards.file]), None, config, 'TRACKS', 'ucsc')['OPTIONS'].get('TRACKS', ""),
            uid = lambda wildcards: "{src}".format(src='TRACKS'+os.sep+SETS.replace(os.sep, '_'))
    shell: "echo -e \"{input.fw}\\n{input.re}\"|python3 {params.bins}/Analysis/GenerateTrackDb.py -i {params.uid} -e 1 -f STDIN -u '' -g {params.gen} {params.options} && touch {input.fw}\.trackdone && touch {input.re}.trackdone 2> {log}"
