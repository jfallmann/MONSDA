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
    if not TRACKS_STRANDED or (TRACKS_STRANDED == 'fr' or TRACKS_STRANDED == 'ISF'):
        rule BamToBed:
            input:  expand("MAPPED/{scombo}/{{file}}_mapped_{{type}}.bam", scombo=scombo)
            output: "BED/{scombo}/{file}_mapped_{type}.bed.gz"
            log:    "LOGS/{scombo}/{file}/TRACKS/ucsc/bam2bed_{type}.log"
            conda:  "bedtools.yaml"
            container: "oras://jfallmann/monsda:bedtools"
            threads: 1
            params: sortmem = get_sortmem
            shell:  "bedtools bamtobed -split -i {input[0]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\\/2$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |sort -S {params.sortmem}G -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output[0]} 2> {log}"

    elif TRACKS_STRANDED and (TRACKS_STRANDED == 'rf' or TRACKS_STRANDED == 'ISR'):
        rule BamToBed:
            input:  expand("MAPPED/{scombo}/{{file}}_mapped_{{type}}.bam", scombo=scombo)
            output: "BED/{scombo}/{file}_mapped_{type}.bed.gz"
            log:    "LOGS/{scombo}/{file}/TRACKS/ucsc/bam2bed_{type}.log"
            conda:  "bedtools.yaml"
            container: "oras://jfallmann/monsda:bedtools"
            threads: 1
            params: sortmem = get_sortmem, flip = r'if($F[3]=~/\/1$/)' if TRACKS_PAIRED == 'paired' else 'if(1)'
            shell:  "bedtools bamtobed -split -i {input[0]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;{params.flip}{{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |sort -S {params.sortmem}G -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output[0]} 2> {log}"

include: "manipulate_genome.smk"

checklist = list()
for file in samplecond(SAMPLES, config):
    for type in ['sorted','unique'] if not rundedup else ['sorted', 'unique', 'sorted_dedup', 'sorted_unique_dedup']:
        for orient in ['fw','rw']:
            checklist.append(os.path.isfile(os.path.abspath('BED/'+scombo+file+'_mapped_'+type+'.'+orient+'.bedg.gz')))


if not all(checklist):
    rule BedToBedg:
        input:  bed = "BED/{combo}/{file}_mapped_{type}.bed.gz",
                fai = expand("{ref}.fa.fai", ref=REFERENCE.replace('.fa.gz', '').replace('.fa.bgz', '').removesuffix('.fa')),
                sizes = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', '').replace('.fa.bgz', '').removesuffix('.fa'))
        output: fw = "TRACKS/{combo}/{file}_mapped_{type}.fw.bedg.gz",
                re = "TRACKS/{combo}/{file}_mapped_{type}.re.bedg.gz"
        log:    "LOGS/{combo}/{file}/TRACKS/ucsc/{type}_ucscbedtobedgraph.log"
        conda:  "bedtools.yaml"
        container: "oras://jfallmann/monsda:bedtools"
        threads: 1
        params: sortmem = get_sortmem
        shell: "export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i {input.bed} -bg -split -strand + -g {input.sizes} | sort --parallel={threads} -S {params.sortmem}G -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.fw} 2> {log} && bedtools genomecov -i {input.bed} -bg -split -strand - -g {input.sizes} |sort --parallel={threads} -S {params.sortmem}G -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.re} 2>> {log}"

### This step normalizes the bedg files for comparison in the browser
rule NormalizeBedg:
    input:  bed = rules.BedToBedg.input.bed,
            fw = rules.BedToBedg.output.fw,
            re = rules.BedToBedg.output.re
    output: fw = "TRACKS/{combo}/{file}_mapped_{type}.fw.norm.bedg.gz",
            re = "TRACKS/{combo}/{file}_mapped_{type}.re.norm.bedg.gz"
    log:    "LOGS/{combo}/{file}/TRACKS/ucsc/{type}_ucscnormalizebedgraph.log"
    conda:  "perl.yaml"
    container: "oras://jfallmann/monsda:perl"
    threads: 1
    params: sortmem = get_sortmem
    shell: "export LC_ALL=C; n=$(zcat {input.bed}|cut -f4|sed 's/\\/[12]$//'|sort -S {params.sortmem}G -T TMP -u|wc -l); if [[ $n -eq 0 ]]; then n=1; fi; export scale=$(bc <<< \"scale=6;$n/1000000\"); if [[ -n \"$(zcat {input.fw} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then perl -wlane '$sc=$ENV{{scale}};print join(\"\t\",@F[0..$#F-1]),\"\t\",$F[-1]/$sc' <(zcat {input.fw})| sort --parallel={threads} -S {params.sortmem}G -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.fw} 2> {log}; else gzip < /dev/null > {output.fw}; echo \"File {input.fw} empty\" >> {log}; fi && if [[ -n \"$(zcat {input.re} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then perl -wlane '$sc=$ENV{{scale}};print join(\"\t\",@F[0..$#F-1]),\"\t\",$F[-1]/$sc' <(zcat {input.re})| sort --parallel={threads} -S {params.sortmem}G -T TMP -t$'\t' -k1,1 -k2,2n|gzip > {output.re} 2> {log}; else gzip < /dev/null > {output.re}; echo \"File {input.re} empty\" >> {log}; fi"

### This step generates bigwig files for bedg which can then be copied to a web-browsable directory and uploaded to TRACKS via the track field
rule BedgToTRACKS:
    input:  fw = rules.NormalizeBedg.output.fw,
            re = rules.NormalizeBedg.output.re,
            sizes = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', '').replace('.fa.bgz', '').removesuffix('.fa'))
    output: fw = "TRACKS/{combo}/{file}_mapped_{type}.fw.bw",
            re = "TRACKS/{combo}/{file}_mapped_{type}.re.bw",
            tfw = temp("TRACKS/{combo}/{file}_mapped_{type}.fw.tmp"),
            tre = temp("TRACKS/{combo}/{file}_mapped_{type}.re.tmp")
    log:    "LOGS/{combo}/{file}/TRACKS/ucsc/{type}_bedgtoucsc.log"
    conda:  "ucsc.yaml"
    container: "oras://jfallmann/monsda:ucsc"
    threads: 1
    priority: 10               # This should be finished before we generate tracks
    shell:  "export LC_ALL=C; if [[ -n \"$(zcat {input.fw} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then zcat {input.fw} > {output.tfw} && bedGraphToBigWig {output.tfw} {input.sizes} {output.fw} 2> {log}; else awk 'NR==1{{print $1\"\t0\t1\t0\"}}' {input.sizes} > {output.tfw} && bedGraphToBigWig {output.tfw} {input.sizes} {output.fw} 2> {log}; echo \"File {input.fw} empty\" >> {log}; fi && if [[ -n \"$(zcat {input.re} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then zcat {input.re} > {output.tre} && bedGraphToBigWig {output.tre} {input.sizes} {output.re} 2>> {log}; else awk 'NR==1{{print $1\"\t0\t1\t0\"}}' {input.sizes} > {output.tre} && bedGraphToBigWig {output.tre} {input.sizes} {output.re} 2>> {log}; echo \"File {input.re} empty\" >> {log}; fi"

rule GenerateTrack:
    input:  fw = rules.BedgToTRACKS.output.fw,
            re = rules.BedgToTRACKS.output.re
    output: "TRACKS/{combo}/{file}_mapped_{type}.fw.bw.trackdone",
            "TRACKS/{combo}/{file}_mapped_{type}.re.bw.trackdone"
    log:    "LOGS/{combo}/{file}/TRACKS/ucsc/track_{type}.log"
    conda:  "base.yaml"
    container: "oras://jfallmann/monsda:base"
    threads: MAXTHREAD
    params: bwdir = lambda wildcards: "TRACKS/{src}".format(src=SETS),
            bins = os.path.abspath(BINS),
            gen = REFDIR,#lambda wildcards: os.path.basename(genomepath(wildcards.file, config)),
            options = lambda wildcards: tool_params(wildcards.file, None, config, 'TRACKS', 'ucsc')['OPTIONS'].get('TRACKS', "") if os.sep in wildcards.file else tool_params(os.sep.join([wildcards.combo[1:], wildcards.file]), None, config, 'TRACKS', 'ucsc')['OPTIONS'].get('TRACKS', ""),
            uid = lambda wildcards: "{src}".format(src='TRACKS'+os.sep+SETS.replace(os.sep, '_'))
    shell: "echo -e \"{input.fw}\\n{input.re}\"|python3 {params.bins}/Analysis/GenerateTrackDb.py -i {params.uid} -e 1 -f STDIN -u '' -g {params.gen} {params.options} && touch {input.fw}\.trackdone && touch {input.re}.trackdone 2> {log}"
