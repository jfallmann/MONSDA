BINS = get_always('BINS')
TRACKSENV = get_always('TRACKSENV')
TRACKSBIN = get_always('TRACKSBIN')
REF = get_always('REF')
REFDIR = get_always('REFDIR')
ANNO = get_always('ANNO')
TRACKSPARAMS = get_always('ucsc_TRACKS_params_UCSC') ?: ''

TRACKBIN = 'ucsc'
TRACKENV = 'ucsc'

include { UnzipGenome; UnzipGenome_no_us } from './manipulate_genome.nf'

process BamToBed{
    conda "bedtools.yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bed.gz") > 0)      "BED/${SCOMBO}/${file(filename).getName()}"                
        else if (filename == ".log")        "LOGS/TRACKS/${SCOMBO}/${file(filename).getName()}_bam2bed.log"
    }

    input:
    path bam

    output:
    path "*.bed.gz", emit: bed
    path "*.log", emit: log

    script: 
    anno = fls[0]
    reads = fls[1]       
    fn = file(reads).getSimpleName()
    fo = fn+'.bed.gz'
    ol = fn+".log"
    sortmem = '30%'
    if (PAIRED == 'paired'){
        pair = "-p"
    }
    else{
        pair= ""
    }
    if (STRANDED == 'rf' || STRANDED == 'ISR'){
        """
        bedtools bamtobed -split -i $bam |sed 's/ /\\_/g'|perl -wl -a -F\'\\t\' -n -e '\$F[0] =~ s/\s/_/g;if(\$F[3]=~/\/2\$/){{if (\$F[5] eq \"+\"){{\$F[5] = \"-\"}}elsif(\$F[5] eq \"-\"){{\$F[5] = \"+\"}}}} print join(\"\t\",@F[0..\$#F])' |sort -S $sortmem -T TMP -t\$'\t' -k1,1 -k2,2n |gzip > $fo 2> $ol"
        """        
    }else{
        """
        bedtools bamtobed -split -i $bam |sed 's/ /\\_/g'|perl -wl -a -F\'\\t\' -n -e '\$F[0] =~ s/\s/_/g;if(\$F[3]=~/\/2\$/){{if (\$F[5] eq \"+\"){{\$F[5] = \"-\"}}elsif(\$F[5] eq \"-\"){{\$F[5] = \"+\"}}}} print join(\"\t\",@F[0..\$#F])' |sort --parallel=$THREADS -S $sortmem -T TMP -t\$'\t' -k1,1 -k2,2n |gzip > $fo 2> $ol
        """     
    }
}


process BedToBedg{
    conda "bedtools.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bedg.gz") > 0)      "TRACKS/${SCOMBO}/${file(filename).getName()}"                
        else if (filename == ".log")        "LOGS/TRACKS/${SCOMBO}/${file(filename).getName()}_ucscbedtobedgraph.log"
    }

    input:
    path bed
    path fai
    path sizes

    output:
    path "*fw.bedg.gz", emit: bedgf
    path "*re.bedg.gz", emit: bedgr
    path "*.log", emit: log

    script: 
    fn = file(bed).getSimpleName()
    fw = fn+'.fw.bedg.gz'
    fr = fn+'.re.bedg.gz'
    ol = fn+".log"
    sortmem = '30%'

    """
    export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i $bed -bg -split -strand + -g $sizes | sort --parallel=$THREADS -S $sortmem -T TMP -t\$'\t' -k1,1 -k2,2n |gzip > $fw 2> $ol && bedtools genomecov -i $bed -bg -split -strand - -g $sizes |sort -S $sortmem -T TMP -t\$'\t' -k1,1 -k2,2n |gzip > $fr 2>> $ol
    """
}


process NormalizeBedg{
    conda "perl.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".norm.bedg.gz") > 0)      "TRACKS/${SCOMBO}/${file(filename).getName()}"                
        else if (filename == ".log")        "LOGS/TRACKS/${SCOMBO}/${file(filename).getName()}_ucscnormalizebedgraph.log"
    }

    input:
    path bedgf
    path bedgr

    output:
    path "*.fw.norm.bedg.gz", emit: bedgf
    path "*.re.norm.bedg.gz", emit: bedgr
    path "*.log", emit: log

    script: 
    fn = file(bedgf).getSimpleName()
    fw = fn+'.norm.fw.bedg.gz'
    fr = fn+'.norm.re.bedg.gz'
    ol = fn+".log"
    sortmem = '30%'
    
    """
    export LC_ALL=C; if [[ -n \"\$(zcat $bedgf | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then scale=\$(bc <<< \"scale=6;\$(zcat $bedgf|cut -f4|perl -wne '{\$x+=\$_;}END{if (\$x == 0){\$x=1} print \$x}')/1000000\") perl -wlane '\$sc=\$ENV{scale};print join(\"\t\",@F[0..\$#F-1]),\"\t\",\$F[-1]/\$sc' <(zcat $bedgf)| sort -S $sortmem -T TMP -t\$'\t' -k1,1 -k2,2n |gzip > $fw 2> $ol; else gzip < /dev/null > $fw; echo \"File $bedgf empty\" >> $ol; fi && if [[ -n \"\$(zcat $bedgr | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then scale=\$(bc <<< \"scale=6;\$(zcat $bedgr|cut -f4|perl -wne '{\$x+=\$_;}END{if (\$x == 0){\$x=1} print \$x}')/1000000\") perl -wlane '\$sc=\$ENV{scale};print join(\"\t\",@F[0..\$#F-1]),\"\t\",\$F[-1]/\$sc' <(zcat $bedgr)| sort -S $sortmem -T TMP -t\$'\t' -k1,1 -k2,2n|gzip > $bedgr 2> $ol; else gzip < /dev/null > $bedgr; echo \"File $bedgr empty\" >> $ol; fi
    """
}

process BedgToTRACKS{
    conda "ucsc.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bw") > 0)      "TRACKS/${SCOMBO}/${file(filename).getName()}"                
        else if (filename == ".log")        "LOGS/TRACKS/${SCOMBO}/${file(filename).getName()}_bedgtoucsc.log"
    }

    input:
    path bedgf
    path bedgr
    path sizes

    output:
    path "*.fw.bw", emit: bwf
    path "*.re.bw", emit: bwr
    path "*.log", emit: log

    script: 
    fn = file(bedgf).getSimpleName()
    fw = fn+'.fw.bw'
    fr = fn+'.re.bw'
    ol = fn+".log"
    
    """
    export LC_ALL=C; if [[ -n \"\$(zcat $bedgf | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then zcat $bedgf > tmp && bedGraphToBigWig tmp $sizes $fw 2> $ol; else gzip < /dev/null > $fw; echo \"File $bedgf empty\" >> $ol; fi && if [[ -n \"\$(zcat $bedgr | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then zcat $bedgr > tmp && bedGraphToBigWig tmp $sizes $fr 2>> $ol; else gzip < /dev/null > $fr; echo \"File $bedgr empty\" >> $ol; fi
    """
}

process GenerateTrack{
    conda "base.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".txt") > 0)      "TRACKS/${file(filename).getName()}"
        else if (filename == ".log")        "LOGS/TRACKS/${SCOMBO}/${file(filename).getName()}_track.log"
    }

    input:
    path bwf
    path bwr

    output:
    path "*.txt", emit: trackdb
    path "*.log", emit: log

    script: 
    fn = file(bwf).getSimpleName()
    ol = fn+".log"
    src='TRACKS'+os.sep+$SETS.replace(os.sep, '_')
    
    """
    echo -e \"$bwf\\n$bwr\"|python3 $BINS/Analysis/GenerateTrackDb.py -i $SETS -e 1 -f STDIN -u '' -g $REFDIR $TRACKSPARAMS 2> $ol
    """
}

workflow DE{ 
    take: collection

    main:
    
    MAPPEDSAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/${COMBO}/"+element+"_mapped_sorted_unique.bam"
    }

    mapsamples_ch = Channel.fromPath(MAPPEDSAMPLES)
    mapsamples_ch.subscribe {  println "MAP: $it \t COMBO: ${COMBO} SCOMBO: ${SCOMBO} LONG: ${LONGSAMPLES}"  }
    annofile = Channel.fromPath(DEANNO)
    //annofile.subscribe {  println "ANNO: $it \t COMBO: ${COMBO} SCOMBO: ${SCOMBO} LONG: ${LONGSAMPLES}"  }

    featurecount_edger(annofile.combine(mapsamples_ch.collate(1)))
    prepare_count_table(featurecount_edger.out.fc_cts.collect())
    run_edger(prepare_count_table.out.counts, prepare_count_table.out.anno, annofile)
    filter_significant(run_edger.out.tbls)
    create_summary_snippet(run_edger.out.tbls.concat(run_edger.out.figs.concat(run_edger.out.session)).collect())
    collect_edger(filter_significant.out.sigtbls.collect())

    emit:
    tbls = run_edger.out.tbls
    sigtbls = filter_significant.out.sigtbls
    figs = run_edger.out.figs
    snps = create_summary_snippet.out.snps
}

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
            log:    "LOGS/TRACKS/{scombo}/{file}bam2bed_{type}.log"
            conda:  "bedtools.yaml"
            threads: 1
            params: sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
            shell:  "bedtools bamtobed -split -i {input[0]} |sed 's/ /\\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/2$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |sort -S {params.sortmem}% -T TMP -t\$'\t' -k1,1 -k2,2n |gzip > {output[0]} 2> {log}"

    elif stranded and (stranded == 'rf' or stranded == 'ISR'):
        rule BamToBed:
            input:  expand("MAPPED/{scombo}/{{file}}_mapped_{{type}}.bam", scombo=scombo)
            output: "BED/{scombo}/{file}_mapped_{type}.bed.gz"
            log:    "LOGS/PEAKS/{scombo}/{file}bam2bed_{type}.log"
            conda:  "bedtools.yaml"
            threads: 1
            params: sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
            shell:  "bedtools bamtobed -split -i {input[0]} |sed 's/ /\\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/1$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |sort -S {params.sortmem}% -T TMP -t\$'\t' -k1,1 -k2,2n |gzip > {output[0]} 2> {log}"

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
        shell: "export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i {input.bed} -bg -split -strand + -g {input.sizes} | sort --parallel={threads} -S {params.sortmem}% -T TMP -t\$'\t' -k1,1 -k2,2n |gzip > {output.fw} 2> {log} && bedtools genomecov -i {input.bed} -bg -split -strand - -g {input.sizes} |sort --parallel={threads} -S {params.sortmem}% -T TMP -t\$'\t' -k1,1 -k2,2n |gzip > {output.re} 2>> {log}"

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
    shell: "export LC_ALL=C; if [[ -n \"$(zcat {input.fw} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then scale=$(bc <<< \"scale=6;$(zcat {input.fw}|cut -f4|perl -wne '{{$x+=$_;}}END{{if ($x == 0){{$x=1}} print $x}}')/1000000\") perl -wlane '$sc=$ENV{{scale}};print join(\"\t\",@F[0..$#F-1]),\"\t\",$F[-1]/$sc' <(zcat {input.fw})| sort --parallel={threads} -S {params.sortmem}% -T TMP -t\$'\t' -k1,1 -k2,2n |gzip > {output.fw} 2> {log}; else gzip < /dev/null > {output.fw}; echo \"File {input.fw} empty\" >> {log}; fi && if [[ -n \"$(zcat {input.re} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then scale=$(bc <<< \"scale=6;$(zcat {input.re}|cut -f4|perl -wne '{{$x+=$_;}}END{{if ($x == 0){{$x=1}} print $x}}')/1000000\") perl -wlane '$sc=$ENV{{scale}};print join(\"\t\",@F[0..$#F-1]),\"\t\",$F[-1]/$sc' <(zcat {input.re})| sort --parallel={threads} -S {params.sortmem}% -T TMP -t\$'\t' -k1,1 -k2,2n|gzip > {output.re} 2> {log}; else gzip < /dev/null > {output.re}; echo \"File {input.re} empty\" >> {log}; fi"

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
