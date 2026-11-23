BINS = get_always('BINS')
TRACKSENV = get_always('TRACKSENV')
TRACKSBIN = get_always('TRACKSBIN')
REF = get_always('REFERENCE')
REFDIR = "${workflow.workDir}/../"+get_always('REFDIR')
SETS = get_always('SETS')
TRACKSPARAMS = get_always('ucsc_params_UCSC') ?: get_always('ucsc_TRACKS_params_UCSC') ?: ''

TRACKBIN = 'ucsc'
TRACKENV = 'ucsc'

include { UnzipGenome; UnzipGenome_no_us } from "manipulate_genome.nf"

process BamToBed{
    conda "bedtools.yaml"
    container "oras://jfallmann/monsda:"+"bedtools"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bed.gz") > 0)      "BED/${SCOMBO}/${CONDITION}/${file(filename).getName()}"                
        else if (filename == ".log")        "LOGS/${SCOMBO}/${CONDITION}/TRACKS/ucsc/${file(filename).getName()}_bam2bed.log"
    }

    input:
    path bam
    
    output:
    path "*.bed.gz", emit: bed
    path "*.log", emit: log
    
    script: 
    fn = file(bam).getSimpleName()
    fo = fn+'.bed.gz'
    ol = fn+".log"
    def sortmem = Math.ceil(task.memory.giga as double) as int 
    
    if (TRACKSSTRANDED == 'rf' || TRACKSSTRANDED == 'ISR'){
        flip = TRACKSPAIRED == 'paired' ? 'if(\$F[3]=~/\\/1\$/)' : 'if(1)'
        """
        bedtools bamtobed -split -i $bam | sed 's/ /_/g' | perl -wl -a -F'\\t' -n -e '\$F[0] =~ s/\\s/_/g;$flip{{if (\$F[5] eq \"+\"){{\$F[5] = \"-\"}}elsif (\$F[5] eq \"-\"){{\$F[5] = \"+\"}}}} print join(\"\\t\",@F[0..\$#F])' | sort --parallel=${task.cpus} -S ${sortmem}G -T . -t\$'\\t' -k1,1 -k2,2n |gzip > $fo 2> $ol
        """
    }else{
        """
        bedtools bamtobed -split -i $bam | sed 's/ /_/g' | perl -wl -a -F'\\t' -n -e '\$F[0] =~ s/\\s/_/g;if(\$F[3]=~/\\/2\$/){{if (\$F[5] eq \"+\"){{\$F[5] = \"-\"}}elsif(\$F[5] eq \"-\"){{\$F[5] = \"+\"}}}} print join(\"\\t\",@F[0..\$#F])' | sort --parallel=${task.cpus} -S ${sortmem}G -T . -t\$'\\t' -k1,1 -k2,2n |gzip > $fo 2> $ol
        """
    }
}


process BedToBedg{
    conda "bedtools.yaml"
    container "oras://jfallmann/monsda:"+"bedtools"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bedg.gz") > 0)      "TRACKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}"                
        else if (filename == ".log")        "LOGS/${SCOMBO}/${CONDITION}/TRACKS/ucsc/${file(filename).getName()}_ucscbedtobedgraph.log"
    }

    input:
    path bedf
    path sizes

    output:
    tuple path(bedf), path("*fw.bedg.gz"), path("*re.bedg.gz"), emit: bedg
    path "*.log", emit: log

    script: 
    bed = bedf

    fn = file(bed).getSimpleName()
    fw = fn+'.fw.bedg.gz'
    fr = fn+'.re.bedg.gz'
    ol = fn+".log"
    def sortmem = Math.ceil(task.memory.giga as double) as int 

    """
    export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i $bed -bg -split -strand + -g $sizes | sort --parallel=${task.cpus} -S ${sortmem}G -T . -t\$'\\t' -k1,1 -k2,2n |gzip > $fw 2> $ol && bedtools genomecov -i $bed -bg -split -strand - -g $sizes |sort -S ${sortmem}G -T . -t\$'\\t' -k1,1 -k2,2n |gzip > $fr 2>> $ol
    """
}


process NormalizeBedg{
    conda "perl.yaml"
    container "oras://jfallmann/monsda:"+"perl"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".norm.bedg.gz") > 0)      "TRACKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}"                
        else if (filename == ".log")        "LOGS/${SCOMBO}/${CONDITION}/TRACKS/ucsc/${file(filename).getName()}_ucscnormalizebedgraph.log"
    }

    input:
    tuple path(bedf), path(bedgf), path(bedgr)

    output:
    tuple path("*.fw.norm.bedg.gz"), path("*.re.norm.bedg.gz"), emit: bedg
    path "*.log", emit: log

    script: 
    fn = file(bedgf).getSimpleName()
    fw = fn+'.fw.norm.bedg.gz'
    fr = fn+'.re.norm.bedg.gz'
    ol = fn+".log"
    def sortmem = Math.ceil(task.memory.giga as double) as int 
    
    """
    export LC_ALL=C; n=\$(zcat $bedf|cut -f4|sed 's/\\/[12]\$//'|sort -S ${sortmem}G -T . -u|wc -l); if [[ \$n -eq 0 ]]; then n=1; fi; export scale=\$(bc <<< \"scale=6;\$n/1000000\"); if [[ -n \"\$(zcat $bedgf | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then perl -wlane '\$sc=\$ENV{scale};print join(\"\\t\",@F[0..\$#F-1]),\"\\t\",\$F[-1]/\$sc' <(zcat $bedgf)| sort -S ${sortmem}G -T . -t\$'\\t' -k1,1 -k2,2n |gzip > $fw 2> $ol; else gzip < /dev/null > $fw; echo \"File $bedgf empty\" >> $ol; fi && if [[ -n \"\$(zcat $bedgr | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then perl -wlane '\$sc=\$ENV{scale};print join(\"\\t\",@F[0..\$#F-1]),\"\\t\",\$F[-1]/\$sc' <(zcat $bedgr)| sort -S ${sortmem}G -T . -t\$'\\t' -k1,1 -k2,2n|gzip > $fr 2> $ol; else gzip < /dev/null > $fr; echo \"File $bedgr empty\" >> $ol; fi
    """
}

process BedgToTRACKS{
    conda "ucsc.yaml"
    container "oras://jfallmann/monsda:"+"ucsc"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bw") > 0)      "TRACKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}"                
        else if (filename == ".log")        "LOGS/${SCOMBO}/${CONDITION}/TRACKS/ucsc/${file(filename).getName()}_bedgtoucsc.log"
    }

    input:
    tuple path(bedgf), path(bedgr)
    path sizes
    
    output:
    tuple path("*.fw.bw"), path("*.re.bw"), emit: bw
    path "*.log", emit: log

    script: 
    fn = file(bedgf).getSimpleName()
    fw = fn+'.fw.bw'
    fr = fn+'.re.bw'
    ol = fn+".log"
    
    """
    export LC_ALL=C; if [[ -n \"\$(zcat $bedgf | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then zcat $bedgf > tmp && bedGraphToBigWig tmp $sizes $fw 2> $ol; else awk 'NR==1{print \$1\"\\t0\\t1\\t0\"}' $sizes > tmp && bedGraphToBigWig tmp $sizes $fw 2> $ol; echo \"File $bedgf empty\" >> $ol; fi && if [[ -n \"\$(zcat $bedgr | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then zcat $bedgr > tmp && bedGraphToBigWig tmp $sizes $fr 2>> $ol; else awk 'NR==1{print \$1\"\\t0\\t1\\t0\"}' $sizes > tmp && bedGraphToBigWig tmp $sizes $fr 2>> $ol; echo \"File $bedgr empty\" >> $ol; fi
    """
}

process GenerateTrack{
    conda "base.yaml"
    container "oras://jfallmann/monsda:"+"base"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".txt") > 0)      "TRACKS/${file(filename).getName()}"
        else if (filename == ".log")        "LOGS/${SCOMBO}/${CONDITION}/TRACKS/ucsc/${file(filename).getName()}"
    }

    input:
    path bw

    output:
    path "*.txt", emit: trackdb
    path "*.log", emit: log

    script: 
    uid= SETS.replace(File.separator, "_")
    ol = uid+"_GenerateTrack.log"
    """
    mkdir -p LOGS;touch LOGS/MONSDA.log; printf '%s\\n' *.bw|python3 $BINS/Analysis/GenerateTrackDb.py -i $uid -e 1 -f STDIN -u "TRACKS/${SCOMBO}/${CONDITION}" -g $REFDIR $TRACKSPARAMS 2>> $ol
    """
}

workflow TRACKS{ 
    take: collection

    main:
    
    MAPPEDSAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/${COMBO}/"+element+"_mapped_sorted{,_*}.bam"
    }

    mapsamples_ch = Channel.fromPath(MAPPEDSAMPLES.sort())
    genomefile = Channel.fromPath(REF)

    UnzipGenome(genomefile)
    chromsizes_ch = UnzipGenome.out.chromsize.first()
    BamToBed(mapsamples_ch)
    BedToBedg(BamToBed.out.bed, chromsizes_ch)
    NormalizeBedg(BedToBedg.out.bedg)
    BedgToTRACKS(NormalizeBedg.out.bedg, chromsizes_ch)
    GenerateTrack(BedgToTRACKS.out.bw.collect())

    emit:
    trackdb = GenerateTrack.out.trackdb
}
