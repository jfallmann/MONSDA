process UnzipGenome{
    conda "samtools.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".fa.fai") > 0)      "${REFDIR}/${file(filename).getName()}"
        else if (filename.indexOf(".fa") > 0)      "${REFDIR}/${file(filename).getName()}"
        else if (filename.indexOf(".chrom.sizes") > 0)      "${REFDIR}/${file(filename).getName()}"
        else if (filename == "log")      "LOGS/${SCOMBO}/${COMBO}_indexfa.log"
    }

    input:
    path ref

    output: 
    path "*.fa", emit: unzipped
    path "*.fa.fai", emit: index
    path "*.chrom.sizes", emit: chromsize
    path "log", emit: log

    script:
    fa = ref.getSimpleName()+".fa"
    fai = ref.getSimpleName()+".fa.fai"
    cs = ref.getSimpleName()+".chrom.sizes"

    """
    zcat $ref |perl -F'\\t' -wane 'if(\$_ =~ /^>/){{chomp(\$F[0]);print \"\\n\".\$F[0].\"\\n\"}} else{{(\$line=\$_)=~s/\\r[\\n]*/\\n/gm; chomp(\$line=\$_); print \$line}}' |tail -n+2 > $fa && $BINS/Preprocessing/indexfa.sh $fa 2> log && cut -f1,2 $fai > $cs
    """
}


process UnzipGenome_no_us{
    conda "samtools.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("_us.fa.fai") > 0)      "${REFDIR}/${file(filename).getName()}"
        else if (filename.indexOf("_us.fa") > 0)      "${REFDIR}/${file(filename).getName()}"
        else if (filename.indexOf("_us.chrom.sizes") > 0)      "${REFDIR}/${file(filename).getName()}"
        else if (filename == "log")      "LOGS/${SCOMBO}/${COMBO}_indexfa_us.log"
    }

    input:
    path ref

    output: 
    path "*.fa", emit: unzipped
    path "*.fa.fai", emit: index
    path "*.chrom.sizes", emit: chromsize
    path "log", emit: log

    script:
    fa = ref.getSimpleName()+"_us.fa"
    fai = ref.getSimpleName()+"_us.fa.fai"
    cs = ref.getSimpleName()+"_us.chrom.sizes"
    
    """
    zcat $ref |perl -F'\\t' -wane 'if(\$_ =~ /^>/){{\$F[0] = \$F[0] =~ /^>chr/ ? \$F[0] : \">chr\".substr(\$F[0],1) =~ s/_/./g;chomp(\$F[0]);print \"\\n\".\$F[0].\"\\n\"}} else{{(\$line=\$_)=~s/\\r[\\n]*/\\n/gm; chomp(\$line=\$_); print \$line}}' |tail -n+2 > $fa && $BINS/Preprocessing/indexfa.sh $fa 2> log && cut -f1,2 $fai > $cs
    """
}   