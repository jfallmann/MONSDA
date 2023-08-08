BINS = get_always('BINS')
PEAKSENV = get_always('PEAKSENV')
PEAKSBIN = get_always('PEAKSBIN')
REF = get_always('REFERENCE')
REFDIR = "${workflow.workDir}/../"+get_always('REFDIR')
SETS = get_always('SETS')
IP = get_always('IP')
SOFTPARAMS = get_always('scyphy_params_SOFTCLIP') ?: ''
PREPARAMS = get_always('scyphy_params_PREPROCESS') ?: ''
PEAKSPARAMS = get_always('scyphy_params_FINDPEAKS') ?: ''

include { UnzipGenome; UnzipGenome_no_us } from "manipulate_genome.nf"

process RemoveSoftclip{
    conda "$PEAKSENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bam.bai") > 0)      "MAPPED/${COMBO}/${CONDITION}/${file(filename).getName()}"                
        else if (filename.indexOf(".bam") > 0)      "MAPPED/${COMBO}/${CONDITION}/${file(filename).getName()}"                
        else if (filename.indexOf("log") > 0)        "LOGS/PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}_removeSoftclip.log"
    }

    input:
    path bam
    
    output:
    path "*_nosoftclip.bam", includeInputs:false, emit: bams
    path "*_nosoftclip.bam.bai", includeInputs:false, emit: bais
    path "*.log", emit: log
    
    script: 
    fn = file(bam).getSimpleName()
    fo = fn+'_nosoftclip.bam'
    fi = fn+'_nosoftclip.bam.bai'
    ol = "log"
    sortmem = '30%'
    
    """
    mkdir -p TMP/$fn; p=\$( dirname \$(realpath \"$bam\") ); ln -s \${p}/${fn}.bam.bai .; python $BINS/Analysis/RemoveSoftClip.py -f $REF -b $bam $SOFTPARAMS -o \'-\' | samtools sort -T TMP/$fn -o $fo --threads ${task.cpus} \'-\' 2>> $ol && samtools index $fo 2>> $ol && mv $ol ${fn}.log; rm -rf TMP/$fn
    """
}

process BamToBed{
    conda "bedtools.yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bed.gz") > 0)      "BED/${COMBO}/${CONDITION}/${file(filename).getName().replaceAll(/\Q_ext.bed.gz\E/,".bed.gz")}"                
        else if (filename.indexOf(".log") > 0)        "LOGS/PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}_bam2bed.log"
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
    sortmem = '30%'
    
    if (STRANDED == 'rf' || STRANDED == 'ISR'){
        """
        p=\$( dirname \$(realpath \"$bam\") ); ln -s \${p}/${fn}.bam.bai .;
        bedtools bamtobed -split -i $bam | sed 's/ /_/g' | perl -wl -a -F'\\t' -n -e '\$F[0] =~ s/\\s/_/g;if(\$F[3]=~/\\/1\$/){{if (\$F[5] eq \"+\"){{\$F[5] = \"-\"}}elsif(\$F[5] eq \"-\"){{\$F[5] = \"+\"}}}} print join(\"\\t\",@F[0..\$#F])' | sort --parallel=${task.cpus} -S $sortmem -T TMP -t\$'\\t' -k1,1 -k2,2n |gzip > $fo 2> $ol
        """
    }else{
        """
        p=\$( dirname \$(realpath \"$bam\") ); ln -s \${p}/${fn}.bam.bai .;
        bedtools bamtobed -split -i $bam | sed 's/ /_/g' | perl -wl -a -F'\\t' -n -e '\$F[0] =~ s/\\s/_/g;if(\$F[3]=~/\\/2\$/){{if (\$F[5] eq \"+\"){{\$F[5] = \"-\"}}elsif(\$F[5] eq \"-\"){{\$F[5] = \"+\"}}}} print join(\"\\t\",@F[0..\$#F])' | sort --parallel=${task.cpus} -S $sortmem -T TMP -t\$'\\t' -k1,1 -k2,2n |gzip > $fo 2> $ol
        """
    }
}


process ExtendBed{
    conda "perl.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bed.gz") > 0)      "BED/${COMBO}/${CONDITION}/${file(filename).getName()}"                
        else if (filename.indexOf(".log") > 0)        "LOGS/PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}_bam2bed.log"
    }

    input:
    path bedf

    output:
    path "*nosoftclip*.bed.gz", emit: nobedext, optional: true
    path "*.bed.gz", includeInputs:false, emit: bedext
    path "*.log", emit: log

    script: 
    bed = bedf[0]
    sizes = bedf[1] 

    
    fn = file(bed).getSimpleName().replaceAll(/\Q_mapped\E/,"_mapped_extended")
    of = fn+'.bed.gz'    
    opt = '-u 0'
    ol = fn+".log"
    sortmem = '30%'

    """        
    $BINS/Universal/ExtendBed.pl $opt -b $bed -o $of -g $sizes 2> $ol
    """    
}

process RevExtendBed{
    conda "perl.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bed.gz") > 0)      "BED/${COMBO}/${CONDITION}/${file(filename).getName()}"                
        else if (filename.indexOf(".log") > 0)        "LOGS/PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}_bam2bed.log"
    }

    input:
    path bedf

    output:
    path "*nosoftclip*.bed.gz", emit: nobedrev, optional: true
    path "*.bed.gz", includeInputs:false, emit: bedrev
    path "*.log", emit: log

    script: 
    bed = bedf[0]
    sizes = bedf[1] 

    
    fn = file(bed).getSimpleName().replaceAll(/\Q_mapped\E/,"_mapped_revtrimmed")
    of = fn+'.bed.gz'    
    opt = '-d 0'
    ol = fn+".log"
    sortmem = '30%'

    """        
    $BINS/Universal/ExtendBed.pl $opt -b $bed -o $of -g $sizes 2> $ol
    """    
}


process BedToBedg{
    conda "bedtools.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bedg.gz") > 0)      "PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}"                
        else if (filename.indexOf(".log") > 0)        "LOGS/PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}_ucscbedtobedgraph.log"
    }

    input:
    path bedf

    output:
    path "*nosoftclip*fw.bedg.gz", emit: nobedgf, optional: true
    path "*nosoftclip*re.bedg.gz", emit: nobedgr, optional: true
    path "*fw.bedg.gz", emit: bedgf
    path "*re.bedg.gz", emit: bedgr
    path "*.log", emit: log

    script: 
    bed = bedf[0]
    fai = bedf[1]
    sizes = bedf[2]

    fn = file(bed).getSimpleName()
     if (file(fn).getName().indexOf("_mapped_extended") > 0){
        fn = fn.replaceAll(/\Q_extended\E/,"")
        fw = fn+'_ext.fw.bedg.gz'
        fr = fn+'_ext.re.bedg.gz'
    } else if (file(fn).getName().indexOf("_mapped_revtrimmed") > 0){
        fn = file(fn).getName().replaceAll(/\Q_revtrimmed\E/,"")
        fw = fn+'_rev.fw.bedg.gz'
        fr = fn+'_rev.re.bedg.gz'
    }
    else{
        fn = file(fn).getName().replaceAll(/\Q_revtrimmed\E/,"").replaceAll(/\Q_extended\E/,"")
        fr = fn+'.re.bedg.gz'
        fw = fn+'.fw.bedg.gz'
    }
    ol = fn+".log"
    sortmem = '30%'

    """
    export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i $bed -bg -split -strand + -g $sizes | sort --parallel=${task.cpus} -S $sortmem -T TMP -t\$'\\t' -k1,1 -k2,2n |gzip > $fw 2> $ol && bedtools genomecov -i $bed -bg -split -strand - -g $sizes |sort -S $sortmem -T TMP -t\$'\\t' -k1,1 -k2,2n |gzip > $fr 2>> $ol
    """
}

process BedToBedgPeak{
    conda "bedtools.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bedg.gz") > 0)      "PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}"                
        else if (filename.indexOf(".log") > 0)        "LOGS/PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}_ucscbedtobedgraph.log"
    }

    input:
    path bedf

    output:
    path "*nosoftclip*fw.bedg.gz", emit: bedgf
    path "*nosoftclip*re.bedg.gz", emit: bedgr
    path "*.log", emit: log

    script: 
    bed = bedf[0]
    fai = bedf[1]
    sizes = bedf[2]


    fn = file(bed).getSimpleName()
    if (file(fn).getName().indexOf("_mapped_extended") > 0){
        fn = fn.replaceAll(/\Q_extended\E/,"")
        fw = fn+'_ext.fw.bedg.gz'
        fr = fn+'_ext.re.bedg.gz'
    } else if (file(fn).getName().indexOf("_mapped_revtrimmed") > 0){
        fn = file(fn).getName().replaceAll(/\Q_revtrimmed\E/,"")
        fw = fn+'_rev.fw.bedg.gz'
        fr = fn+'_rev.re.bedg.gz'
    }
    else{
        fn = file(fn).getName().replaceAll(/\Q_revtrimmed\E/,"").replaceAll(/\Q_extended\E/,"")
        fr = fn+'.re.bedg.gz'
        fw = fn+'.fw.bedg.gz'
    }
    ol = fn+".log"
    sortmem = '30%'

    """
    export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i $bed -bg -split -strand + -g $sizes | perl -wlane 'print join(\"\\t\",@F[0..2],\".\",\$F[3],\"+\")' > tosrt && cat tosrt| sort --parallel=${task.cpus} -S $sortmem -T TMP -t\$'\\t' -k1,1 -k2,2n |gzip > $fw 2> $ol && bedtools genomecov -i $bed -bg -split -strand - -g $sizes | perl -wlane 'print join(\"\\t\",@F[0..2],\".\",\$F[3],\"-\")' > tosrt && cat tosrt|sort -S $sortmem -T TMP -t\$'\\t' -k1,1 -k2,2n |gzip > $fr 2>> $ol
    """
}

process PreprocessPeaks{
    conda "perl.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bed.gz") > 0)      "PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}"                
        else if (filename.indexOf(".log") > 0)        "LOGS/PEAKS/${SCOMBO}/${CONDITION}/prepeak_${file(filename).getName()}.log"
    }

    input:
    path bedgf
    path bedgr

    output:
    path "*_prepeak*.bed.gz", emit: prepeak
    path "*.log", emit: log

    script: 
    of = file(bedgf).getSimpleName().replaceAll(/\Q_mapped\E/,"_prepeak").replaceAll(/\Q_rev\E/,"").replaceAll(/\Q_ext\E/,"")+".bed.gz"
    ol = file(bedgf).getSimpleName()+".log"
    sortmem = '30%'

    """        
    perl $BINS/Analysis/PreprocessPeaks.pl -p <(zcat $bedgf $bedgr| sort --parallel=${task.cpus} -S $sortmem -T TMP -t\$'\\t' -k1,1 -k3,3n -k2,2n -k6,6) $PREPARAMS | sort --parallel=${task.cpus} -S $sortmem -T TMP -t\$'\\t' -k1,1 -k3,3n -k2,2n -k6,6 |gzip > $of 2> $ol
    """    
}

process FindPeaks{
    conda "$PEAKSENV"+".yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bed.gz") > 0)      "PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}"                
        else if (filename.indexOf(".log") > 0)        "LOGS/PEAKS/${SCOMBO}/${CONDITION}/findpeaks_${file(filename).getName()}.log"
    }

    input:
    path bed

    output:
    path "*_peak*.bed.gz", emit: peak
    path "*.log", emit: log

    script: 
    of = file(bed).getName().replaceAll(/\Q_prepeak\E/,"_peak").replaceAll(/\Q_nosoftclip\E/,"").replaceAll(/\Q_rev.bed.gz\E/,".bed.gz").replaceAll(/\Q_ext.bed.gz\E/,".bed.gz")
    ol = file(of).getSimpleName()+".log"
    sortmem = '30%'

    """  
    export LC_ALL=C; if [[ -n \"\$(zcat $bed | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then $PEAKSBIN $PEAKSPARAMS <(zcat $bed|sort -t\$\'\\t\' -k1,1 -k3,3n -k2,2n -k6,6) 2> $ol|tail -n+2| sort --parallel=${task.cpus} -S $sortmem -T TMP -t\$\'\\t\' -k1,1 -k2,2n |grep -v \'nan\'| gzip > $of 2>> $ol; else gzip < /dev/null > $of; echo \"File $bed empty\" >> $ol; fi
    """    
}

process AddSequenceToPeak{
    conda "bedtools.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bed.gz") > 0)      "PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}"                
        else if (filename.indexOf(".log") > 0)        "LOGS/PEAKS/${SCOMBO}/${CONDITION}/findpeaks_${file(filename).getName()}.log"
    }

    input:
    path bed

    output:
    path "*_peak_seq*.bed.gz", includeInputs:false, emit: peak
    path "*.log", emit: log

    script: 
    pk = bed[0]
    fa = bed[1]
    of = file(pk).getName().replaceAll(/\Q_peak\E/,"_peak_seq")
    ol = file(pk).getSimpleName()+".log"
    sortmem = '30%'

    """  
    export LC_ALL=C; mkdir -p TMP; if [[ -n \"\$(zcat $pk | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then  zcat $pk | perl -wlane '\$F[0] = \$F[0] =~ /^chr/ ? \$F[0] : \"chr\".\$F[0]; print join(\"\\t\",@F[0..5])' > pktmp && bedtools getfasta -fi $fa -bed pktmp -name -tab -s -fullHeader -fo pkseqtmp && cut -d\$'\\t' -f2 pkseqtmp|sed 's/t/u/ig'|paste -d\$'\\t' <(zcat pktmp) - |sort --parallel=${task.cpus} -S $sortmem -T TMP -t\$'\\t' -k1,1 -k2,2n |gzip  > $of 2> $ol; else gzip < /dev/null > $of; echo \"File $pk empty\" >> $ol; fi
    """    
}

process PeakToBedg{
    conda "perl.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bedg.gz") > 0)      "PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}"                
        else if (filename.indexOf(".log") > 0)        "LOGS/PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}_peaktobedgraph.log"
    }

    input:
    path bedf

    output:
    path "*fw.bedg.gz", emit: bedgf
    path "*re.bedg.gz", emit: bedgr
    path "*.log", emit: log

    script: 
    bed = bedf[0]
    sizes = bedf[1]

    fn = file(bed).getSimpleName()
    fw = fn+'.fw.bedg.gz'
    fr = fn+'.re.bedg.gz'
    ol = fn+".log"
    sortmem = '30%'

    """
    mkdir -p TMP/fw TMP/re ; perl $BINS/Universal/Bed2Bedgraph.pl -f <(zcat $bed) -c $sizes -p peak -x tmp.fw.gz -y tmp.re.gz -a track 2>> $ol && zcat tmp.fw.gz | sort --parallel=${task.cpus} -S $sortmem -T TMP/fw -t\$'\\t' -k1,1 -k2,2n |gzip > $fw 2>> $ol && zcat tmp.re.gz |sort -S $sortmem -T TMP/re -t\$'\\t' -k1,1 -k2,2n |gzip > $fr 2>> $ol
    """
}


process NormalizeBedg{
    conda "perl.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".norm.bedg.gz") > 0)      "PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}"                
        else if (filename.indexOf(".log") > 0)        "LOGS/PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}_ucscpeaknormalizebedgraph.log"
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
    fw = fn+'.fw.norm.bedg.gz'
    fr = fn+'.re.norm.bedg.gz'
    ol = fn+".log"
    sortmem = '30%'
    
    """
    export LC_ALL=C; if [[ -n \"\$(zcat $bedgf | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then scale=\$(bc <<< \"scale=6;\$(zcat $bedgf|cut -f4|perl -wne '{\$x+=\$_;}END{if (\$x == 0){\$x=1} print \$x}')/1000000\") perl -wlane '\$sc=\$ENV{scale};print join(\"\\t\",@F[0..\$#F-1]),\"\\t\",\$F[-1]/\$sc' <(zcat $bedgf)| sort -S $sortmem -T TMP -t\$'\\t' -k1,1 -k2,2n |gzip > $fw 2> $ol; else gzip < /dev/null > $fw; echo \"File $bedgf empty\" >> $ol; fi && if [[ -n \"\$(zcat $bedgr | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then scale=\$(bc <<< \"scale=6;\$(zcat $bedgr|cut -f4|perl -wne '{\$x+=\$_;}END{if (\$x == 0){\$x=1} print \$x}')/1000000\") perl -wlane '\$sc=\$ENV{scale};print join(\"\\t\",@F[0..\$#F-1]),\"\\t\",\$F[-1]/\$sc' <(zcat $bedgr)| sort -S $sortmem -T TMP -t\$'\\t' -k1,1 -k2,2n|gzip > $fr 2> $ol; else gzip < /dev/null > $fr; echo \"File $bedgr empty\" >> $ol; fi
    """
}

process NormalizePeakBedg{
    conda "perl.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".norm.bedg.gz") > 0)      "PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}"                
        else if (filename.indexOf(".log") > 0)        "LOGS/PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}_ucscpeaknormalizebedgraph.log"
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
    fw = fn+'.fw.norm.bedg.gz'
    fr = fn+'.re.norm.bedg.gz'
    ol = fn+".log"
    sortmem = '30%'
    
    """
    export LC_ALL=C; if [[ -n \"\$(zcat $bedgf | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then scale=\$(bc <<< \"scale=6;\$(zcat $bedgf|cut -f4|perl -wne '{\$x+=\$_;}END{if (\$x == 0){\$x=1} print \$x}')/1000000\") perl -wlane '\$sc=\$ENV{scale};print join(\"\\t\",@F[0..\$#F-1]),\"\\t\",\$F[-1]/\$sc' <(zcat $bedgf)| sort -S $sortmem -T TMP -t\$'\\t' -k1,1 -k2,2n |gzip > $fw 2> $ol; else gzip < /dev/null > $fw; echo \"File $bedgf empty\" >> $ol; fi && if [[ -n \"\$(zcat $bedgr | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then scale=\$(bc <<< \"scale=6;\$(zcat $bedgr|cut -f4|perl -wne '{\$x+=\$_;}END{if (\$x == 0){\$x=1} print \$x}')/1000000\") perl -wlane '\$sc=\$ENV{scale};print join(\"\\t\",@F[0..\$#F-1]),\"\\t\",\$F[-1]/\$sc' <(zcat $bedgr)| sort -S $sortmem -T TMP -t\$'\\t' -k1,1 -k2,2n|gzip > $fr 2> $ol; else gzip < /dev/null > $fr; echo \"File $bedgr empty\" >> $ol; fi
    """
}

process PeakToTRACKS{
    conda "ucsc.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bw") > 0)      "TRACKS/PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}"                
        else if (filename.indexOf(".log") > 0)        "LOGS/TRACKS/PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}_peaktoucsc.log"
    }

    input:
    path bedgf
    path rest
    
    output:
    path "*.fw.bw", emit: bwf
    path "*.re.bw", emit: bwr
    path "*.log", emit: log

    script: 
    bedgr = rest[0]
    sizes = rest[1]
    fn = file(bedgf).getSimpleName()
    fw = fn+'.fw.bw'
    fr = fn+'.re.bw'
    ol = fn+".log"
    
    """
    export LC_ALL=C; if [[ -n \"\$(zcat $bedgf | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then zcat $bedgf > tmp && bedGraphToBigWig tmp $sizes $fw 2> $ol; else gzip < /dev/null > $fw; echo \"File $bedgf empty\" >> $ol; fi && if [[ -n \"\$(zcat $bedgr | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then zcat $bedgr > tmp && bedGraphToBigWig tmp $sizes $fr 2>> $ol; else gzip < /dev/null > $fr; echo \"File $bedgr empty\" >> $ol; fi
    """
}

process BedgToTRACKS{
    conda "ucsc.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bw") > 0)      "TRACKS/PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}"                
        else if (filename.indexOf(".log") > 0)        "LOGS/TRACKS/PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}_peaktoucsc.log"
    }

    input:
    path bedgf
    path rest
    
    output:
    path "*.fw.bw", emit: bwf
    path "*.re.bw", emit: bwr
    path "*.log", emit: log

    script: 
    bedgr = rest[0]
    sizes = rest[1]
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
        if (filename.indexOf(".txt") > 0)      "TRACKS/PEAKS/${file(filename).getName()}"
        else if (filename.indexOf(".log") > 0)        "LOGS/TRACKS/PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}"
    }

    input:
    path pkbw
    path mapbw

    output:
    path "*.txt", emit: trackdb
    path "*.log", emit: log

    script: 
    uid= SETS.replace(File.separator, "_")
    ol = uid+"_GenerateTrack_peaks.log"
    opt = '-n Peaks_'+"$PEAKSENV"+' -s peaks -l TRACKS_peaks_'+"$PEAKSENV"+' -b TRACKS_'+"$PEAKSENV"
    """
    bf=($pkbw[0]); br=($pkbw[1]); mf=($mapbw[0]); mr=($mapbw[1]); for i in \"\${!bf[@]}\";do fp=\${bf[\$i]}; rp=\${br[\$i]}; fm=\${mf[\$i]}; rm=\${mr[\$i]}; echo -e \"\$fp\\n\$rp\\n\$fm\\n\$rm\"|python3 $BINS/Analysis/GenerateTrackDb.py -i $uid -e 1 -f STDIN -u \"\" -g $REFDIR $opt 2>> $ol;done
    """
}


workflow PEAKS{ 
    take: collection

    main:
    
    MAPPEDSAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/${COMBO}/"+element+"_mapped_sorted*.bam" 
    }
    //MAPPEDSAMPLES = MAPPEDSAMPLES*.removeIf({it -> it ~= /nosoftclip/})
    //BAMINDICES = MAPPEDSAMPLES*.replaceAll(/.bam$/, ".bam.bai")

    mapsamples_ch = Channel.fromPath(MAPPEDSAMPLES.sort()).filter({ it=~/sorted.bam$|sorted_unique.bam$|sorted_dedup.bam$|sorted_unique_dedup.bam$/ })
    genomefile = Channel.fromPath(REF)

    //mapsamples_ch.subscribe {  println "BAM: $it"  }
    
    UnzipGenome(genomefile)
    RemoveSoftclip(mapsamples_ch)
    BamToBed(mapsamples_ch.concat(RemoveSoftclip.out.bams))
    ExtendBed(BamToBed.out.bed.combine(UnzipGenome.out.chromsize))
    RevExtendBed(BamToBed.out.bed.combine(UnzipGenome.out.chromsize))
    BedToBedg(ExtendBed.out.bedext.concat(RevExtendBed.out.bedrev).concat(BamToBed.out.bed).combine(UnzipGenome.out.index).combine(UnzipGenome.out.chromsize))
    NormalizeBedg(BedToBedg.out.bedgf, BedToBedg.out.bedgr)
    BedgToTRACKS(NormalizeBedg.out.bedgf, NormalizeBedg.out.bedgr.combine(UnzipGenome.out.chromsize))
    
    if (IP == 'iCLIP' || IP == 'CLIP'){
        BedToBedgPeak(ExtendBed.out.nobedext.combine(UnzipGenome.out.index.combine(UnzipGenome.out.chromsize)))
    } else if (IP == 'revCLIP'){
        BedToBedgPeak(RevExtendBed.out.nobedrev.combine(UnzipGenome.out.index.combine(UnzipGenome.out.chromsize)))
    }
    else{
        BedToBedgPeak(ExtendBed.out.nobedext.concat(RevExtendBed.out.nobedrev).combine(UnzipGenome.out.index.combine(UnzipGenome.out.chromsize)))
    }
    PreprocessPeaks(BedToBedgPeak.out.bedgf, BedToBedgPeak.out.bedgr)
    FindPeaks(PreprocessPeaks.out.prepeak)
    AddSequenceToPeak(FindPeaks.out.peak.combine(UnzipGenome.out.chromsize))
    PeakToBedg(FindPeaks.out.peak.combine(UnzipGenome.out.chromsize))
    NormalizePeakBedg(PeakToBedg.out.bedgf, PeakToBedg.out.bedgr)
    PeakToTRACKS(NormalizePeakBedg.out.bedgf, NormalizePeakBedg.out.bedgr.combine(UnzipGenome.out.chromsize))

    GenerateTrack(PeakToTRACKS.out.bwf.merge(PeakToTRACKS.out.bwr).collect(), BedgToTRACKS.out.bwf.merge(BedgToTRACKS.out.bwr).collect())

    emit:
    trackdb = GenerateTrack.out.trackdb
}
