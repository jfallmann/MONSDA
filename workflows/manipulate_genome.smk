rule UnzipGenome:
    input:  ref = REFERENCE,
    output: fa = expand("{ref}.fa", ref=REFERENCE.replace('.fa.gz', '')),
            fai = expand("{ref}.fa.fai", ref=REFERENCE.replace('.fa.gz', '')),
            fas = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
    log:    expand("LOGS/PEAKS/{combo}/indexfa.log", combo=combo)
    conda:  "samtools.yaml"
    threads: 1
    params: bins = BINS
    shell:  "set +o pipefail; zcat {input[0]} |perl -F\\\\040 -wane 'if($_ =~ /^>/){{chomp($F[0]);print \"\\n\".$F[0].\"\\n\"}} else{{($line=$_)=~s/\\r[\\n]*/\\n/gm; chomp($line=$_); print $line}}' |tail -n+2 > {output.fa} && {params.bins}/Preprocessing/indexfa.sh {output.fa} 2> {log} && cut -f1,2 {output.fai} > {output.fas}"
    #shell:  "set +o pipefail; zcat {input[0]} |perl -F\\\\040 -wane 'if($_ =~ /^>/){{$F[0] = $F[0] =~ /^>chr/ ? $F[0] : \">chr\".substr($F[0],1);chomp($F[0]);print \"\\n\".$F[0].\"\\n\"}} else{{($line=$_)=~s/\\r[\\n]*/\\n/gm; chomp($line=$_); print $line}}' |tail -n+2 > {output.fa} && {params.bins}/Preprocessing/indexfa.sh {output.fa} 2> {log} && cut -f1,2 {output.fai} > {output.fas}"

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
    #shell:  "set +o pipefail; zcat {input[0]} |perl -F\\\\040 -wane 'if($_ =~ /^>/){{$F[0] = $F[0] =~ /^>chr/ ? $F[0] : \">chr\".substr($F[0],1))=~ s/\_/\./g;chomp($F[0]);print \"\\n\".$F[0].\"\\n\"}} else{{($line=$_)=~s/\\r[\\n]*/\\n/gm; chomp($line=$_); print $line}}' |tail -n+2 > {output.fa} && {params.bins}/Preprocessing/indexfa.sh {output.fa} 2> {log} && cut -f1,2 {output.fai} > {output.fas}"
