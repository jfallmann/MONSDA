rule themall:
    input: expand("BED/{combo}/{file}_anno_seq_{type}_merged.bed.gz", file=samplecond(SAMPLES, config), type=['sorted','unique'])

checklist = list()
checklist2 = list()
for file in samplecond(SAMPLES, config):
    for type in ['sorted','unique']:
        checklist.append(os.path.isfile(os.path.abspath('TRACKS/'+file+'_mapped_'+type+'.bed.gz')) and not os.path.islink(os.path.abspath('TRACKS/'+file+'_mapped_'+type+'.bed.gz')))
        checklist2.append(os.path.isfile(os.path.abspath('PEAKS/'+file+'_mapped_'+type+'.bed.gz')) and not os.path.islink(os.path.abspath('PEAKS/'+file+'_mapped_'+type+'.bed.gz')))

if all(checklist):
    rule BamToBed:
        input:  "TRACKS/{combo}/{file}_mapped_{type}.bed.gz"
        output: "BED/{combo}/{file}_mapped_{type}.bed.gz"
        log:    "LOGS/Bed/linkbed{file}_{type}.log"
        conda:  "base.yaml"
        threads: 1
        params: abs = lambda wildcards: os.path.abspath('TRACKS/'+wildcards.file+'_mapped_'+wildcards.type+'.bed.gz')
        shell:  "ln -s {params.abs} {output}"
elif all(checklist2):
    rule BamToBed:
        input:  "PEAKS/{combo}/{file}_mapped_{type}.bed.gz"
        output: "BED/{combo}/{file}_mapped_{type}.bed.gz"
        log:    "LOGS/Bed/linkbed{file}_{type}.log"
        conda:  "base.yaml"
        threads: 1
        params: abs = lambda wildcards: os.path.abspath('PEAKS/'+wildcards.file+'_mapped_'+wildcards.type+'.bed.gz')
        shell:  "ln -s {params.abs} {output}"
else:
    if not stranded or stranded == 'fr':
        rule BamToBed:
            input:  "MAPPED/{combo}/{file}_mapped_sorted.bam",
                    "MAPPED/{combo}/{file}_mapped_sorted_unique.bam"
            output: "BED/{combo}/{file}_mapped_sorted.bed.gz",
                    "BED/{combo}/{file}_mapped_unique.bed.gz"
            log:    "LOGS/Bed/createbed{file}.log"
            conda:  "bedtools.yaml"
            threads: 1
            shell:  "bedtools bamtobed -split -i {input[0]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/2$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |gzip > {output[0]} 2> {log} && bedtools bamtobed -split -i {input[1]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/2$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |gzip > {output[1]} 2>> {log}"
    elif stranded and stranded == 'rf':
        rule BamToBed:
        rule BamToBed:
            input:  "MAPPED/{combo}/{file}_mapped_sorted.bam",
                    "MAPPED/{combo}/{file}_mapped_sorted_unique.bam"
            output: "BED/{combo}/{file}_mapped_sorted.bed.gz",
                    "BED/{combo}/{file}_mapped_unique.bed.gz"
            log:    "LOGS/Bed/createbed{file}.log"
            conda:  "bedtools.yaml"
            threads: 1
            shell:  "bedtools bamtobed -split -i {input[0]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/1$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |gzip > {output[0]} 2> {log} && bedtools bamtobed -split -i {input[1]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/1$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |gzip > {output[1]} 2>> {log}"

rule AnnotateBed:
    input:  rules.BamToBed.output
    output: "BED/{combo}/{file}_anno_{type}.bed.gz"
    log:    "LOGS/Bed/annobeds_{type}_{file}.log"
    conda:  "perl.yaml"
    threads: 1
    params: bins=BINS,
            anno = lambda wildcards: tool_params(wildcards.file, None, config, 'ANNOTATE').get('ANNOTATION', ""),
            annop = lambda wildcards: tool_params(wildcards.file, None, config, 'ANNOTATE')['OPTIONS'].get('ANNOTATE', ""),
            annof = lambda wildcards: tool_params(wildcards.file, None, config, 'ANNOTATE').get('ANNOFEATURE', "")
    shell:  "perl {params.bins}/Universal/AnnotateBed.pl -b {input[0]} -a {params.anno} {params.annof} {params.annop} |gzip > {output[0]}"

rule UnzipGenome:
    input:  ref = REFERENCE,
    output: fa = expand("{ref}.fa", ref=REFERENCE.replace('.fa.gz', '')),
            fai = expand("{ref}.fa.fai", ref=REFERENCE.replace('.fa.gz', '')),
            fas = expand("{ref}.chrom.sizes", ref=REFERENCE.replace('.fa.gz', ''))
    log:    expand("LOGS/PEAKS/{combo}/indexfa.log", combo=combo)
    conda:  "samtools.yaml"
    threads: 1
    params: bins = BINS
    shell:  "set +o pipefail; zcat {input[0]} |perl -F\\\\040 -wane 'if($_ =~ /^>/){{$F[0] = $F[0] =~ /^>chr/ ? $F[0] : \">chr\".substr($F[0],1);chomp($F[0]);print \"\\n\".$F[0].\"\\n\"}} else{{($line=$_)=~s/\\r[\\n]*/\\n/gm; chomp($line=$_); print $line}}' |tail -n+2 > {output.fa} && {params.bins}/Preprocessing/indexfa.sh {output.fa} 2> {log} && cut -f1,2 {output.fai} > {output.fas}"

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
    

rule AddSequenceToBed:
    input:  bd = rules.AnnotateBed.output,
            fa = expand("{ref}.fa", ref=REFERENCE.replace('.fa.gz', '')),
    output: bed = "BED/{combo}/{file}_anno_seq_{type}.bed.gz",
            bt = temp("BED/{combo}/{file}_bed_chr_{type}.tmp"),
            bs = temp("BED/{combo}/{file}_bed_seq_{type}.tmp")
    log:    "LOGS/BED/seq2bed{type}_{file}.log"
    conda:  "bedtools.yaml"
    threads: 1
    params: bins=BINS,
            sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
    shell:  "export LC_ALL=C; zcat {input.bd} | perl -wlane '$F[0] = $F[0] =~ /^chr/ ? $F[0] : \"chr\".$F[0]; print join(\"\\t\",@F[0..5])' > {output.bt} && bedtools getfasta -fi {input.fa} -bed {output.bt} -name+ -tab -s -fullHeader -fo {output.bs} && cut -d$'\t' -f2 {output.bs}|sed 's/t/u/ig'|paste -d$'\t' <(zcat {input.bd}) - |sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.bed}"

rule MergeAnnoBed:
    input:  rules.AddSequenceToBed.output
    output: "BED/{combo}/{file}_anno_seq_{type}_merged.bed.gz"
    log:    "LOGS/Bed/mergebeds_{type}_{file}.log"
    conda:  "bedtools.yaml"
    threads: 1
    params: sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
    shell:  "export LC_ALL=C; zcat {input[0]}|perl -wlane 'print join(\"\t\",@F[0..6],$F[-3],$F[-2])' |bedtools merge -s -c 7,8,9 -o distinct -delim \"|\" |sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k1,1 -k2,2n|gzip > {output[0]}"
