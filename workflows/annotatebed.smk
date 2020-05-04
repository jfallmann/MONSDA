rule themall:
    input: expand("BED/{file}_anno_seq_{type}_merged.bed.gz",file=samplecond(SAMPLES,config), type=['sorted','unique'])
           #expand("DONE/BED/{file}_{type}",file=samplecond(SAMPLES,config), type=['sorted','unique'])

checklist = list()
checklist2 = list()
for file in samplecond(SAMPLES,config):
    for type in ['sorted','unique']:
        checklist.append(os.path.isfile(os.path.abspath('UCSC/'+file+'_mapped_'+type+'.bed.gz')) and not os.path.islink(os.path.abspath('UCSC/'+file+'_mapped_'+type+'.bed.gz')))
        checklist2.append(os.path.isfile(os.path.abspath('PEAKS/'+file+'_mapped_'+type+'.bed.gz')) and not os.path.islink(os.path.abspath('PEAKS/'+file+'_mapped_'+type+'.bed.gz')))

if all(checklist):
    rule BamToBed:
        input:  "UCSC/{file}_mapped_{type}.bed.gz"
        output: "BED/{file}_mapped_{type}.bed.gz"
        log:    "LOGS/Bed/linkbed{file}_{type}.log"
        conda:  "snakes/envs/base.yaml"
        threads: 1
        params: abs = lambda wildcards: os.path.abspath('UCSC/'+wildcards.file+'_mapped_'+wildcards.type+'.bed.gz')
        shell:  "ln -s {params.abs} {output}"
elif all(checklist2):
    rule BamToBed:
        input:  "PEAKS/{file}_mapped_{type}.bed.gz"
        output: "BED/{file}_mapped_{type}.bed.gz"
        log:    "LOGS/Bed/linkbed{file}_{type}.log"
        conda:  "snakes/envs/base.yaml"
        threads: 1
        params: abs = lambda wildcards: os.path.abspath('PEAKS/'+wildcards.file+'_mapped_'+wildcards.type+'.bed.gz')
        shell:  "ln -s {params.abs} {output}"
else:
    if not stranded or stranded == 'fr':
        rule BamToBed:
            input:  "MAPPED/{file}_mapped_sorted.bam",
                    "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
            output: "BED/{file}_mapped_sorted.bed.gz",
                    "BED/{file}_mapped_unique.bed.gz"
            log:    "LOGS/Bed/createbed{file}.log"
            conda:  "snakes/envs/bedtools.yaml"
            threads: 1
            shell:  "bedtools bamtobed -split -i {input[0]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/2$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |gzip > {output[0]} 2> {log} && bedtools bamtobed -split -i {input[1]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/2$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |gzip > {output[1]} 2>> {log}"
    elif stranded and stranded == 'rf':
        rule BamToBed:
        rule BamToBed:
            input:  "MAPPED/{file}_mapped_sorted.bam",
                    "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
            output: "BED/{file}_mapped_sorted.bed.gz",
                    "BED/{file}_mapped_unique.bed.gz"
            log:    "LOGS/Bed/createbed{file}.log"
            conda:  "snakes/envs/bedtools.yaml"
            threads: 1
            shell:  "bedtools bamtobed -split -i {input[0]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/1$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |gzip > {output[0]} 2> {log} && bedtools bamtobed -split -i {input[1]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/1$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |gzip > {output[1]} 2>> {log}"

rule AnnotateBed:
    input:  rules.BamToBed.output
    output: "BED/{file}_anno_{type}.bed.gz"
    log:    "LOGS/Bed/annobeds_{type}_{file}.log"
    conda:  "snakes/envs/perl.yaml"
    threads: 1
    params: bins=BINS,
            anno = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'ANNOTATE')['ANNOTATION']]),
            annop = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'ANNOTATE')['OPTIONS'][0].items()),
            annof = lambda wildcards: tool_params(wildcards.file, None, config, 'ANNOTATE')['ANNOFEATURE']
    shell:  "perl {params.bins}/Universal/AnnotateBed.pl -b {input[0]} -a {params.anno} {params.annof} {params.annop} |gzip > {output[0]}"

rule UnzipGenome:
    input:  expand("{ref}/{{org}}/{{gen}}{{name}}.fa.gz",ref=REFERENCE),
    output: fa = expand("{ref}/{{org}}/{{gen}}{{name}}_fastafrombed.fa",ref=REFERENCE)
    log:    "LOGS/Peaks/{org}/{gen}{name}/indexfa.log"
    conda:  "snakes/envs/samtools.yaml"
    threads: 1
    params: bins = BINS
    shell:  "zcat {input[0]} |perl -F\\\\040 -wlane 'if($_ =~ /^>/){{($F[0] = $F[0] =~ /^>chr/ ? $F[0] : \">chr\".substr($F[0],1))=~ s/\_/\./g;print $F[0]}}else{{print}}' > {output.fa} && {params.bins}/Preprocessing/indexfa.sh {output.fa} 2> {log}"

rule AddSequenceToBed:
    input:  bd = rules.AnnotateBed.output,
            fa = lambda wildcards: "{ref}/{gen}{name}_fastafrombed.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config))
    output: bed = "BED/{file}_anno_seq_{type}.bed.gz",
            bt = temp("BED/{file}_bed_chr_{type}.tmp"),
            bs = temp("BED/{file}_bed_seq_{type}.tmp")
    log:    "LOGS/BED/seq2bed{type}_{file}.log"
    conda:  "snakes/envs/bedtools.yaml"
    threads: 1
    params: bins=BINS
    shell:  "export LC_ALL=C; zcat {input.bd} | perl -wlane '$F[0] = $F[0] =~ /^chr/ ? $F[0] : \"chr\".$F[0]; print join(\"\\t\",@F[0..5])' > {output.bt} && fastaFromBed -fi {input.fa} -bed {output.bt} -name+ -tab -s -fullHeader -fo {output.bs} && cut -d$'\t' -f2 {output.bs}|sed 's/t/u/ig'|paste -d$'\t' <(zcat {input.bd}) - |sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.bed}"  # NEED TO GET RID OF SPACES AND WHATEVER IN HEADER

#rule AddSequenceToBed:
#    input:  rules.AnnotateBed.output
#    output: "BED/{file}_anno_seq_{type}.bed.gz",
#            temp("BED/{file}_anno_seq_{type}.tmp")
#    log:    "LOGS/Bed/seq2bed_{type}_{file}.log"
#    conda:  "snakes/envs/bedtools.yaml"
#    threads: 1
#    params: fasta = lambda wildcards: "{ref}/{gen}{name}.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config)),
#    shell:  "export LC_ALL=C; if [ ! -f \"{params.fasta}\" ] && [ -f \"{params.fasta}.gz\" ];then zcat {params.fasta}.gz > {params.fasta};fi && fastaFromBed -fi {params.fasta} -bed <(zcat {input[0]}|cut -d$'\t' -f 1-6) -name+ -tab -s -fullHeader -fo {output[1]} && cut -d$'\t' -f2 {output[1]}|sed 's/t/u/ig'|paste -d$'\t' <(zcat {input[0]}) -|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip  > {output[0]}"

rule MergeAnnoBed:
    input:  rules.AddSequenceToBed.output
    output: "BED/{file}_anno_seq_{type}_merged.bed.gz"
    log:    "LOGS/Bed/mergebeds_{type}_{file}.log"
    conda:  "snakes/envs/bedtools.yaml"
    threads: 1
    shell:  "export LC_ALL=C; zcat {input[0]}|perl -wlane 'print join(\"\t\",@F[0..6],$F[-3],$F[-2])' |bedtools merge -s -c 7,8,9 -o distinct -delim \"|\" |sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n|gzip > {output[0]}"
