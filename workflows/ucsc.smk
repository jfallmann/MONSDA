wildcard_constraints:
    type="sorted|unique"

rule all:
    input:  expand("DONE/UCSC/{file}_{type}_tracks",file=samplecond(SAMPLES,config),type=['sorted','unique'])

checklist = list()
checklist2 = list()
for file in samplecond(SAMPLES,config):
    for type in ['sorted','unique']:
        checklist.append(os.path.isfile(os.path.abspath('BED/'+file+'_mapped_'+type+'.bed.gz')) and not os.path.islink(os.path.abspath('BED/'+file+'_mapped_'+type+'.bed.gz')))
        checklist2.append(os.path.isfile(os.path.abspath('PEAKS/'+file+'_mapped_'+type+'.bed.gz')) and not os.path.islink(os.path.abspath('PEAKS/'+file+'_mapped_'+type+'.bed.gz')))

if all(checklist):
    rule BamToBed:
        input:  "BED/{file}_mapped_{type}.bed.gz"
        output: "UCSC/{file}_mapped_{type}.bed.gz"
        log:    "LOGS/UCSC/linkbed{file}_{type}.log"
        conda:  "snakes/envs/base.yaml"
        threads: 1
        params: abs = lambda wildcards: os.path.abspath('BED/'+wildcards.file+'_mapped_'+wildcards.type+'.bed.gz')
        shell:  "ln -s {params.abs} {output}"
elif all(checklist2):
    rule BamToBed:
        input:  "PEAKS/{file}_mapped_{type}.bed.gz"
        output: "UCSC/{file}_mapped_{type}.bed.gz"
        log:    "LOGS/UCSC/linkbed{file}_{type}.log"
        conda:  "snakes/envs/base.yaml"
        threads: 1
        params: abs = lambda wildcards: os.path.abspath('PEAKS/'+wildcards.file+'_mapped_'+wildcards.type+'.bed.gz')
        shell:  "ln -s {params.abs} {output}"
else:
    if not stranded or stranded == 'fr':
        rule BamToBed:
            input:  "SORTED_MAPPED/{file}_mapped_sorted.bam",
                    "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
            output: "UCSC/{file}_mapped_sorted.bed.gz",
                    "UCSC/{file}_mapped_unique.bed.gz"
            log:    "LOGS/UCSC/{file}_ucscbamtobed.log"
            conda:  "snakes/envs/bedtools.yaml"
            threads: 1
            shell:  "bedtools bamtobed -split -i {input[0]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/2$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |gzip > {output[0]} 2> {log} && bedtools bamtobed -split -i {input[1]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/2$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |gzip > {output[1]} 2>> {log}"
    elif stranded and stranded == 'rf':
        rule BamToBed:
            input:  "SORTED_MAPPED/{file}_mapped_sorted.bam",
                    "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
            output: "UCSC/{file}_mapped_sorted.bed.gz",
                    "UCSC/{file}_mapped_unique.bed.gz"
            log:    "LOGS/UCSC/{file}_ucscbamtobed.log"
            conda:  "snakes/envs/bedtools.yaml"
            threads: 1
            shell:  "bedtools bamtobed -split -i {input[0]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/1$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |gzip > {output[0]} 2> {log} && bedtools bamtobed -split -i {input[1]} |sed 's/ /\_/g'|perl -wl -a -F\'\\t\' -n -e '$F[0] =~ s/\s/_/g;if($F[3]=~/\/1$/){{if ($F[5] eq \"+\"){{$F[5] = \"-\"}}elsif($F[5] eq \"-\"){{$F[5] = \"+\"}}}} print join(\"\t\",@F[0..$#F])' |gzip > {output[1]} 2>> {log}"
        #        shell:  "bedtools bamtobed -split -i {input[0]} |gzip > {output[0]} && bedtools bamtobed -split -i {input[1]} |gzip > {output[1]}"
rule index_fa:
    input:  expand("{ref}/{{org}}/{{gen}}{{name}}.fa.gz",ref=REFERENCE),
    output: expand("{ref}/{{org}}/{{gen}}{{name}}.fa.fai",ref=REFERENCE)
    log:    "LOGS/UCSC/{org}/{gen}{name}_ucscindexfa.log"
    conda:  "snakes/envs/samtools.yaml"
    params: bins = BINS
    threads: 1
    shell:  "for i in {input};do {params.bins}/Preprocessing/indexfa.sh $i 2> {log};done"

rule get_chromsize_genomic:
    input: expand("{ref}/{{org}}/{{gen}}{{name}}.fa.fai",ref=REFERENCE)
    output: expand("{ref}/{{org}}/{{gen}}{{name}}.chrom.sizes",ref=REFERENCE)
    log:    "LOGS/UCSC/{org}/{gen}{name}_ucscgetchrom.log"
    conda:  "snakes/envs/samtools.yaml"
    threads: 1
    params: bins = BINS
    shell:  "cut -f1,2 {input} > {output} 2> {log}"

checklist = list()
checklist2 = list()
for file in samplecond(SAMPLES,config):
    for type in ['sorted','unique']:
        for orient in ['fw','rw']:
            checklist.append(os.path.isfile(os.path.abspath('BED/'+file+'_mapped_'+type+'.'+orient+'.bedg.gz')))
            checklist2.append(os.path.isfile(os.path.abspath('PEAKS/'+file+'_mapped_'+type+'.'+orient+'.bedg.gz')))

if all(checklist):
    rule BedToBedg:
        input:  "BED/{file}_mapped_{type}.{orient}.bedg.gz",
                sizes = lambda wildcards: "{ref}/{gen}{name}.chrom.sizes".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config))
        output: "UCSC/{file}_mapped_{type}.{orient}.bedg.gz"
        log:    "LOGS/UCSC/{file}_{type}_{orient}_ucscbedtobedgraph.log"
        conda:  "snakes/envs/ucsc.yaml"
        #    conda:  "snakes/envs/perl.yaml"
        threads: 1
        params: abs = lambda wildcards: os.path.abspath('BED/'+wildcards.file+'_mapped_'+wildcards.type+'.'+wildcards.orient+'.bedg.gz')
        shell:  "ln -s {params.abs} {output}"

elif all(checklist2):
    rule BedToBedg:
        input:  "PEAKS/{file}_mapped_{type}.{orient}.bedg.gz",
                sizes = lambda wildcards: "{ref}/{gen}{name}.chrom.sizes".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config))
        output: "UCSC/{file}_mapped_{type}.{orient}.bedg.gz"
        log:    "LOGS/UCSC/{file}_{type}_{orient}_ucscbedtobedgraph.log"
        conda:  "snakes/envs/ucsc.yaml"
        threads: 1
        params: abs = lambda wildcards: os.path.abspath('PEAKS/'+wildcards.file+'_mapped_'+wildcards.type+'.'+wildcards.orient+'.bedg.gz')
        shell:  "ln -s {params.abs} {output}"

else:
    rule BedToBedg:
        input:  bed = "UCSC/{file}_mapped_{type}.bed.gz",
                fai = lambda wildcards: "{ref}/{gen}{name}.fa.fai".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config)),
                sizes = lambda wildcards: "{ref}/{gen}{name}.chrom.sizes".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config))
        output: fw = "UCSC/{file}_mapped_{type}.fw.bedg.gz",
                re = "UCSC/{file}_mapped_{type}.re.bedg.gz"
        log:    "LOGS/UCSC/{file}_{type}_ucscbedtobedgraph.log"
        conda:  "snakes/envs/bedtools.yaml"
        #    conda:  "snakes/envs/perl.yaml"
        threads: 1
        shell: "export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i {input.bed} -bg -split -strand + -g {input.sizes} | sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.fw} 2> {log} && bedtools genomecov -i {input.bed} -bg -split -strand - -g {input.sizes} |sort --parallel={threads} -S 25% -T TMP -t$'\t' -k1,1 -k2,2n |gzip > {output.re} 2>> {log}"
        #        shell:  "awk '{{if($6==\"+\") print}}' {input[0]} | bedItemOverlapCount {params.genome} -chromSize={params.sizes} stdin |sort -k1,1 -k2,2n|gzip > {output[0]} 2> {log} && awk '{{if($6==\"-\") print}}' {input[0]} | bedItemOverlapCount {params.genome} -chromSize={params.sizes} stdin |sort -k1,1 -k2,2n|gzip > {output[1]} 2>> {log} && awk '{{if($6==\"+\") print}}' {input[1]} | bedItemOverlapCount {params.genome} -chromSize={params.sizes} stdin |sort -k1,1 -k2,2n|gzip > {output[2]} 2>> {log} && awk '{{if($6==\"-\") print}}' {input[1]} | bedItemOverlapCount {params.genome} -chromSize={params.sizes} stdin |sort -k1,1 -k2,2n|gzip > {output[3]} 2>> {log}"
        #    shell:  "perl {params.bins}/Universal/Bed2Bedgraph.pl -f {input[0]} -c {params.sizes} -v on -x {output[0]} -y {output[1]} -a track 2> {log} && perl {params.bins}/Universal/Bed2Bedgraph.pl -f {input[1]} -c {params.sizes} -v on -x {output[2]} -y {output[3]} -a track 2>> {log}"

### This step normalized the bedg files for comparison in the browser
rule NormalizeBedg:
    input:  fw = rules.BedToBedg.output.fw,
            re = rules.BedToBedg.output.re
    output: fw = "UCSC/{file}_mapped_{type}.fw.norm.bedg.gz",
            re = "UCSC/{file}_mapped_{type}.re.norm.bedg.gz"
    log:    "LOGS/UCSC/{file}_{type}_ucscnormalizebedgraph.log"
    conda:  "snakes/envs/perl.yaml"
    threads: 1
    shell: "export LC_ALL=C; if [[ -n \"$(zcat {input.fw} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then scale=$(bc <<< \"scale=6;1000000/$(zcat {input.fw}|cut -f4|sort -u|wc -l)\") perl -wlane '$sc=$ENV{{scale}};print join(\"\t\",@F[0..$#F-1]),\"\t\",$F[-1]/$sc' <(zcat {input.fw}) |gzip > {output.fw} 2> {log}; else gzip < /dev/null > {output.fw}; echo \"File {input.fw} empty\" >> {log}; fi && if [[ -n \"$(zcat {input.re} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then scale=$(bc <<< \"scale=6;1000000/$(zcat {input.re}|cut -f4|sort -u|wc -l)\") perl -wlane '$sc=$ENV{{scale}};print join(\"\t\",@F[0..$#F-1]),\"\t\",$F[-1]/$sc' <(zcat {input.re})|gzip > {output.re} 2> {log}; else gzip < /dev/null > {output.re}; echo \"File {input.re} empty\" >> {log}; fi"

### This step generates bigwig files for bedg which can then be copied to a web-browsable directory and uploaded to UCSC via the track field
rule BedgToUCSC:
    input:  fw = rules.NormalizeBedg.output.fw,
            re = rules.NormalizeBedg.output.re,
    output: fw = "UCSC/{file}_mapped_{type}.fw.bw",
            re = "UCSC/{file}_mapped_{type}.re.bw",
            t1 = temp("UCSC/{file}_mapped_{type}.fw.tmp"),
            t2 = temp("UCSC/{file}_mapped_{type}.re.tmp")
    log:    "LOGS/UCSC/{file}_{type}_bedgtoucsc.log"
    conda:  "snakes/envs/ucsc.yaml"
    threads: 1
    priority: 100               # This should be finished before we generate tracks
    params: sizes = lambda wildcards: "{ref}/{gen}{name}.chrom.sizes".format(ref=REFERENCE,gen=genomepath(wildcards.file,config),name=namefromfile(wildcards.file, config))
    shell:  "export LC_ALL=C; if [[ -n \"$(zcat {input.fw} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then zcat {input.fw} > {output.t1} && bedGraphToBigWig {output.t1} {params.sizes} {output.fw} 2> {log}; else touch {output.t1}; gzip < /dev/null > {output.fw}; echo \"File {input.fw} empty\" >> {log}; fi && if [[ -n \"$(zcat {input.re} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then zcat {input.re} > {output.t2} && bedGraphToBigWig {output.t2} {params.sizes} {output.re} 2>> {log}; else touch {output.t2}; gzip < /dev/null > {output.re}; echo \"File {input.re} empty\" >> {log}; fi"

rule GenerateTrack:
    input:  fw = rules.BedgToUCSC.output.fw,
            re = rules.BedgToUCSC.output.re
    output: "UCSC/{file}_mapped_{type}.fw.bw.trackdone",
            "UCSC/{file}_mapped_{type}.re.bw.trackdone"
    log:    "LOGS/UCSC/{file}_track_{type}.log"
    conda:  "snakes/envs/base.yaml"
    threads: MAXTHREAD
    params: bwdir = lambda wildcards: "UCSC/{src}".format(src=source_from_sample(wildcards.file,config)),
            bins = os.path.abspath(BINS),
            gen = lambda wildcards: os.path.basename(genomepath(wildcards.file,config)),
            options = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'UCSC')['OPTIONS'][0].items()),
            uid = lambda wildcards: "{src}".format(src='_'.join(source_from_sample(wildcards.file,config).split(os.sep)))
    shell: "echo -e \"{input.fw}\\n{input.re}\"|python3 {params.bins}/Analysis/GenerateTrackDb.py -i {params.uid} -e 1 -f STDIN -u '' -g {params.gen} {params.options} && touch {input.fw}\.trackdone && touch {input.re}.trackdone 2> {log}"

rule themall:
    input:  rules.GenerateTrack.output
    output: "DONE/UCSC/{file}_{type}_tracks"
    run:
        for f in output:
            with open(f, "w") as out:
                out.write("DONE")
onsuccess:
    print("Workflow finished, no error")
onerror:
	print("ERROR: "+str({log}))
