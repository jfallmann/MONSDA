include: "header.smk"
#include: "cutadapt.smk"

rule all:
    input: expand("DONE/{gen}_prepared", ref=config["REFERENCE"], gen=pathstogenomes(SAMPLES, config))

rule filter_mt:
    input:  expand("{ref}/{{org}}/{{gen}}.fa", ref=config["REFERENCE"])#, gen=pathstogenomes(SAMPLES, config))
    output: genomic = expand("{ref}/{{org}}/{{gen}}.genomic.fa", ref=config["REFERENCE"]),
            mt = expand("{ref}/{{org}}/{{gen}}.genomic_onlyMT.fa", ref=config["REFERENCE"]),#, gen=pathstogenomes(SAMPLES, config)),
            nomt = expand("{ref}/{{org}}/{{gen}}.genomic_withoutMT.fa", ref=config["REFERENCE"])#, gen=pathstogenomes(SAMPLES, config))
    log:    "LOGS/PrepareGenomes/{org}/{gen}/filter_mt.log"
    conda: "../envs/base.yaml"
    params: bins = BINS        # HIER WEITER
#            pos = lambda wildcards, input: input.index(wildards.input)
#            comp = i.replace('.genomic.fa','.genomic_onlyMT.fa'),
#            opos = output.mt.index(comp),
#            comp2 = i.replace('.genomic.fa','.genomic_withoutMT.fa'),
#            opos2 = output.nomt.index(comp2)
    threads: 1
    shell:  "{params.bins}/Preprocessing/transfer_multi.sh {input} {output.mt} {output.nomt} 2> {log}"

rule index_fa:
    input: expand("{ref}/{{org}}/{{gen}}.genomic_onlyMT.fa", ref=config["REFERENCE"]),# gen=pathstogenomes(SAMPLES, config)),
           expand("{ref}/{{org}}/{{gen}}.genomic_withoutMT.fa", ref=config["REFERENCE"]),# gen=pathstogenomes(SAMPLES, config)),
           expand("{ref}/{{org}}/{{gen}}.genomic.fa", ref=config["REFERENCE"])#, gen=pathstogenomes(SAMPLES, config))
    output: expand("{ref}/{{org}}/{{gen}}.genomic_onlyMT.fa.fai", ref=config["REFERENCE"]),# gen=pathstogenomes(SAMPLES, config)),
            expand("{ref}/{{org}}/{{gen}}.genomic_withoutMT.fa.fai", ref=config["REFERENCE"]),# gen=pathstogenomes(SAMPLES, config)),
            expand("{ref}/{{org}}/{{gen}}.genomic.fa.fai", ref=config["REFERENCE"])#, gen=pathstogenomes(SAMPLES, config))
    log:    "LOGS/{org}/{gen}/indexfa.log"
    conda:  "../envs/samtools.yaml"
    params: bins = BINS
    threads: 1
    shell:  "for i in {input};do {params.bins}/Preprocessing/indexfa.sh $i 2> {log};done"

rule get_chromsize_genomic:
    input:  expand("{ref}/{{org}}/{{gen}}.genomic.fa.fai", ref=config["REFERENCE"])
    output: expand("{ref}/{{org}}/{{gen}}.genomic_chrom.sizes", ref=config["REFERENCE"])
    log:    "LOGS/PrepareGenomes/{org}/{gen}/chromsize.log"
    conda:  "../envs/samtools.yaml"
    params: bins = BINS
    threads: 1
    shell:  "cut -f1,2 {input} > {output} 2> {log}"

rule scan_nuclear:
    input:  expand("{ref}/{{org}}/{{gen}}.genomic_withoutMT.fa", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.genomic_withoutMT.fa.fai", ref=config["REFERENCE"])
    output: expand("{ref}/{{org}}/{{gen}}.tRNAscan_withoutMT.fa", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.tRNAscan_withoutMT.bed", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.tRNAscan_withoutMT.out", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.tRNAscan_withoutMT.stat", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.tRNAscan_withoutMT.ss", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.tRNAscan_withoutMT.log", ref=config["REFERENCE"])
    log:    "LOGS/PrepareGenomes/{org}/{gen}/scan_nuclear.log"
    conda: "../envs/trnascan.yaml"
    params:  tscan = lambda w: trnascan_params(w.gen, "genomic", config)
    threads: 20
    shell:  "tRNAscan-SE {params.tscan} --thread {threads} -m {output[3]} -o {output[2]} -f {output[4]} -b {output[1]} -l {output[5]} -a {output[0]} {input[0]} 2> {log}"

rule scan_mitochondrial:
    input:  expand("{ref}/{{org}}/{{gen}}.genomic_onlyMT.fa", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.genomic_onlyMT.fa.fai", ref=config["REFERENCE"])
    output: expand("{ref}/{{org}}/{{gen}}.tRNAscan_onlyMT.fa", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.tRNAscan_onlyMT.bed", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.tRNAscan_onlyMT.out", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.tRNAscan_onlyMT.stat", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.tRNAscan_onlyMT.ss", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.tRNAscan_onlyMT.log", ref=config["REFERENCE"])
    log:    "LOGS/PrepareGenomes/{org}/{gen}/scan_mt.log"
    conda: "../envs/trnascan.yaml"
    params: gen = lambda w: w.gen+".genomic.fa",
            org = lambda w: w.org,
            ref = config["REFERENCE"],
            tscan = lambda w: trnascan_params(w.gen, "genomic", config)
    threads: 20
    shell:  "tRNAscan-SE {params.tscan} --thread {threads} -m {output[3]} -o {output[2]} -f {output[4]} -b {output[1]} -l {output[5]} -a {output[0]} {input[0]} 2> {log}"

rule remove_pseudogenes:
    input:  expand("{ref}/{{org}}/{{gen}}.tRNAscan_withoutMT.bed", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.tRNAscan_withoutMT.out", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.tRNAscan_onlyMT.bed", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.tRNAscan_onlyMT.out", ref=config["REFERENCE"])
    output: expand("{ref}/{{org}}/{{gen}}.tRNAscan_all_withoutPG.bed", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.tRNAscan_withoutMT_withoutPG.bed", ref=config["REFERENCE"])
    log:    "LOGS/PrepareGenomes/{org}/{gen}/rm_pseudo.log"
    conda: "../envs/perl.yaml"
    threads: 1
    params: bins = BINS,
            org = lambda w: w.org
    shell:  "perl {params.bins}/Mapping/1b_removePseudo.pl {input[1]} {input[0]} > {output[1]} 2> {log} && cat {input[2]} {output[1]} > {output[0]}"

rule check_distance:
    input:  expand("{ref}/{{org}}/{{gen}}.tRNAscan_all_withoutPG.bed", ref=config["REFERENCE"])
    output: expand("{ref}/{{org}}/{{gen}}.tRNAscan_all_withoutPG_sorted.bed", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.tRNAscan_all_withoutPG_sorted.dist", ref=config["REFERENCE"])
    log:    "LOGS/PrepareGenomes/{org}/{gen}/checkdist.log"
    conda: "../envs/perl.yaml"
    threads: 1
    params: bins = BINS
#            org = lambda w: w.org
    shell: "sort -k2 -n {input[0]} |uniq > {output[0]} && perl {params.bins}/Universal/gettRNAdist.pl {output[0]} {output[1]} 2> {log}"

rule mask_fasta:
    input:  expand("{ref}/{{org}}/{{gen}}.fa", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.tRNAscan_all_withoutPG_sorted.bed", ref=config["REFERENCE"])
    output: expand("{ref}/{{org}}/{{gen}}.genomic_masked.fa", ref=config["REFERENCE"])
    log:    "LOGS/PrepareGenomes/{org}/{gen}/maskfasta.log"
    conda: "../envs/bedtools.yaml"
    threads: 1
    params: bins = BINS
    shell: "bedtools maskfasta -fi {input[0]} -fo {output[0]} -bed {input[1]} -fullHeader 2> {log}"

rule create_pretrnas: ###create pre-tRNA library ##add 10 nt 5' and 3' flanking regions
    input:  expand("{ref}/{{org}}/{{gen}}.tRNAscan_all_withoutPG_sorted.bed", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.genomic_chrom.sizes", ref=config["REFERENCE"])
    output: expand("{ref}/{{org}}/{{gen}}.tRNAscan_pre-tRNAs.bed", ref=config["REFERENCE"])
    log:    "LOGS/PrepareGenomes/{org}/{gen}/create_pretrnas.log"
    conda: "../envs/perl.yaml"
    threads: 1
    params: bins = BINS,
            org = lambda w: w.org,
            flank = config["FLANK"]
#    shell: "perl {params.bins}/Universal/modBed12.pl {input[0]} {output[0]} 10"#Add extendbed.pl here!
    shell: "perl {params.bins}/Universal/ExtendBed.pl -g {input[1]} -b {input[0]} -l {params.flank} -r {params.flank} -o {output} -i ON 2> {log}"

rule rm_introns: ##remove introns, make fasta from bed
    input:  expand("{ref}/{{org}}/{{gen}}.tRNAscan_pre-tRNAs.bed", ref=config["REFERENCE"])
    output: expand("{ref}/{{org}}/{{gen}}.tRNAscan_pre-tRNAs.fa", ref=config["REFERENCE"])
    log:    "LOGS/PrepareGenomes/{org}/{gen}/rm_introns.log"
    conda: "../envs/bedtools.yaml"
    threads: 1
    params: gen = lambda w: w.gen+".genomic.fa",
            org = lambda w: w.org,
            ref = config["REFERENCE"]
    shell: "bedtools getfasta -name -split -s -fi {params.ref}/{params.org}/{params.gen} -bed {input[0]} -fo {output[0]} 2> {log}"

rule add_pretrnas:##add pre-tRNAs as extra chromosoms to the genome (get the artificial genome)
    input:  expand("{ref}/{{org}}/{{gen}}.genomic_masked.fa", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.tRNAscan_pre-tRNAs.fa", ref=config["REFERENCE"])
    output: expand("{ref}/{{org}}/{{gen}}.genomic_artificial.fa", ref=config["REFERENCE"])
    log:    "LOGS/PrepareGenomes/{org}/{gen}/add_pretrnas.log"
    conda: "../envs/base.yaml"
    threads: 1
    shell: "cat {input[0]} {input[1]} > {output[0]} 2> {log}"

rule index_artificial: ##indexing artificial genome
    input:  expand("{ref}/{{org}}/{{gen}}.genomic_artificial.fa", ref=config["REFERENCE"])
    output: expand("{ref}/{{org}}/{{gen}}.genomic_artificial.fa.fai", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.genomic_artificial.idx", ref=config["REFERENCE"])
    log:    "LOGS/PrepareGenomes/{org}/{gen}/index_artificial.log"
    conda: "../envs/index.yaml"
    threads: 20
    params: bins = BINS
    shell: "{params.bins}/Preprocessing/indexfa.sh {input} && segemehl.x -x {output[1]} -d {input} --threads {threads} 2> {log}"

rule create_mature: ###create mature tRNA library ##remove introns, make fasta from bed
    input:  expand("{ref}/{{org}}/{{gen}}.genomic.fa", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.tRNAscan_all_withoutPG.bed", ref=config["REFERENCE"])
    output: expand("{ref}/{{org}}/{{gen}}.tRNAscan_mature.fa", ref=config["REFERENCE"])
    log:    "LOGS/PrepareGenomes/{org}/{gen}/create_mature.log"
    conda: "../envs/bedtools.yaml"
    threads: 1
    shell: "bedtools getfasta -name -split -s -fi {input[0]} -bed {input[1]} -fo {output[0]} 2> {log}"

rule add_cca:##add CCA tail to tRNA chromosomes
    input:  expand("{ref}/{{org}}/{{gen}}.tRNAscan_mature.fa", ref=config["REFERENCE"])
    output: expand("{ref}/{{org}}/{{gen}}.tRNAscan_mature_CCA.fa", ref=config["REFERENCE"])
    log:    "LOGS/PrepareGenomes/{org}/{gen}/add_cca.log"
    conda: "../envs/perl.yaml"
    threads: 1
    params: bins = BINS
    shell: "perl {params.bins}/Universal/addCCA.pl {input} {output} 2> {log}"

rule cluster_mature:###mature tRNA clustering ##only identical tRNAs were clustered
    input:  expand("{ref}/{{org}}/{{gen}}.tRNAscan_mature_CCA.fa", ref=config["REFERENCE"])
    output: expand("{ref}/{{org}}/{{gen}}.tRNAscan_cluster.fa", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.tRNAscan_clusterInfo.fa", ref=config["REFERENCE"])
    log:    "LOGS/PrepareGenomes/{org}/{gen}/cluster_mature.log"
    conda: "../envs/perl.yaml"
    threads: 1
    params: bins = BINS
    shell: "perl {params.bins}/Universal/clustering.pl {input} {output[0]} {output[1]} 2> {log}"

rule index_cluster: ##indexing cluster
    input:  expand("{ref}/{{org}}/{{gen}}.tRNAscan_cluster.fa", ref=config["REFERENCE"])
    output: expand("{ref}/{{org}}/{{gen}}.tRNAscan_cluster.fa.fai", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.tRNAscan_cluster.idx", ref=config["REFERENCE"])
    log:    "LOGS/PrepareGenomes/{org}/{gen}/index_cluster.log"
    conda: "../envs/index.yaml"
    threads: 20
    params: bins = BINS
    shell: "{params.bins}/Preprocessing/indexfa.sh {input} && segemehl.x -x {output[1]} -d {input} --threads {threads} 2> {log}"

rule themall:
    input:  expand("{ref}/{{org}}/{{gen}}.tRNAscan_cluster.fa", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.tRNAscan_cluster.idx", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.tRNAscan_mature_CCA.fa", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.tRNAscan_mature.fa", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.genomic_artificial.idx", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.tRNAscan_all_withoutPG_sorted.dist", ref=config["REFERENCE"]),
            expand("{ref}/{{org}}/{{gen}}.genomic_masked.fa", ref=config["REFERENCE"])
    output: "DONE/{org}/{gen}_prepared"
    run:
        for f in output:
            with open(f, "w") as out:
                out.write("DONE")
