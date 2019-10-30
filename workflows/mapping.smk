rule mapall:
    input: expand("DONE/{file}_mapped",file=samplecond(SAMPLES,config))

rule gzipsam:
    input:  mapps = rules.mapping.output.mapped
    output: gzipped = report("MAPPED/{file}_mapped.sam.gz", category="ZIPIT")
    log:    "LOGS/{file}/gzipsam.log"
    conda:  "snakes/envs/base.yaml"
    threads: MAXTHREAD
    shell: "pigz -k -p {threads} -f {input.mapps} > {output.gzipped} 2> {log}"

rule sortsam:
    input:  gzipped = rules.gzipsam.output
    output: sortedsam = report("SORTED_MAPPED/{file}_mapped_sorted.sam.gz", category="SORTING"),
            tmphead = temp("SORTED_MAPPED/{file}_mapped_header.gz"),
            tmpfile = temp("SORTTMP/{file}")
    log:    "LOGS/{file}/sortsam.log"
    conda: "snakes/envs/samtools.yaml"
    threads: MAXTHREAD
    shell: "set +o pipefail;samtools view -H {input.gzipped}|grep -P '^@HD' |pigz -p {threads} -f > {output.tmphead} ; samtools view -H {input.gzipped}|grep -P '^@SQ'|sort -t$'\t' -k1,1 -k2,2V |pigz -p {threads} -f >> {output.tmphead} ; samtools view -H {input.gzipped}|grep -P '^@RG'|pigz -p {threads} -f >> {output.tmphead} ; samtools view -H {input.gzipped}|grep -P '^@PG'|pigz -p {threads} -f >> {output.tmphead} ; export LC_ALL=C;zcat {input.gzipped} | grep -v \"^@\"|sort --parallel={threads} -S 25% -T SORTTMP -t$'\t' -k3,3V -k4,4n - |pigz -p {threads} -f > {output.tmpfile} ; cat {output.tmphead} {output.tmpfile} > {output.sortedsam} 2> {log}"

rule sam2bam:
    input:  sortedsam = rules.sortsam.output.sortedsam
    output: bam = report("SORTED_MAPPED/{file}_mapped_sorted.bam", category="2BAM"),
            bamindex = "SORTED_MAPPED/{file}_mapped_sorted.bam.bai"
    log:    "LOGS/{file}/sam2bam.log"
    conda: "snakes/envs/samtools.yaml"
    threads: MAXTHREAD
    params: bins = BINS,
            fn = lambda wildcards: "{fn}".format(fn=str(sample_from_path(wildcards.file)))
    shell: "zcat {input.sortedsam} | samtools view -bS - | samtools sort -T {params.fn} -o {output.bam} --threads {threads} && samtools index {output.bam} 2> {log}"

rule uniqsam:
    input: sortedsam = rules.sortsam.output.sortedsam,
           bam = rules.sam2bam.output
    output: uniqsam = report("UNIQUE_MAPPED/{file}_mapped_sorted_unique.sam.gz", category="UNIQUE")
    log: "LOGS/{file}/uniqsam.log"
    conda: "snakes/envs/base.yaml"
    threads: MAXTHREAD
    params: bins=BINS
    shell:  "{params.bins}/Shells/UniqueSam_woPicard.sh {input.sortedsam} {output.uniqsam} {threads} 2> {log}"

rule sam2bamuniq:
    input: uniqsam = rules.uniqsam.output,
           bam = rules.sam2bam.output
    output:  uniqbam = report("UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam", category="2BAM"),
             uniqbamindex = "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam.bai"
    log:     "LOGS/{file}/sam2bamuniq.log"
    conda:   "snakes/envs/samtools.yaml"
    threads: MAXTHREAD
    params: bins=BINS,
            fn = lambda wildcards: "{fn}".format(fn=sample_from_path(wildcards.file))
    shell: "zcat {input.uniqsam} | samtools view -bS - | samtools sort -T {params.fn} -o {output.uniqbam} --threads {threads} && samtools index {output.uniqbam} 2> {log}"

rule themall_mapping:
    input:  rules.sam2bamuniq.output
    output: "DONE/{file}_mapped"
    run:
        for f in output:
            with open(f, "w") as out:
                out.write("DONE")

onsuccess:
    print("Workflow finished, no error")
