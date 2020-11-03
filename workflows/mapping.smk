rule sortsam:
    input:  mapps = rules.mapping.output.mapped
    output: sortedsam = report("MAPPED/{file}_mapped_sorted.sam.gz", category="SORTING"),
            tmphead = temp("MAPPED/{file}_mapped_header.gz"),
            tmpfile = temp("TMP/{file}")
    log:    "LOGS/{file}/sortsam.log"
    conda: "nextsnakes/envs/samtools.yaml"
    threads: MAXTHREAD
    params: linkto = lambda wildcards, output: os.path.basename(output.sortedsam)
    shell: "set +o pipefail;samtools view -H {input.mapps}|grep -P '^@HD' |pigz -p {threads} -f > {output.tmphead} ; samtools view -H {input.mapps}|grep -P '^@SQ'|sort -t$'\t' -k1,1 -k2,2V |pigz -p {threads} -f >> {output.tmphead} ; samtools view -H {input.mapps}|grep -P '^@RG'|pigz -p {threads} -f >> {output.tmphead} ; samtools view -H {input.mapps}|grep -P '^@PG'|pigz -p {threads} -f >> {output.tmphead} ; export LC_ALL=C;samtools view -h {input.mapps} | grep -v \"^@\"|sort --parallel={threads} -S 25% -T TMP -t$'\t' -k3,3V -k4,4n - |pigz -p {threads} -f > {output.tmpfile} ; cat {output.tmphead} {output.tmpfile} > {output.sortedsam} 2> {log}"# && rm -f {input.mapps} && touch {input.mapps}"

rule sam2bam:
    input:  sortedsam = rules.sortsam.output.sortedsam
    output: bam = report("MAPPED/{file}_mapped_sorted.bam", category="2BAM"),
            bamindex = "MAPPED/{file}_mapped_sorted.bam.bai"
    log:    "LOGS/{file}/sam2bam.log"
    conda: "nextsnakes/envs/samtools.yaml"
    threads: MAXTHREAD
    params: bins = BINS,
            fn = lambda wildcards: "{fn}".format(fn=str(sample_from_path(wildcards.file)))
    shell: "zcat {input.sortedsam} | samtools view -bS - | samtools sort -T {params.fn} -o {output.bam} --threads {threads} && samtools index {output.bam} 2> {log}"

rule uniqsam:
    input:  sortedsam = rules.sortsam.output.sortedsam,
            bam = rules.sam2bam.output
    output: uniqsam = report("MAPPED/{file}_mapped_sorted_unique.sam.gz", category="UNIQUE")
    log: "LOGS/{file}/uniqsam.log"
    conda: "nextsnakes/envs/base.yaml"
    threads: MAXTHREAD
    params: bins=BINS
    shell:  "{params.bins}/Shells/UniqueSam_woPicard.sh {input.sortedsam} {output.uniqsam} {threads} 2> {log}"

rule sam2bamuniq:
    input: uniqsam = rules.uniqsam.output,
           bam = rules.sam2bam.output
    output:  uniqbam = report("MAPPED/{file}_mapped_sorted_unique.bam", category="2BAM"),
             uniqbamindex = "MAPPED/{file}_mapped_sorted_unique.bam.bai"
    log:     "LOGS/{file}/sam2bamuniq.log"
    conda:   "nextsnakes/envs/samtools.yaml"
    threads: MAXTHREAD
    params: bins = BINS,
            fn = lambda wildcards: "{fn}".format(fn=sample_from_path(wildcards.file))
    shell: "zcat {input.uniqsam} | samtools view -bS - | samtools sort -T {params.fn} -o {output.uniqbam} --threads {threads} && samtools index {output.uniqbam} 2> {log}"

rule dedupbam:
    input:  bam = rules.sam2bam.output.bam
    output: bam = report("MAPPED/{file}_mapped_sorted_dedup.bam", category="DEDUP")
    log:    "LOGS/{file}/dedupbam.log"
    conda:  "nextsnakes/envs/"+DEDUPENV+".yaml"
    threads: MAXTHREAD
    params: dpara = lambda wildcards: ' '.join("{!s}={!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "DEDUP")['OPTIONS'][1].items()),
            dedup = DEDUPBIN
    shell: "{params.dedup} dedup {params.dpara} --stdin={input.bam} --log={log} --stdout={output.bam} 2>> {log}"

rule dedupuniqbam:
    input:  bam = rules.sam2bamuniq.output.uniqbam,
            check = rules.dedupbam.output.bam
    output: bam = report("MAPPED/{file}_mapped_sorted_unique_dedup.bam", category="DEDUP")
    log:    "LOGS/{file}/dedupuniqbam.log"
    conda:  "nextsnakes/envs/"+DEDUPENV+".yaml"
    threads: MAXTHREAD
    params: dpara = lambda wildcards: ' '.join("{!s}={!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "DEDUP")['OPTIONS'][1].items()),
            dedup = DEDUPBIN
    shell: "{params.dedup} dedup {params.dpara} --stdin={input.bam} --log={log} --stdout={output.bam} 2>> {log}"
