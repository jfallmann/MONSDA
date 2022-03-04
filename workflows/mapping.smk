rule sortsam:
    input:  mapps = rules.mapping.output.mapped
    output: sortedsam = report("MAPPED/{combo}/{file}_mapped_sorted.sam.gz", category="SORTING"),
            tmphead = temp("MAPPED/{combo}/{file}_mapped_header.gz"),
            tmpfile = temp("TMP/{combo}/{file}")
    log:    "LOGS/{combo}/{file}/sortsam.log"
    conda: "samtools.yaml"
    threads: MAXTHREAD
    priority: 100
    params: linkto = lambda wildcards, output: os.path.basename(output.sortedsam),
            sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
    shell: "set +o pipefail;samtools view -H {input.mapps}|grep -P '^@HD' |pigz -p {threads} -f > {output.tmphead} ; samtools view -H {input.mapps}|grep -P '^@SQ'|sort -t$'\t' -k1,1 -k2,2V |pigz -p {threads} -f >> {output.tmphead} ; samtools view -H {input.mapps}|grep -P '^@RG'|pigz -p {threads} -f >> {output.tmphead} ; samtools view -H {input.mapps}|grep -P '^@PG'|pigz -p {threads} -f >> {output.tmphead} ; export LC_ALL=C;samtools view -h {input.mapps} | grep -v \"^@\"|sort --parallel={threads} -S {params.sortmem}% -T TMP -t$'\t' -k3,3V -k4,4n - |pigz -p {threads} -f > {output.tmpfile} ; cat {output.tmphead} {output.tmpfile} > {output.sortedsam} 2> {log}"# && rm -f {input.mapps} && touch {input.mapps}"

rule sam2bam:
    input:  sortedsam = rules.sortsam.output.sortedsam
    output: bam = report("MAPPED/{combo}/{file}_mapped_sorted.bam", category="2BAM"),
            bamindex = "MAPPED/{combo}/{file}_mapped_sorted.bam.bai"
    log:    "LOGS/{combo}/{file}/sam2bam.log"
    conda: "samtools.yaml"
    threads: MAXTHREAD
    params: bins = BINS
    shell: "zcat {input.sortedsam} | samtools view -bS - > {output.bam} && samtools index {output.bam} 2> {log}"

rule uniqsam:
    input:  sortedsam = rules.sortsam.output.sortedsam,
            bam = rules.sam2bam.output
    output: uniqsam = report("MAPPED/{combo}/{file}_mapped_sorted_unique.sam.gz", category="UNIQUE")
    log: "LOGS/{combo}/{file}/uniqsam.log"
    conda: "base.yaml"
    threads: MAXTHREAD
    params: bins=BINS
    shell:  "{params.bins}/Shells/UniqueSam_woPicard.sh {input.sortedsam} {output.uniqsam} {threads} 2> {log}"

rule sam2bamuniq:
    input: uniqsam = rules.uniqsam.output,
           bam = rules.sam2bam.output
    output:  uniqbam = report("MAPPED/{combo}/{file}_mapped_sorted_unique.bam", category="2BAM"),
             uniqbamindex = "MAPPED/{combo}/{file}_mapped_sorted_unique.bam.bai"
    log:     "LOGS/{combo}/{file}/sam2bamuniq.log"
    conda:   "samtools.yaml"
    threads: MAXTHREAD
    priority: 50
    params: bins = BINS
    shell: "zcat {input.uniqsam} | samtools view -bS - > {output.uniqbam} && samtools index {output.uniqbam} 2> {log}"
