rule sortsam:
    input:  mapps = rules.mapping.output.mapped
    output: sortedsam = report("MAPPED/{combo}/{file}_mapped_sorted.sam.gz", category="SORTING"),
            tmphead = temp("MAPPED/{combo}/{file}_mapped_header.gz"),
            tmpfile = temp("TMP/{combo}/{file}")
    log:    "LOGS/{combo}/{file}/sortsam.log"
    conda: "samtools.yaml"
    container: "oras://jfallmann/monsda:samtools"
    threads: MAXTHREAD
    priority: 100
    params: linkto = lambda wildcards, output: os.path.basename(output.sortedsam),
            sortmem = lambda w, resources: int(int(resources.mem_mb) / 1024)
    shell: "set +o pipefail; export LC_ALL=C; samtools view -H {input.mapps}|grep '^@HD' |pigz -p 1 -f > {output.tmphead} 2> {log}; samtools view -H {input.mapps}|grep '^@SQ'|sort -t$'\t' -k1,1 -k2,2V |pigz -p 1 -f >> {output.tmphead} 2>> {log}; samtools view -H {input.mapps}|grep '^@RG'|pigz -p 1 -f >> {output.tmphead} 2>> {log}; samtools view -H {input.mapps}|grep -P '^@PG'|pigz -p {threads} -f >> {output.tmphead} 2>> {log}; samtools view -h {input.mapps} | grep -v \"^@\"|sort --parallel={threads} -S {params.sortmem}G -T TMP -t$'\t' -k3,3V -k4,4n - |pigz -p {threads} -f > {output.tmpfile} 2>> {log}; cat {output.tmphead} {output.tmpfile} > {output.sortedsam} 2>> {log}"

rule sam2bam:
    input:  sortedsam = rules.sortsam.output.sortedsam
    output: bam = report("MAPPED/{combo}/{file}_mapped_sorted.bam", category="2BAM"),
            bamindex = "MAPPED/{combo}/{file}_mapped_sorted.bam.bai"
    log:    "LOGS/{combo}/{file}/sam2bam.log"
    conda: "samtools.yaml"
    container: "oras://jfallmann/monsda:samtools"
    threads: MAXTHREAD
    params: bins = BINS
    shell: "zcat {input.sortedsam} | samtools view -bS - > {output.bam} && samtools index {output.bam} 2> {log}"

rule uniqsam:
    input:  sortedsam = rules.sortsam.output.sortedsam,
            bam = rules.sam2bam.output
    output: uniqsam = report("MAPPED/{combo}/{file}_mapped_sorted_unique.sam.gz", category="UNIQUE")
    log: "LOGS/{combo}/{file}/uniqsam.log"
    conda: "samtools.yaml"
    container: "oras://jfallmann/monsda:samtools"
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
    container: "oras://jfallmann/monsda:samtools"
    threads: MAXTHREAD
    priority: 50
    params: bins = BINS
    shell: "zcat {input.uniqsam} | samtools view -bS - > {output.uniqbam} && samtools index {output.uniqbam} 2> {log}"
