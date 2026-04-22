DEDUPBIN, DEDUPENV = env_bin_from_config(config, 'DEDUP')

wildcard_constraints:
    type = "sorted|sorted_unique"

if paired == 'paired':
    rule extract:
        input:  r1 = lambda wildcards: "FASTQ/{rawfile}_R1.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
                r2 = lambda wildcards: "FASTQ/{rawfile}_R2.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0])
        output: o1 = "DEDUP_FASTQ/{combo}/{file}_R1_dedup.fastq.gz",
                o2 = "DEDUP_FASTQ/{combo}/{file}_R2_dedup.fastq.gz",
                ubam = "TMP/FGEX/{combo}/{file}_extracted.bam",
                td = temp(directory("TMP/FGEX/{combo}/{file}"))
        log:   "LOGS/{combo}/{file}_dedup_extract.log"
        conda: ""+DEDUPENV+".yaml"
        container: "oras://jfallmann/monsda:"+DEDUPENV+""
        threads: 1
        params: epara = lambda wildcards: tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'].get('EXTRACT', ""),
                dedup = DEDUPBIN,
                sname = lambda wildcards: os.path.basename(wildcards.file)
        shell:  "mkdir -p {output.td} && {params.dedup} extract {params.epara} --inputs {input.r1} {input.r2} --sample {params.sname} --library {params.sname} --output {output.ubam} > {log} 2>&1 && samtools fastq -n -1 {output.o1} -2 {output.o2} -0 /dev/null -s /dev/null {output.ubam} >> {log} 2>&1"
else:
    rule extract:
        input:  r1 = lambda wildcards: "FASTQ/{rawfile}.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0])
        output: o1 = "DEDUP_FASTQ/{combo}/{file}_dedup.fastq.gz",
                ubam = "TMP/FGEX/{combo}/{file}_extracted.bam",
                td = temp(directory("TMP/FGEX/{combo}/{file}"))
        log:   "LOGS/{combo}/{file}_dedup_extract.log"
        conda: ""+DEDUPENV+".yaml"
        container: "oras://jfallmann/monsda:"+DEDUPENV+""
        threads: 1
        params: epara = lambda wildcards: tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'].get('EXTRACT', ""),
                dedup = DEDUPBIN,
                sname = lambda wildcards: os.path.basename(wildcards.file)
        shell:  "mkdir -p {output.td} && {params.dedup} extract {params.epara} --inputs {input.r1} --sample {params.sname} --library {params.sname} --output {output.ubam} > {log} 2>&1 && samtools fastq -n {output.ubam} | gzip -c > {output.o1} && echo done >> {log}"

if paired == 'paired':
    rule dedupbam:
        input:  bam = "MAPPED/{combo}/{file}_mapped_{type}.bam",
                ubam = "TMP/FGEX/{combo}/{file}_extracted.bam"
        output: bam = report("MAPPED/{combo}/{file}_mapped_{type}_dedup.bam", category="DEDUP"),
                bai = report("MAPPED/{combo}/{file}_mapped_{type}_dedup.bam.bai", category="DEDUP"),
                td = temp(directory("TMP/UMIDD/{combo}/{file}_{type}"))
        log:    "LOGS/{combo}/{file}_{type}/dedupbam.log"
        conda:  ""+DEDUPENV+".yaml"
        container: "oras://jfallmann/monsda:"+DEDUPENV+""
        threads: 1
        priority: 0               # This should be done after all mapping is done
        params: dpara = lambda wildcards: tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'].get('DEDUP', ""),
                dedup = DEDUPBIN
        shell: """mkdir -p {output.td}
{params.dedup} zipper --unmapped {input.ubam} --aligned {input.bam} --output {output.td}/zippered.bam > {log} 2>&1
{params.dedup} sort --order template-coordinate --input {output.td}/zippered.bam --output {output.td}/sorted.bam >> {log} 2>&1
{params.dedup} dedup {params.dpara} --input {output.td}/sorted.bam --output {output.bam} >> {log} 2>&1
samtools index {output.bam} >> {log} 2>&1
rm {input.ubam}"""
else:
    rule dedupbam:
        input:  bam = "MAPPED/{combo}/{file}_mapped_{type}.bam",
                ubam = "TMP/FGEX/{combo}/{file}_extracted.bam"
        output: bam = report("MAPPED/{combo}/{file}_mapped_{type}_dedup.bam", category="DEDUP"),
                bai = report("MAPPED/{combo}/{file}_mapped_{type}_dedup.bam.bai", category="DEDUP"),
                td = temp(directory("TMP/UMIDD/{combo}/{file}_{type}"))
        log:    "LOGS/{combo}/{file}_{type}/dedupbam.log"
        conda:  ""+DEDUPENV+".yaml"
        container: "oras://jfallmann/monsda:"+DEDUPENV+""
        threads: 1
        priority: 0               # This should be done after all mapping is done
        params: dpara = lambda wildcards: tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'].get('DEDUP', ""),
                dedup = DEDUPBIN
        shell: """mkdir -p {output.td}
{params.dedup} zipper --unmapped {input.ubam} --aligned {input.bam} --output {output.td}/zippered.bam > {log} 2>&1
{params.dedup} sort --order template-coordinate --input {output.td}/zippered.bam --output {output.td}/sorted.bam >> {log} 2>&1
{params.dedup} dedup {params.dpara} --input {output.td}/sorted.bam --output {output.bam} >> {log} 2>&1
samtools index {output.bam} >> {log} 2>&1
rm {input.ubam}"""
