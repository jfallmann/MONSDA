QCBIN, QCENV = env_bin_from_config2(SAMPLES,config,'QC')

if paired == 'paired':
    log.info('Running paired mode QC')
    rule qc_raw:
        input:  r1 = expand("FASTQ/{rawfile}_{read}.fastq.gz", rawfile=list(SAMPLES), read=['R1','R2'])
        output: o1 = report(expand("QC/FASTQC/{rawfile}_{read}_fastqc.zip", rawfile=list(SAMPLES), read=['R1','R2']), category="QC")
        log:    expand("LOGS/{rawfile}/fastqc_{read}_raw.log", rawfile=list(SAMPLES), read=['R1','R2'])
        conda:  "nextsnakes/envs/"+QCENV+".yaml"
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None ,config, 'QC')['OPTIONS'][0].items())
        shell: "arr=({input.r1}); orr=({output.o1});alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do OUT=$(dirname ${{orr[$i]}}); fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq  ${{arr[$i]}} 2> {log};done"

    rule qc_trimmed:
        input:  r1 = expand("TRIMMED_FASTQ/{file}_{read}_trimmed.fastq.gz", file=samplecond(SAMPLES,config), read=["R1", "R2"])
        output: o1 = report(expand("QC/FASTQC/{file}_{read}_trimmed_fastqc.zip", file=samplecond(SAMPLES,config), read=['R1','R2']), category="QC")
        log:    expand("LOGS/{file}/fastqc_{read}_trimmed.log", file=samplecond(SAMPLES,config), read=['R1','R2'])
        conda:  "nextsnakes/envs/"+QCENV+".yaml"
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None ,config, 'QC')['OPTIONS'][0].items())
        shell: "arr=({input.r1}); orr=({output.o1});alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do OUT=$(dirname ${{orr[$i]}}); fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq  ${{arr[$i]}} 2> {log};done"

    rule qc_mapped:
        input:  r1 = "MAPPED/{file}_mapped_sorted.sam.gz"
        output: o1 = report("QC/FASTQC/{file}_mapped_sorted_fastqc.zip", category="QC")
        log:    "LOGS/{file}/fastqc_mapped.log"
        conda:  "nextsnakes/envs/"+QCENV+".yaml"
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'QC')['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f sam_mapped {input.r1} 2> {log}"

    rule qc_uniquemapped:
        input:  r1 = "MAPPED/{file}_mapped_sorted_unique.bam",
                r2 = "MAPPED/{file}_mapped_sorted_unique.bam.bai"
        output: o1 = report("QC/FASTQC/{file}_mapped_sorted_unique_fastqc.zip", category="QC")
        log:    "LOGS/{file}/fastqc_uniquemapped.log"
        conda:  "nextsnakes/envs/"+QCENV+".yaml"
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'QC')['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f bam {input.r1} 2> {log}"

else:
    rule qc_raw:
        input:  r1 = expand("FASTQ/{rawfile}.fastq.gz", rawfile=list(SAMPLES))
        output: o1 = report(expand("QC/FASTQC/{rawfile}_fastqc.zip", rawfile=list(SAMPLES)), category="QC")
        log:    expand("LOGS/{rawfile}/fastqc_raw.log", rawfile=list(SAMPLES))
        conda:  "nextsnakes/envs/"+QCENV+".yaml"
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None ,config, 'QC')['OPTIONS'][0].items())
        shell: "arr=({input.r1}); orr=({output.o1});alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do OUT=$(dirname ${{orr[$i]}}); fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq  ${{arr[$i]}} 2>> {log};done"

    rule qc_trimmed:
        input:  r1 = expand("TRIMMED_FASTQ/{file}_trimmed.fastq.gz", file=samplecond(SAMPLES,config))
        output: o1 = report(expand("QC/FASTQC/{file}_trimmed_fastqc.zip", file=samplecond(SAMPLES,config)), category="QC")
        log:    expand("LOGS/{file}/fastqc_trimmed.log", file=samplecond(SAMPLES,config))
        conda:  "nextsnakes/envs/"+QCENV+".yaml"
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None ,config, 'QC')['OPTIONS'][0].items())
        shell: "arr=({input.r1}); orr=({output.o1});alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do OUT=$(dirname ${{orr[$i]}}); fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq  ${{arr[$i]}} 2>> {log};done"

    rule qc_mapped:
        input:   r1 = "MAPPED/{file}_mapped_sorted.sam.gz"
        output:  o1 = report("QC/FASTQC/{file}_mapped_sorted_fastqc.zip", category="QC")
        log:     "LOGS/{file}/fastqc_mapped.log"
        conda:  "nextsnakes/envs/"+QCENV+".yaml"
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'QC')['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f sam_mapped {input.r1} 2> {log}"

    rule qc_uniquemapped:
        input:  r1 = "MAPPED/{file}_mapped_sorted_unique.bam",
                r2 = "MAPPED/{file}_mapped_sorted_unique.bam.bai"
        output: o1 = report("QC/FASTQC/{file}_mapped_sorted_unique_fastqc.zip", category="QC")
        log:    "LOGS/{file}/fastqc_uniquemapped.log"
        conda:  "nextsnakes/envs/"+QCENV+".yaml"
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'QC')['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f bam {input.r1} 2> {log}"
