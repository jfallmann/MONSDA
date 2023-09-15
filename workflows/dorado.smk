CALLERBIN, CALLERENV = env_bin_from_config(config,'BASECALL')

wildcard_constraints:
    rawfile = '|'.join(SAMPLES)

rule themall:
    input: expand("FASTQ/{rawfile}.fastq.gz", rawfile = SAMPLES)

rule call_base:
    input:  p5 = "RAW/{rawfile}.pod5"
    output: fq = "FASTQ/{rawfile}.fastq.gz",
            bam = temp("FASTQ/{rawfile}.bam")
    log:    "LOGS/BASECALL/{rawfile}_dorado.log"
    conda:  ""+CALLERENV+".yaml"
    threads: MAXTHREAD
    params: caller = CALLERBIN,
            cpara = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'BASECALL', CALLERENV)['OPTIONS'].get('BASECALL', ""),
            cmodel = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'BASECALL', CALLERENV)['OPTIONS'].get('MODEL', ""),
            p5dir = lambda wildcards, input: os.path.dirname(input.p5),
            p5file = lambda wildcards, input: os.path.basename(input.p5),
            fqdir = lambda wildcards, output: os.path.dirname(output.fq)
    shell: "{params.caller} download --directory {params.p5dir} --model {params.cmodel} &> {log} && {params.caller} basecaller {params.cpara} {params.p5dir}/{params.cmodel} {params.p5dir}/ 2>> {log} 1> {output.bam} && samtools view -h {output.bam}|samtools fastq -n - | pigz > {output.fq}"