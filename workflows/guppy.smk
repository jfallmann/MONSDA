CALLERBIN, CALLERENV = env_bin_from_config(config,'BASECALL')

wildcard_constraints:
    rawfile = '|'.join(SAMPLES)

rule themall:
    input: expand("FASTQ/{rawfile}.fastq.gz", rawfile = SAMPLES)

rule call_base:
    input:  f5 = "RAW/{rawfile}.fast5"
    output: fq = "FASTQ/{rawfile}.fastq.gz",
            summary = "FASTQ/{rawfile}_summary.txt",
            telemetry = "FASTQ/{rawfile}_telemetry.js"
    log:    "LOGS/BASECALL/{rawfile}_guppy.log"
    conda:  ""+CALLERENV+".yaml"
    threads: MAXTHREAD
    params: caller = CALLERBIN,
            cpara = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'BASECALL', CALLERENV)['OPTIONS'].get('BASECALL', ""),
            cmodel = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'BASECALL', CALLERENV)['OPTIONS'].get('MODEL', ""),
            f5dir = lambda wildcards, input: os.path.dirname(input.f5),
            f5file = lambda wildcards, input: os.path.basename(input.f5),
            fqdir = lambda wildcards, output: os.path.dirname(output.fq)
    shell: " echo \"{params.f5file}\" > {params.f5dir}/f5list && {params.caller} {params.cpara} -c {params.cmodel} --compress_fastq -i {params.f5dir} --input_file_list {params.f5dir}/f5list -s {params.f5dir}/BASECALL 2> {log} && cat {params.f5dir}/BASECALL/pass/fastq_runid_*.fastq.gz > {output.fq} && rm -f {params.f5dir}/BASECALL/pass/fastq_runid_*.fastq.gz && cat {params.f5dir}/BASECALL/*.log >> {log} && rm -f {params.f5dir}/BASECALL/*.log && mv -f {params.f5dir}/BASECALL/sequencing_summary.txt {output.summary} &&  mv -f {params.f5dir}/BASECALL/sequencing_telemetry.js {output.telemetry} && rm -f {params.f5dir}/f5list"