CALLERBIN, CALLERENV = env_bin_from_config2(SAMPLES,config,'BASECALL')

wildcard_constraints:
    rawfile = '|'.join(list(SAMPLES))

rule themall:
    input: expand("FASTQ/{rawfile}.fastq.gz", rawfile = SAMPLES)

rule call_base:
    input:  f5 = "RAW/{rawfile}.fast5"
    output: fq = "FASTQ/{rawfile}.fastq.gz",
            summary = "FASTQ/{rawfile}_summary.txt",
            telemetry = "FASTQ/{rawfile}_telemetry.js"
    log:    "LOGS/BASECALL/{rawfile}_guppy.log"
    conda:  "nextsnakes/envs/"+CALLERENV+".yaml"
    threads: int(MAXTHREAD/4) if int(MAXTHREAD/4) >= 1 else 1
    params: caller = CALLERBIN,
            cpara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None, config, 'BASECALL')['OPTIONS'][0].items()),
            f5dir = lambda wildcards, input: os.path.dirname(input.f5),
            f5file = lambda wildcards, input: os.path.basename(input.f5),
            fqdir = lambda wildcards, output: os.path.dirname(output.fq),
            callers = lambda wildcards, threads: int(MAXTHREAD/threads)
    shell: " echo \"{params.f5file}\" > f5list && {params.caller} {params.cpara} --cpu_threads_per_caller {threads} --num_callers {params.callers} --compress_fastq -i {params.f5dir} --input_file_list f5list -s {params.fqdir} 2> {log} && cat {params.fqdir}/fastq_runid_*.fastq.gz > {output.fq} && rm -f {params.fqdir}/fastq_runid_*.fastq.gz && cat {params.fqdir}/*.log >> {log} && rm -f {params.fqdir}/*.log && mv -f {params.fqdir}/sequencing_summary.txt {output.summary} &&  mv -f {params.fqdir}/sequencing_telemetry.js {output.telemetry} && rm -f f5list"
