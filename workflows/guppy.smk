CALLERBIN, CALLERENV = env_bin_from_config2(SAMPLES,config,'BASECALL')

wildcard_constraints:
    rawfile = '|'.join(list(SAMPLES)),

rule call_base:
    input:  f5 = "RAW/{rawfile}.fast5"
    output: fq = "FASTQ/{rawfile}.fastq.gz"
    log:    expand("LOGS/{sets}/{cape}.call.log", sets=SETS, cape=CALLERENV)
    conda:  "nextsnakes/envs/"+CALLERENV+".yaml"
    threads: int(MAXTHREAD/2)
    params: caller = CALLERBIN,
            cpara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None, config, 'CALLING')['OPTIONS'][0].items()),
            f5dir = lambda wildcards, input: os.path.dirname(input.f5),
            fqdir = lambda wildcards, input: os.path.dirname(output.fq)
    shell: "{params.caller} {params.cpara} --recursive --cpu_threads_per_caller {threads} --num_callers {threads} --verbose_logs --compress_fastq -i {params.f5dir} -s {params.fqdir} 2> {log}"
