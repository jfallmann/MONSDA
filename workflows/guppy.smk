CALLERBIN, CALLERENV = env_bin_from_config2(SAMPLES,config,'BASECALL')

rule call_base:
    input:  f5 = "RAW/{rawfile}.fast5"
    output: fq = "FASTQ/{rawfile}.fastq.gz"
    log:    expand("LOGS/{sets}/{cape}.call.log", sets=SETS, cape=CALLERENV)
    conda:  "nextsnakes/envs/"+CALLERENV+".yaml"
    threads: MAXTHREAD
    params: caller = CALLERBIN,
            cpara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None, config, 'CALLING')['OPTIONS'][0].items()),
            f5dir = lambda wildcards, input: os.path.dirname(input.f5)
    shell: "{params.caller} {params.cpara} --verbose_logs --compress_fastq -i {input.f5} -s {params.f5dir} {output.fq} 2> {log}"
