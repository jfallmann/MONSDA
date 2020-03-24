RAWBIN, RAWENV = env_bin_from_config2(SAMPLES,config,'RAW')

def parse_sra_list:             # Do this


rule get_from_sra:
    input:  list = expand("{ref}/{{dir}}/{{gen}}{{name}}.sra.gz", ref=REFERENCE)
    output: idx = expand("{ref}/{{dir}}/{map}/{{gen}}{{name}}_{{extension}}_{map}.idx", ref=REFERENCE, map=MAPPERENV)
    log:    expand("LOGS/{{dir}}/{{gen}}{{name}}_{{extension}}_{map}.idx.log", map=MAPPERENV)
    conda:  "snakes/envs/"+MAPPERENV+".yaml"
    threads: MAXTHREAD
    params: indexer = MAPPERBIN.split(' ')[0],
            ipara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None, config, 'MAPPING')['OPTIONS'][0].items())
    shell:  "fastq-dump --split-files --origfmt --gzip "
