rule mapping:
    input:  expand("TRIMMED_FASTQ/{rawfile}_trimmed.fastq.gz",rawfile=SAMPLES)
    output: report("MAPPED/{file}_mapped.sam", category="MAPPING"),
            "UNMAPPED/{file}_unmapped.fastq.gz"
    log:    "LOGS/{file}/mapping.log"
    conda:  "../envs/"+MAPPERENV+".yaml"
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: 20
    params: mpara = lambda wildcards: ''.join(mapping_params(wildcards.file, None ,config)[0]),
            index = lambda wildcards: "{ref}/{gen}{name}.idx".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=''.join(mapping_params(wildcards.file, None ,config)[1])),
            ref = lambda wildcards: check_ref("{ref}/{gen}{name}.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=NAME['genomic'])),
            mapp=MAPPERBIN
    shell: "{params.mapp} {params.mpara} -d {params.ref} -i {params.index} -q {input[0]} --threads {threads} -o {output[0]} -u {output[1]} 2> {log}"
