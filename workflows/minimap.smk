rule mapping:
    input:  expand("TRIMMED_FASTQ/{rawfile}_trimmed.fastq.gz",rawfile=SAMPLES)
    output: report("MAPPED/{file}_mapped.sam", category="MAPPING"),
            "UNMAPPED/{file}_unmapped.fastq.gz"
    log:    "LOGS/{file}/mapping.log"
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: 20
    params: mpara = lambda wildcards: ''.join(mapping_params(wildcards.file, None ,config)[0]),
            index = lambda wildcards: "{ref}/{gen}{name}.idx".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=''.join(mapping_params(wildcards.file, None ,config)[1])),
            ref = lambda wildcards: check_ref("{ref}/{gen}{name}.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=NAME['genomic'])),
            mapp=MAPPERBIN
    shell: "{params.mapp} -t {threads} {params.mpara} {params.index} {params.ref} {input[0]} | tee >(grep -v -P '\t4\t' > {output[0]}) >(grep -P '[^@|\t4\t]' |samtools fastq -n - | pigz > {output[1]}) 1>/dev/null 2>> {log}"
