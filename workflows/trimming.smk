subworkflow sampleqc:
    snakefile: lambda wildcards: create_subworkflow(wildcards.file, config, "TRIMMING")
    configfile: lambda wildcards: create_subconfig(wildcards.file, config, "TRIMMING")

rule all_trim:
    input: sampleqc("TRIMMED_FASTQ/{file}_trimmed.fastq.gz")
