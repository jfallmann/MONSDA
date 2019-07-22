rule generate_index:
    input: fa = expand("{ref}/{{gen}}{{name}}.fa.gz", ref=config["REFERENCE"])
    output: idx = expand("{ref}/{{gen}}{{name}}_{map}.idx", ref=config["REFERENCE"], map=MAPPERBIN)
    log:    expand("LOGS/{ref}/{{gen}}{{name}}_idx.log",ref=config["REFERENCE"])
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: 20
    params: mapp = MAPPERBIN,
            anno = lambda wildcards: "{annotation}".format(annotation=os.path.join(config["REFERENCE"],os.path.dirname(wildcards.gen),config["ANNOTATION"][os.path.dirname(wildcards.gen)])),
            genpath = lambda wildards: "{refe}/{genname}".format(refe=config["REFERENCE"],genname=os.path.dirname(wildcards.gen))
    shell: "mkdir -p STARTMP && {params.mapp} --runThreadN {threads} --runMode genomeGenerate --outTmpDir STARTMP --genomeDir {params.genpath} --genomeFastaFiles {input.fa} --sjdbGTFfile {params.anno} --sjdbGTFtagExonParentTranscript Parent  2> {log}"

rule mapping:
    input:  query = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz",
            index = lambda wildcards: "{ref}/{gen}{name}_{map}.idx".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=''.join(tool_params(wildcards.file, None ,config, "MAPPING")[1]), map=MAPPERBIN)
    output: mapped = report("MAPPED/{file}_mapped.sam", category="MAPPING"),
            unmapped = "UNMAPPED/{file}_unmapped.fastq.gz"
    log:    "LOGS/{file}/mapping.log"
    conda:  "../envs/"+MAPPERENV+".yaml"
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: 20
    params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "MAPPING")[0].items()),
            ref = lambda wildcards: check_ref("{ref}/{gen}{name}.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config))),
            mapp=MAPPERBIN,
            genpath = lambda wildards: "{refe}/{genname}".format(refe=config["REFERENCE"],genname=os.path.dirname(wildcards.gen))
    shell: "{params.mapp} {params.mpara}  --runThreadN {threads} --genomeDir {params.genpath} --readFilesCommand zcat --readFilesIn {input.query} 2> {log}"
