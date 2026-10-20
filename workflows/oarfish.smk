COUNTBIN, COUNTENV = env_bin_from_config(config,'COUNTING')
keydict = sub_dict(tool_params(SAMPLES[0], None, config, 'COUNTING', COUNTENV)['OPTIONS'], ['INDEX'])
keydict["REF"] = REFERENCE
keydict["DECOY"] = DECOY
keydict["ENV"] = COUNTENV
unik = get_dict_hash(keydict)

rule themall:
    input:  expand("COUNTS/{combo}/{file}_counts.sf.gz", combo=combo, file=samplecond(SAMPLES, config))

rule decompress_annotation:
    input:  anno = ANNOTATION
    output: anno = temp("COUNTS/{combo}/annotation.gtf")
    log:    "LOGS/{combo}/COUNTING/oarfish/decompress_annotation.log"
    conda:  ""+COUNTENV+".yaml"
    container: "oras://jfallmann/monsda:"+COUNTENV+""
    threads: 1
    shell:  "zcat -f {input.anno} > {output.anno} 2> {log}; if [ ! -s {output.anno} ]; then echo 'ERROR: decompressed annotation from {input.anno} is empty, check source file and filesystem' >> {log}; exit 1; fi"

rule mapping:
    input:  bam = expand("MAPPED/{scombo}/{{file}}_mapped_sorted.bam", scombo=scombo),
            anno = rules.decompress_annotation.output.anno,
            fa = REFERENCE
    output: cnts = report("COUNTS/{combo}/{file}_counts.sf.gz", category="COUNTING"),
            ctsdir = report(directory("COUNTS/{combo}/{file}"), category="COUNTING")
    log:    "LOGS/{combo}/{file}/COUNTING/oarfish/oarfishquant.log"
    conda:  ""+COUNTENV+".yaml"
    container: "oras://jfallmann/monsda:"+COUNTENV+""
    threads: MAXTHREAD
    params: cpara = lambda wildcards: tool_params(wildcards.file, None, config, 'COUNTING', COUNTENV)['OPTIONS'].get('COUNT', ""),
            mapp = COUNTBIN,
            nbam = lambda wildcards, output: os.path.join(str(output.ctsdir), "namecollated.bam"),
            prefix = lambda wildcards, output: os.path.join(str(output.ctsdir), "oarfish"),
            linksf = lambda wildcards, output: str(os.path.abspath(output.ctsdir)),
            bins = BINS
    shell: "set +euo pipefail; mkdir -p {output.ctsdir}; samtools collate -@ {threads} -o {params.nbam} {input.bam} &>> {log}; {params.mapp} --threads {threads} --genome-alignments {params.nbam} --annotation {input.anno} --genome-fasta {input.fa} {params.cpara} --output {params.prefix} &>> {log} && rm -f {params.nbam} && python3 {params.bins}/Analysis/oarfish_to_salmon.py -i {params.prefix}.quant -o {output.ctsdir}/quant.sf.gz &>> {log} ; ln -fs {params.linksf}/quant.sf.gz {output.cnts} &>> {log}"
