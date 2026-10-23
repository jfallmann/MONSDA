MAPPERBIN, MAPPERENV = env_bin_from_config(config,'MAPPING')
keydict = sub_dict(tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'], ['INDEX'])
keydict["REF"] = REFERENCE
keydict["DECOY"] = DECOY
keydict["ENV"] = MAPPERENV
unik = get_dict_hash(keydict)

rule generate_index:
    input:  fa = REFERENCE
    output: idx = directory(INDEX),
            uidx = directory(expand("{refd}/INDICES/{mape}_{unikey}", refd=REFDIR, mape=MAPPERENV, unikey=unik)),
            idxfile = expand("{refd}/INDICES/{mape}_{unikey}/{pref}", refd=REFDIR, mape=MAPPERENV, unikey=unik, pref = PREFIX),
            tmp = temp(expand("TMP/{mape}/ref.fa", mape=MAPPERENV))
    log:    expand("LOGS/{sets}/MAPPING/{mape}/idx.log", sets=SETS, mape=MAPPERENV)
    conda:  ""+MAPPERENV+".yaml"
    container: "oras://jfallmann/monsda:"+MAPPERENV+""
    threads: MAXTHREAD
    params: indexer = MAPPERBIN.split(' ')[0],
            ipara = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('INDEX', ""),
            lnkidx = lambda wildcards, output: str(os.path.abspath(output.uidx[0])),
            pref = PREFIX
    shell:  "if [[ -f \"{output.idxfile}\" ]]; then ln -fs {params.lnkidx} {output.idx} && mkdir -p {output.uidx} && touch {output.tmp} {output.idxfile} && echo \"Found bowtie2 index, continue with mapping\" ; else zcat {input.fa} > {output.tmp} && mkdir -p {output.uidx} ;{params.indexer}-build {params.ipara} --threads {threads} {output.tmp} {output.uidx}/{params.pref} 2> {log} && ln -fs {params.lnkidx} {output.idx} && touch {output.uidx} && touch {output.idxfile};fi"

# NOTE: unlike hisat2, bowtie2 has no --rna-strandness/stranded-protocol flag,
# so MAPPING_STRANDED is intentionally unused here.

if MAPPING_PAIRED == 'paired':
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_R1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{combo}/{file}_R2_trimmed.fastq.gz",
                uidx = rules.generate_index.output.uidx[0],
                dummy = rules.generate_index.output.idxfile,
                ref = REFERENCE
        output: mapped = temp(report("MAPPED/{combo}/{file}_mapped.sam.gz", category="MAPPING")),
                unmapped_r1 = "UNMAPPED/{combo}/{file}_R1_unmapped.fastq.gz",
                unmapped_r2 = "UNMAPPED/{combo}/{file}_R2_unmapped.fastq.gz"
        log:    "LOGS/{combo}/{file}/MAPPING/bowtie2/mapping.log"
        conda:  ""+MAPPERENV+".yaml"
        container: "oras://jfallmann/monsda:"+MAPPERENV+""
        threads: MAXTHREAD
        params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', ""),
                mapp=MAPPERBIN,
                um = lambda wildcards, output: output.unmapped_r1.replace('_R1_unmapped.fastq.gz', '.unmapped.gz'),
                umn = lambda wildcards, output: output.unmapped_r1.replace('_R1_unmapped.fastq.gz', ''),
                muz = lambda wildcards, output: output.mapped.replace('.gz', ''),
                pref = PREFIX
        shell: "set +euo pipefail; {params.mapp} {params.mpara} -p {threads} -x {input.uidx}/{params.pref} -1 {input.r1} -2 {input.r2} --un-conc-gz {params.umn}.unmapped -S {params.muz} &> {log} && gzip {params.muz} && touch {params.umn}.unmapped.1.gz {params.umn}.unmapped.2.gz; rename 's/.unmapped.([1|2]).gz/_R$1_unmapped.fastq.gz/' {params.umn}.unmapped.*.gz; touch {output.unmapped_r1} {output.unmapped_r2} &>> {log}"

else:
    rule mapping:
        input:  query = ("TRIMMED_FASTQ/{combo}/{file}_R1_trimmed.fastq.gz" if (MAPPING_PAIRED != 'paired' and paired == 'paired') else "TRIMMED_FASTQ/{combo}/{file}_trimmed.fastq.gz"),
                uidx = rules.generate_index.output.uidx[0],
                dummy = rules.generate_index.output.idxfile,
                ref = REFERENCE
        output: mapped = temp(report("MAPPED/{combo}/{file}_mapped.sam.gz", category="MAPPING")),
                unmapped = "UNMAPPED/{combo}/{file}_unmapped.fastq.gz"
        log:    "LOGS/{combo}/{file}/MAPPING/bowtie2/mapping.log"
        conda:  ""+MAPPERENV+".yaml"
        container: "oras://jfallmann/monsda:"+MAPPERENV+""
        threads: MAXTHREAD
        params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', ""),
                mapp=MAPPERBIN,
                muz = lambda wildcards, output: output.mapped.replace('.gz', ''),
                pref = PREFIX
        shell: "set +euo pipefail; {params.mapp} {params.mpara} -p {threads} -x {input.uidx}/{params.pref} -U {input.query} --un-gz {output.unmapped} -S {params.muz} &> {log} && gzip {params.muz} && touch {output.unmapped}"
