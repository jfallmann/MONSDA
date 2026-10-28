MAPPERBIN, MAPPERENV = env_bin_from_config(config,'MAPPING')
keydict = sub_dict(tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'], ['INDEX'])
keydict["REF"] = REFERENCE
keydict["DECOY"] = DECOY
keydict["ENV"] = MAPPERENV
unik = get_dict_hash(keydict)

# piscem needs a sub-command. MONSDA maps bulk libraries, single-cell mapping
# is only used when the user asks for it by passing a --geometry to MAP.
PISCEMSUB = "map-sc" if "--geometry" in str(tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', "")) else "map-bulk"

rule generate_index:
    input:  fa = REFERENCE
    output: idx = directory(INDEX),
            uidx = directory(expand("{refd}/INDICES/{mape}_{unikey}", refd=REFDIR, mape=MAPPERENV, unikey=unik)),
            idxfile = expand("{refd}/INDICES/{mape}_{unikey}/{pref}.sshash", refd=REFDIR, mape=MAPPERENV, unikey=unik, pref=PREFIX),
            tmpfa = temp(expand("TMP/{mape}/ref.fa", mape=MAPPERENV))
    log:    expand("LOGS/{sets}/MAPPING/{mape}/idx.log", sets=SETS, mape=MAPPERENV)
    conda:  ""+MAPPERENV+".yaml"
    container: "oras://jfallmann/monsda:"+MAPPERENV+""
    threads: MAXTHREAD
    params: indexer = MAPPERBIN,
            ipara = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('INDEX', ""),
            pref = PREFIX,
            lnkidx = lambda wildcards, output: str(os.path.abspath(output.uidx[0]))
    shell:  "if [[ -f \"{output.idxfile}\" ]]; then touch {output.idxfile} && ln -fs {params.lnkidx} {output.idx} && echo \"Found index, continue with mapping\" ; else zcat {input.fa} > {output.tmpfa} && mkdir -p {output.uidx} && {params.indexer} build -t {threads} {params.ipara} -s {output.tmpfa} -o {output.uidx}/{params.pref} &> {log} && touch {output.idxfile} && ln -fs {params.lnkidx} {output.idx};fi"

if PISCEMSUB == "map-sc":  # -o names a directory

    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_R1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{combo}/{file}_R2_trimmed.fastq.gz",
                uidx = rules.generate_index.output.uidx[0],
                dummy = rules.generate_index.output.idxfile
        output: rad = directory("MAPPED/{combo}/{file}"),
                radfile = "MAPPED/{combo}/{file}/map.rad"
        log:    "LOGS/{combo}/{file}/MAPPING/piscem/mapping.log"
        conda:  ""+MAPPERENV+".yaml"
        container: "oras://jfallmann/monsda:"+MAPPERENV+""
        threads: MAXTHREAD
        params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', ""),
                mapp = MAPPERBIN,
                pref = PREFIX,
                idxpref = lambda wildcards, input: os.path.join(str(input.uidx), PREFIX),
                subcmd = PISCEMSUB
        shell: "{params.mapp} {params.subcmd} -t {threads} {params.mpara} -i {params.idxpref} -1 {input.r1} -2 {input.r2} -o {output.rad} &> {log}"

else:  # map-bulk, -o is a file stem and writes <stem>.rad and <stem>.map_info.json

    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_R1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{combo}/{file}_R2_trimmed.fastq.gz",
                uidx = rules.generate_index.output.uidx[0],
                dummy = rules.generate_index.output.idxfile
        output: radfile = "MAPPED/{combo}/{file}.rad",
                mapinfo = "MAPPED/{combo}/{file}.map_info.json"
        log:    "LOGS/{combo}/{file}/MAPPING/piscem/mapping.log"
        conda:  ""+MAPPERENV+".yaml"
        container: "oras://jfallmann/monsda:"+MAPPERENV+""
        threads: MAXTHREAD
        params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', ""),
                mapp = MAPPERBIN,
                pref = PREFIX,
                idxpref = lambda wildcards, input: os.path.join(str(input.uidx), PREFIX),
                subcmd = PISCEMSUB,
                stem = lambda wildcards, output: output.radfile[:-len(".rad")]
        shell: "mkdir -p $(dirname {params.stem}) && {params.mapp} {params.subcmd} -t {threads} {params.mpara} -i {params.idxpref} -1 {input.r1} -2 {input.r2} -o {params.stem} &> {log}"
