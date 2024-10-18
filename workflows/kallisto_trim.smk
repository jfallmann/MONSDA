COUNTBIN, COUNTENV = env_bin_from_config(config,'COUNTING')
keydict = sub_dict(tool_params(SAMPLES[0], None, config, 'COUNTING', COUNTENV)['OPTIONS'], ['INDEX'])
keydict["REF"] = REFERENCE
keydict["DECOY"] = DECOY
keydict["ENV"] = COUNTENV
unik = get_dict_hash(keydict)

rule themall:
    input:  expand("COUNTS/{combo}/{file}_counts.sf.gz", combo=combo, file=samplecond(SAMPLES, config)),

if paired == 'paired':
    rule simulate_trim:
        input:  r1 = lambda wildcards: "FASTQ/{rawfile}_R1.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]) if not prededup else "DEDUP_FASTQ/{combo}/{file}_R1_dedup.fastq.gz",
                r2 = lambda wildcards: "FASTQ/{rawfile}_R2.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]) if not prededup else "DEDUP_FASTQ/{combo}/{file}_R2_dedup.fastq.gz"
        output: r1 = "TRIMMED_FASTQ/{scombo}/{file}_R1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{scombo}/{file}_R2_trimmed.fastq.gz"
        threads: 1
        params: filetolink = lambda w, input: "{r}".format(r=os.path.abspath(input.r1)),
                filetolink2 = lambda w, input: "{r}".format(r=os.path.abspath(input.r2))
        shell:  "ln -s {params.filetolink} {output.r1} && ln -s {params.filetolink2} {output.r2}"

else:
    rule simulate_trim:
        input:  r1 = lambda wildcards: "FASTQ/{rawfile}.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]) if not prededup else "DEDUP_FASTQ/{combo}/{file}_dedup.fastq.gz"
        output: r1 = "TRIMMED_FASTQ/{scombo}/{file}_trimmed.fastq.gz"
        threads: 1
        params: filetolink = lambda w, input: "{r}".format(r=os.path.abspath(input.r1))
        shell:  "ln -s {params.filetolink} {output.r1}"

rule kallisto_index:
    input:  fa = REFERENCE
    output: idx = INDEX,
            uidx = expand("{refd}/INDICES/{mape}_{unikey}.idx", refd=REFDIR, mape=COUNTENV, unikey=unik)
    log:    expand("LOGS/{sets}/{cape}.idx.log", sets=SETS, cape=COUNTENV)
    conda:  ""+COUNTENV+".yaml"
    container: "docker://jfallmann/monsda:COUNTENV"
    threads: MAXTHREAD
    params: mapp = COUNTBIN,
            ipara = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'COUNTING', COUNTENV)['OPTIONS'].get('INDEX', ""),
            decoy = f"-d {os.path.abspath(DECOY)}" if DECOY else '',
            linkidx = lambda wildcards, output: str(os.path.abspath(output.uidx[0]))
    shell:  "set +euo pipefail; {params.mapp} index {params.ipara} {params.decoy} -t {threads} -i {output.uidx} {input.fa} &>> {log} && ln -fs {params.linkidx} {output.idx}"


if paired == 'paired':
    rule mapping:
        input:  r1 = expand("TRIMMED_FASTQ/{scombo}/{{file}}_R1_trimmed.fastq.gz", scombo=scombo),
                r2 = expand("TRIMMED_FASTQ/{scombo}/{{file}}_R2_trimmed.fastq.gz", scombo=scombo),
                uidx = rules.kallisto_index.output.uidx[0]
        output: cnts = report("COUNTS/{combo}/{file}_counts.gz", category="COUNTING"),
                ctsdir = report(directory("COUNTS/{combo}/{file}"), category="COUNTING")
        log:    "LOGS/{combo}/{file}/kallistoquant.log"
        conda:  ""+COUNTENV+".yaml"
        container: "docker://jfallmann/monsda:COUNTENV"
        threads: MAXTHREAD
        params: cpara = lambda wildcards: tool_params(wildcards.file, None, config, 'COUNTING', COUNTENV)['OPTIONS'].get('COUNT', ""),
                mapp=COUNTBIN,
                stranded = lambda x: '--fr-stranded' if (stranded == 'fr' or stranded == 'ISF') else '-rf-stranded' if (stranded == 'rf' or stranded == 'ISR') else '',
                linksf = lambda wildcards, output: str(os.path.abspath(output.ctsdir))
        shell: "set +euo pipefail; {params.mapp} quant -t {threads} -i {input.uidx} {params.stranded} {params.cpara} -o {output.ctsdir} {input.r1} {input.r2} &>> {log} && gzip {output.ctsdir}/abundance.tsv && ln -fs {params.linksf}/abundance.tsv.gz {output.cnts} &>> {log}"

else:
    rule mapping:
        input:  r1 = expand("TRIMMED_FASTQ/{scombo}/{{file}}_trimmed.fastq.gz", scombo=scombo),
                uidx = rules.kallisto_index.output.uidx[0]
        output: cnts = report("COUNTS/{combo}/{file}_counts.gz", category="COUNTING"),
                ctsdir = report(directory("COUNTS/{combo}/{file}"), category="COUNTING")
        log:    "LOGS/{combo}/{file}/kallistoquant.log"
        conda:  ""+COUNTENV+".yaml"
        container: "docker://jfallmann/monsda:COUNTENV"
        threads: MAXTHREAD
        params: cpara = lambda wildcards: tool_params(wildcards.file, None, config, 'COUNTING', COUNTENV)['OPTIONS'].get('COUNT', ""),
                mapp=COUNTBIN,
                stranded = lambda x: '--fr-stranded' if (stranded == 'fr' or stranded == 'ISF') else '-rf-stranded' if (stranded == 'rf' or stranded == 'ISR') else '',
                linksf = lambda wildcards, output: str(os.path.abspath(output.ctsdir))
        shell: "set +euo pipefail; {params.mapp} quant -t {threads} -i {input.uidx} {params.stranded} {params.cpara} -o {output.ctsdir} --single {input.r1} &>> {log} && gzip {output.ctsdir}/abundance.tsv && ln -fs {params.linksf}/abundance.tsv.gz {output.cnts} &>> {log}"
