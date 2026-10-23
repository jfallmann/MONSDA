if COUNTTOOL == 'salmon':
    rule salmon_index:
        input:  fa = REFERENCE
        output: idx = directory(INDEX),
                uidx = directory(expand("{refd}/INDICES/{mape}_{unikey}", refd=REFDIR, mape=COUNTENV, unikey=unik))
        log:    expand("LOGS/{sets}/DTU/{cape}/idx.log", sets=SETS, cape=COUNTENV)
        conda:  ""+COUNTENV+".yaml"
        container: "oras://jfallmann/monsda:"+COUNTENV+""
        threads: MAXTHREAD
        params: mapp = COUNTBIN,
                ipara = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'DTU', DTUENV)['OPTIONS'].get('INDEX', ""),
                decoy = f"-d {os.path.abspath(DECOY)}" if DECOY else '',
                linkidx = lambda wildcards, output: str(os.path.abspath(output.uidx[0]))
        shell:  "set +euo pipefail; {params.mapp} index {params.ipara} {params.decoy} -p {threads} -t {input.fa} -i {output.uidx} &>> {log} && ln -fs {params.linkidx} {output.idx}"

    if DTU_PAIRED == 'paired':
        rule simulate_trim:
            input:  r1 = lambda wildcards: "FASTQ/{rawfile}_R1.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]) if not usededup else "DEDUP_FASTQ/{scombo}/{file}_R1_dedup.fastq.gz".format(scombo=scombo, file=wildcards.file),
                    r2 = lambda wildcards: "FASTQ/{rawfile}_R2.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]) if not usededup else "DEDUP_FASTQ/{scombo}/{file}_R2_dedup.fastq.gz".format(scombo=scombo, file=wildcards.file)
            output: r1 = "TRIMMED_FASTQ/{scombo}/{file}_R1_trimmed.fastq.gz",
                    r2 = "TRIMMED_FASTQ/{scombo}/{file}_R2_trimmed.fastq.gz"
            threads: 1
            params: filetolink = lambda w, input: "{r}".format(r=os.path.abspath(input.r1)),
                    filetolink2 = lambda w, input: "{r}".format(r=os.path.abspath(input.r2))
            shell:  "ln -s {params.filetolink} {output.r1} && ln -s {params.filetolink2} {output.r2}"

    else:
        rule simulate_trim:
            input:  r1 = lambda wildcards: "FASTQ/{rawfile}.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]) if not usededup else "DEDUP_FASTQ/{scombo}/{file}_dedup.fastq.gz".format(scombo=scombo, file=wildcards.file)
            output: r1 = "TRIMMED_FASTQ/{scombo}/{file}_trimmed.fastq.gz"
            threads: 1
            params: filetolink = lambda w, input: "{r}".format(r=os.path.abspath(input.r1))
            shell:  "ln -s {params.filetolink} {output.r1}"

    if DTU_PAIRED == 'paired':
        rule mapping:
            input:  r1 = expand("TRIMMED_FASTQ/{scombo}/{{file}}_R1_trimmed.fastq.gz", scombo=scombo),
                    r2 = expand("TRIMMED_FASTQ/{scombo}/{{file}}_R2_trimmed.fastq.gz", scombo=scombo),
                    index = rules.salmon_index.output.idx,
                    uix = rules.salmon_index.output.uidx
            output: cnts = report("DTU/{combo}/salmon/{file}_counts.sf.gz", category="COUNTING"),
                    ctsdir = report(directory("DTU/{combo}/salmon/{file}"), category="COUNTING")
            log:    "LOGS/{combo}/{file}/DTU/salmon/salmonquant.log"
            conda:  ""+COUNTENV+".yaml"
            container: "oras://jfallmann/monsda:"+COUNTENV+""
            threads: MAXTHREAD
            params: cpara = lambda wildcards: tool_params(wildcards.file, None, config, 'DTU', DTUENV)['OPTIONS'].get('QUANT', ""),
                    mapp=COUNTBIN,
                    stranded = lambda x: '-l ISF' if (DTU_STRANDED == 'fr' or DTU_STRANDED == 'ISF') else '-l ISR' if (DTU_STRANDED == 'rf' or DTU_STRANDED == 'ISR') else '-l IU',
                    linksf = lambda wildcards, output: str(os.path.abspath(output.ctsdir))
            shell: "set +euo pipefail; {params.mapp} quant -p {threads} -i {input.index} {params.stranded} {params.cpara} -o {output.ctsdir} -1 {input.r1} -2 {input.r2} &>> {log} && gzip {output.ctsdir}/quant.sf && ln -fs {params.linksf}/quant.sf.gz {output.cnts} &>> {log}"

    else:
        rule mapping:
            input:  r1 = expand("TRIMMED_FASTQ/{scombo}/{{file}}_trimmed.fastq.gz", scombo=scombo),
                    index = rules.salmon_index.output.idx,
                    uix = rules.salmon_index.output.uidx
            output: cnts = report("DTU/{combo}/salmon/{file}_counts.sf.gz", category="COUNTING"),
                    ctsdir = report(directory("DTU/{combo}/salmon/{file}"), category="COUNTING")
            log:    "LOGS/{combo}/{file}/DTU/salmon/salmonquant.log"
            conda:  ""+COUNTENV+".yaml"
            container: "oras://jfallmann/monsda:"+COUNTENV+""
            threads: MAXTHREAD
            params: cpara = lambda wildcards: tool_params(wildcards.file, None, config, 'DTU', DTUENV)['OPTIONS'].get('QUANT', ""),
                    mapp=COUNTBIN,
                    stranded = lambda x: '-l SF' if (DTU_STRANDED == 'fr' or DTU_STRANDED == 'SF') else '-l SR' if (DTU_STRANDED == 'rf' or DTU_STRANDED == 'SR') else '-l U',
                    linksf = lambda wildcards, output: str(os.path.abspath(output.ctsdir))
            shell: "set +euo pipefail; {params.mapp} quant -p {threads} -i {input.index} {params.stranded} {params.cpara} -o {output.ctsdir} -r {input.r1} &>> {log} && gzip {output.ctsdir}/quant.sf ; ln -fs {params.linksf}/quant.sf.gz {output.cnts} &>> {log}"

else:
    # Long-read quantifier path (e.g. oarfish): quantify directly from the
    # already generated genome-mapped, sorted BAM. No trimming/index step
    # is needed, output is normalized to the salmon quant.sf schema so
    # downstream DTU tooling (tximport type="salmon", build_DTU_table.py,
    # SPIT.py) does not need to change.
    rule decompress_annotation:
        input:  anno = ANNOTATION
        output: anno = temp("DTU/{combo}/"+COUNTENV+"/annotation.gtf")
        log:    "LOGS/{combo}/DTU/"+COUNTENV+"/decompress_annotation.log"
        conda:  ""+COUNTENV+".yaml"
        container: "oras://jfallmann/monsda:"+COUNTENV+""
        threads: 1
        shell:  "zcat -f {input.anno} > {output.anno} 2> {log}; if [ ! -s {output.anno} ]; then echo 'ERROR: decompressed annotation from {input.anno} is empty, check source file and filesystem' >> {log}; exit 1; fi"

    rule mapping:
        input:  bam = expand("MAPPED/{scombo}/{{file}}_mapped_sorted.bam", scombo=scombo),
                anno = rules.decompress_annotation.output.anno,
                fa = REFERENCE
        output: cnts = report("DTU/{combo}/"+COUNTENV+"/{file}_counts.sf.gz", category="COUNTING"),
                ctsdir = report(directory("DTU/{combo}/"+COUNTENV+"/{file}"), category="COUNTING")
        log:    "LOGS/{combo}/{file}/DTU/"+COUNTENV+"/quant.log"
        conda:  ""+COUNTENV+".yaml"
        container: "oras://jfallmann/monsda:"+COUNTENV+""
        threads: MAXTHREAD
        params: cpara = lambda wildcards: tool_params(wildcards.file, None, config, 'DTU', DTUENV)['OPTIONS'].get('QUANT', ""),
                mapp = COUNTBIN,
                nbam = lambda wildcards, output: os.path.join(str(output.ctsdir), "namecollated.bam"),
                prefix = lambda wildcards, output: os.path.join(str(output.ctsdir), COUNTENV),
                linksf = lambda wildcards, output: str(os.path.abspath(output.ctsdir)),
                bins = BINS
        shell: "set +euo pipefail; mkdir -p {output.ctsdir}; samtools collate -@ {threads} -o {params.nbam} {input.bam} &>> {log}; {params.mapp} --threads {threads} --genome-alignments {params.nbam} --annotation {input.anno} --genome-fasta {input.fa} {params.cpara} --output {params.prefix} &>> {log} && rm -f {params.nbam} && python3 {params.bins}/Analysis/oarfish_to_salmon.py -i {params.prefix}.quant -o {output.ctsdir}/quant.sf.gz &>> {log} ; ln -fs {params.linksf}/quant.sf.gz {output.cnts} &>> {log}"


if runterminus:
    rule terminus_group:
        input:  ctsdir = rules.mapping.output.ctsdir
        output: grp = "DTU/{combo}/terminus/{file}/groups.txt"
        log:    "LOGS/{combo}/{file}/DTU/salmon/terminus_group.log"
        conda:  ""+TERMINUSENV+".yaml"
        container: "oras://jfallmann/monsda:"+TERMINUSENV+""
        threads: MAXTHREAD
        params: tpara = termpara,
                outdir = lambda wildcards, output: os.path.dirname(os.path.dirname(str(output.grp)))
        shell:  "terminus group {params.tpara} -d {input.ctsdir} -o {params.outdir} &> {log}"

    rule terminus_collapse:
        input:  grps = expand("DTU/{combo}/terminus/{file}/groups.txt", combo=combo, file=samplecond(SAMPLES, config)),
                dirs = expand(rules.mapping.output.ctsdir, combo=combo, file=samplecond(SAMPLES, config))
        output: cnts = expand("DTU/{combo}/terminus/{file}/quant.sf.gz", combo=combo, file=samplecond(SAMPLES, config))
        log:    expand("LOGS/{combo}/DTU/salmon/terminus_collapse.log", combo=combo)
        conda:  ""+TERMINUSENV+".yaml"
        container: "oras://jfallmann/monsda:"+TERMINUSENV+""
        threads: MAXTHREAD
        params: outdir = lambda wildcards, input: os.path.dirname(os.path.dirname(str(input.grps[0])))
        shell:  "terminus collapse -d {input.dirs} -o {params.outdir} &> {log}; for d in {input.dirs}; do b=$(basename $d); gzip -c {params.outdir}/$b/quant.sf > {params.outdir}/$b/quant.sf.gz; done 2>> {log}"

    ctsource = expand("DTU/{combo}/terminus/{file}", combo=combo, file=samplecond(SAMPLES, config))
else:
    ctsource = expand(rules.mapping.output.ctsdir, combo=combo, file=samplecond(SAMPLES, config))

rule create_annotation_table:
    input:  dir  = ctsource,
    output: anno = expand("DTU/{combo}/Tables/{scombo}_ANNOTATION.gz", combo=combo, scombo=scombo)
    log:    expand("LOGS/{combo}/DTU/salmon/create_DTU_table.log", combo=combo)
    conda:  ""+COUNTENV+".yaml"
    container: "oras://jfallmann/monsda:"+COUNTENV+""
    threads: 1
    params: dereps = lambda wildcards, input: get_reps(input.dir, config, 'DTU'),
            bins = BINS
    shell:  "python3 {params.bins}/Analysis/build_DTU_table.py {params.dereps} --anno {output.anno} --loglevel DEBUG 2> {log}"
