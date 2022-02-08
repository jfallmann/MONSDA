CBIN, CENV = env_bin_from_config3(config, 'CIRCS')

if not 'bwa' in combo or not 'bwa' in scombo:
        log.warning('Ciri2 needs BWA input, can only be used with BWA in mapping step')

if not rundedup:
    rule themall:
        input:  expand("CIRCS/{combo}/{file}_circs", combo=combo, file=samplecond(SAMPLES, config))
else:
    rule themall:
        input:  expand("CIRCS/{combo}/{file}_{type}", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_dedup'])

rule FindCircs:
    input:  sam = expand("MAPPED/{scombo}/{{file}}_mapped_sorted.sam.gz", scombo=scombo),
            ref = REFERENCE,
            anno = ANNOTATION
    output: circs = "CIRCS/{combo}/{file}_circs",
            tmp = temp(directory("CIRCS/{combo}/{file}_tmp")),
            ts = temp("CIRCS/{combo}/{file}_tmp.sam"),
            ta = temp("CIRCS/{combo}/{file}_tmp.gtf"),
            tf = temp("CIRCS/{combo}/{file}_tmp.fa")
    log:    "LOGS/CIRCS/{combo}/{file}_ciri2.log"
    conda:  ""+CENV+".yaml"
    threads: MAXTHREAD
    params: cpara = lambda wildcards: tool_params(wildcards.file, None, config, "CIRCS", CENV)['OPTIONS'].get('CIRC', ""),
            circ = CBIN
    shell:  "set +o pipefail; export LC_ALL=C; if [[ -n \"$(zcat {input.sam} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then mkdir -p {output.tmp} && zcat {input.sam}|samtools sort -n -@ {threads} -u -m 20G -O sam -T {output.tmp} > {output.ts} && zcat {input.anno} > {output.ta} && zcat {input.ref} > {output.tf} && perl {params.circ} -I {output.ts} -O {output.circs} -F {output.tf} -T {threads} -A {output.ta} -G {log} {params.cpara} &>> {log}; else gzip < /dev/null > {output.circs}; echo \"File {input.sam} empty\" >> {log}; fi; touch CIRIerror.log && cat CIRIerror.log >> {log} && echo '' > CIRIerror.log && touch {output.circs} && mkdir -p {output.tmp}"
