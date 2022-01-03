CBIN, CENV = env_bin_from_config3(config, 'CIRCS')

if not 'bwa' in combo or not 'bwa' in scombo:
        log.warning('Ciri2 needs BWA input, can only be used with BWA in mapping step')


if not rundedup:
    rule themall:
        input:  expand("CIRCS/{combo}/CIRI2/{file}_circs", combo=combo, file=samplecond(SAMPLES, config))
else:
    rule themall:
        input:  expand("CIRCS/{combo}/CIRI2/{file}_{type}", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_dedup'])

rule FindCircs:
    input:  sam = expand("MAPPED/{scombo}/{{file}}_mapped_sorted.sam.gz", scombo=scombo),
            ref = REFERENCE,
            anno = ANNOTATION
    output: circs = "CIRCS/{combo}/CIRI2/{file}_circs",
            tmp = temp(directory("CIRCS/{combo}/CIRI2/TMP/{file}")),
            ts = temp("CIRCS/{combo}/CIRI2/{file}_tmp.sam")
    log:    "LOGS/CIRCS/{combo}/{file}_ciri2.log"
    conda:  ""+CENV+".yaml"
    threads: MAXTHREAD
    params: cpara = lambda wildcards: tool_params(wildcards.file, None, config, "CIRCS", CENV)['OPTIONS'].get('CIRC', ""),
            circ = CBIN,
            outdir = lambda wildcards, output: os.path.dirname(output.peak),
            outname = lambda wildcards: os.path.basename(wildcards.file)+'_peak_'+wildcards.type,
    shell:  "set +o pipefail; export LC_ALL=C; if [[ -n \"$(samtools view {input.sam} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then gunzip -c {input.sam}|samtools sort -n -@ {threads} -u -m 20G -O sam -T {output.tmp} > {output.ts} && perl {params.circ} -I {output.ts} -O {output.circs} -F <(zcat {input.ref}) -T {threads} -A <(zcat {input.anno}) -G {log} {params.cpara}; else gzip < /dev/null > {output.circs}; echo \"File {input.sam} empty\" >> {log}; fi"
