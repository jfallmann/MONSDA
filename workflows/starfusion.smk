FBIN, FENV = env_bin_from_config(config, 'FUSIONS')

if not 'star' in combo or not 'star' in scombo:
        log.warning('STAR-Fusion needs STAR chimeric output, can only be used with STAR in mapping step')

fusopts = tool_params(SAMPLES[0], None, config, "FUSIONS", FENV)['OPTIONS']
fastqmode = str(fusopts.get('FASTQ', "")).lower() in ("1", "true", "yes")
CTATLIB = fusopts.get('INDEX', "") or os.path.join(REFDIR, "CTAT", FENV)

if not rundedup:
    rule themall:
        input:  expand("FUSIONS/{combo}/{file}_fusions", combo=combo, file=samplecond(SAMPLES, config))
else:
    rule themall:
        input:  expand("FUSIONS/{combo}/{file}_{type}", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_dedup'])

rule generate_ctat_lib:
    input:  fa = REFERENCE,
            anno = ANNOTATION
    output: lib = directory(CTATLIB),
            tmpfa = temp(expand("TMP/{fenv}/ctat_ref.fa", fenv=FENV)),
            tmpanno = temp(expand("TMP/{fenv}/ctat_ref.gtf", fenv=FENV))
    log:    expand("LOGS/{sets}/{fenv}.ctat.log", sets=SETS, fenv=FENV)
    conda:  ""+FENV+".yaml"
    container: "oras://jfallmann/monsda:"+FENV+""
    threads: MAXTHREAD
    params: sf = FBIN,
            bpara = lambda wildcards: tool_params(SAMPLES[0], None, config, "FUSIONS", FENV)['OPTIONS'].get('BUILD', "")
    shell:  "if [[ -f \"{output.lib}/ref_annot.gtf\" || -f \"{output.lib}/ctat_genome_lib_build_dir/ref_annot.gtf\" ]]; then echo \"Found CTAT genome lib, continue with STAR-Fusion\" &> {log} && mkdir -p {output.lib} && touch {output.tmpfa} {output.tmpanno}; else ( zcat {input.fa} > {output.tmpfa} && zcat {input.anno} > {output.tmpanno} && mkdir -p {output.lib} && prep_genome_lib.pl --genome_fa {output.tmpfa} --gtf {output.tmpanno} --output_dir {output.lib} --CPU {threads} {params.bpara} ) &> {log}; fi"

if not fastqmode:
    rule starfusion:
        input:  junction = expand("MAPPED/{scombo}/{{file}}.Chimeric.out.junction", scombo=scombo),
                lib = rules.generate_ctat_lib.output.lib
        output: fusions = "FUSIONS/{combo}/{file}_fusions",
                outdir = temp(directory("FUSIONS/{combo}/{file}_tmp"))
        log:    "LOGS/FUSIONS/{combo}/{file}_starfusion.log"
        conda:  ""+FENV+".yaml"
        container: "oras://jfallmann/monsda:"+FENV+""
        threads: MAXTHREAD
        params: fpara = lambda wildcards: tool_params(wildcards.file, None, config, "FUSIONS", FENV)['OPTIONS'].get('FUSION', ""),
                sf = FBIN
        shell:  "set +o pipefail; GLIB={input.lib}; if [[ -f \"$GLIB/ctat_genome_lib_build_dir/ref_annot.gtf\" ]]; then GLIB=$GLIB/ctat_genome_lib_build_dir; fi; if [[ -s \"{input.junction}\" ]]; then {params.sf} --genome_lib_dir $GLIB -J {input.junction} --output_dir {output.outdir} --CPU {threads} {params.fpara} &> {log} && cp {output.outdir}/star-fusion.fusion_predictions.tsv {output.fusions} 2>> {log}; else mkdir -p {output.outdir}; echo \"File {input.junction} empty, no chimeric STAR output found\" >> {log}; fi; touch {output.fusions}"
else:
    if paired == 'paired':
        rule starfusion:
            input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_R1_trimmed.fastq.gz",
                    r2 = "TRIMMED_FASTQ/{combo}/{file}_R2_trimmed.fastq.gz",
                    lib = rules.generate_ctat_lib.output.lib
            output: fusions = "FUSIONS/{combo}/{file}_fusions",
                    outdir = temp(directory("FUSIONS/{combo}/{file}_tmp"))
            log:    "LOGS/FUSIONS/{combo}/{file}_starfusion.log"
            conda:  ""+FENV+".yaml"
            container: "oras://jfallmann/monsda:"+FENV+""
            threads: MAXTHREAD
            params: fpara = lambda wildcards: tool_params(wildcards.file, None, config, "FUSIONS", FENV)['OPTIONS'].get('FUSION', ""),
                    sf = FBIN
            shell:  "GLIB={input.lib}; if [[ -f \"$GLIB/ctat_genome_lib_build_dir/ref_annot.gtf\" ]]; then GLIB=$GLIB/ctat_genome_lib_build_dir; fi; {params.sf} --genome_lib_dir $GLIB --left_fq {input.r1} --right_fq {input.r2} --output_dir {output.outdir} --CPU {threads} {params.fpara} &> {log} && cp {output.outdir}/star-fusion.fusion_predictions.tsv {output.fusions} 2>> {log}; touch {output.fusions}"
    else:
        rule starfusion:
            input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_trimmed.fastq.gz",
                    lib = rules.generate_ctat_lib.output.lib
            output: fusions = "FUSIONS/{combo}/{file}_fusions",
                    outdir = temp(directory("FUSIONS/{combo}/{file}_tmp"))
            log:    "LOGS/FUSIONS/{combo}/{file}_starfusion.log"
            conda:  ""+FENV+".yaml"
            container: "oras://jfallmann/monsda:"+FENV+""
            threads: MAXTHREAD
            params: fpara = lambda wildcards: tool_params(wildcards.file, None, config, "FUSIONS", FENV)['OPTIONS'].get('FUSION', ""),
                    sf = FBIN
            shell:  "GLIB={input.lib}; if [[ -f \"$GLIB/ctat_genome_lib_build_dir/ref_annot.gtf\" ]]; then GLIB=$GLIB/ctat_genome_lib_build_dir; fi; {params.sf} --genome_lib_dir $GLIB --left_fq {input.r1} --output_dir {output.outdir} --CPU {threads} {params.fpara} &> {log} && cp {output.outdir}/star-fusion.fusion_predictions.tsv {output.fusions} 2>> {log}; touch {output.fusions}"
