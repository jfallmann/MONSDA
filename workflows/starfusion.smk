FBIN, FENV = env_bin_from_config(config, 'FUSIONS')

if not 'star' in combo or not 'star' in scombo:
        log.warning('STAR-Fusion needs STAR chimeric output, can only be used with STAR in mapping step')

fusopts = tool_params(SAMPLES[0], None, config, "FUSIONS", FENV)['OPTIONS']
fastqmode = str(fusopts.get('FASTQ', "")).lower() in ("1", "true", "yes")
CTATLIB = fusopts.get('INDEX', "") or os.path.join(REFDIR, "CTAT", FENV)

if os.path.isfile(os.path.join(CTATLIB, "ref_annot.gtf")):
    CTATGENOMEDIR = CTATLIB
elif os.path.isfile(os.path.join(CTATLIB, "ctat_genome_lib_build_dir", "ref_annot.gtf")):
    CTATGENOMEDIR = os.path.join(CTATLIB, "ctat_genome_lib_build_dir")
else:
    CTATGENOMEDIR = None

if not rundedup:
    rule themall:
        input:  expand("FUSIONS/{combo}/{file}/star-fusion.fusion_predictions.tsv", combo=combo, file=samplecond(SAMPLES, config))
else:
    rule themall:
        input:  expand("FUSIONS/{combo}/{file}/star-fusion.fusion_predictions.tsv", combo=combo, file=samplecond(SAMPLES, config))

rule generate_ctat_lib:
    input:  fa = REFERENCE,
            anno = ANNOTATION
    output: lib = directory(CTATLIB),
            tmpfa = temp(expand("TMP/{fenv}/ctat_ref.fa", fenv=FENV)),
            tmpanno = temp(expand("TMP/{fenv}/ctat_ref.gtf", fenv=FENV))
    log:    expand("LOGS/{sets}/FUSIONS/{fenv}/ctat.log", sets=SETS, fenv=FENV)
    conda:  ""+FENV+".yaml"
    container: "oras://jfallmann/monsda:"+FENV+""
    threads: MAXTHREAD
    params: sf = FBIN,
            bpara = lambda wildcards: tool_params(SAMPLES[0], None, config, "FUSIONS", FENV)['OPTIONS'].get('BUILD', "")
    shell:  "( zcat {input.fa} > {output.tmpfa} && zcat {input.anno} > {output.tmpanno} && mkdir -p {output.lib} && prep_genome_lib.pl --genome_fa {output.tmpfa} --gtf {output.tmpanno} --output_dir {output.lib} --CPU {threads} {params.bpara} ) &> {log}"

if CTATGENOMEDIR is not None:
    ctat_lib_input = CTATGENOMEDIR
else:
    ctat_lib_input = rules.generate_ctat_lib.output.lib


if not fastqmode:
    rule starfusion:
        input:  junction = expand("MAPPED/{scombo}/{{file}}.Chimeric.out.junction", scombo=scombo),
                lib = ctat_lib_input
        output: preds = "FUSIONS/{combo}/{file}/star-fusion.fusion_predictions.tsv",
                abridged = "FUSIONS/{combo}/{file}/star-fusion.fusion_predictions.abridged.tsv",
                normjunc = "FUSIONS/{combo}/{file}.Chimeric.norm.junction"
        log:    "LOGS/{combo}/{file}/FUSIONS/starfusion/starfusion.log"
        conda:  ""+FENV+".yaml"
        container: "oras://jfallmann/monsda:"+FENV+""
        threads: MAXTHREAD
        params: fpara = lambda wildcards: tool_params(wildcards.file, None, config, "FUSIONS", FENV)['OPTIONS'].get('FUSION', ""),
                sf = FBIN
        shell:  "set +o pipefail; OUTDIR=$(dirname {output.preds}); if [[ -s \"{input.junction}\" ]]; then ctat_chr=$(grep -v '^#' {input.lib}/ref_annot.gtf 2>/dev/null | head -1 | cut -f1 | grep -c '^chr' || true); junc_chr=$(awk '!/^#/ && $1!=\"chr_donorA\"{{print $1; exit}}' {input.junction} | grep -c '^chr' || true); if [[ \"$ctat_chr\" == \"1\" && \"$junc_chr\" == \"0\" ]]; then echo \"Adding chr prefix to junction to match CTAT lib\" >> {log}; awk 'BEGIN{{OFS=\"\\t\"}} /^#/{{print;next}} $1==\"chr_donorA\"{{print;next}} {{if($1!~/^chr/)$1=\"chr\"$1; if($4!~/^chr/)$4=\"chr\"$4; if($1==\"chrMT\")$1=\"chrM\"; if($4==\"chrMT\")$4=\"chrM\"; print}}' {input.junction} > {output.normjunc}; elif [[ \"$ctat_chr\" == \"0\" && \"$junc_chr\" == \"1\" ]]; then echo \"Stripping chr prefix from junction to match CTAT lib\" >> {log}; awk 'BEGIN{{OFS=\"\\t\"}} /^#/{{print;next}} $1==\"chr_donorA\"{{print;next}} {{sub(/^chr/,\"\",$1); sub(/^chr/,\"\",$4); if($1==\"M\")$1=\"MT\"; if($4==\"M\")$4=\"MT\"; print}}' {input.junction} > {output.normjunc}; else cp {input.junction} {output.normjunc}; fi; {params.sf} --genome_lib_dir {input.lib} -J {output.normjunc} --output_dir $OUTDIR --CPU {threads} {params.fpara} &> {log}; else mkdir -p $OUTDIR; touch {output.normjunc}; echo \"File {input.junction} empty, no chimeric STAR output found\" >> {log}; fi; touch {output.preds} {output.abridged}"
else:
    if paired == 'paired':
        rule starfusion:
            input:  r1 = expand("TRIMMED_FASTQ/{scombo}/{{file}}_R1_trimmed.fastq.gz", scombo=scombo),
                    r2 = expand("TRIMMED_FASTQ/{scombo}/{{file}}_R2_trimmed.fastq.gz", scombo=scombo),
                    lib = ctat_lib_input
            output: preds = "FUSIONS/{combo}/{file}/star-fusion.fusion_predictions.tsv",
                    abridged = "FUSIONS/{combo}/{file}/star-fusion.fusion_predictions.abridged.tsv"
            log:    "LOGS/{combo}/{file}/FUSIONS/starfusion/starfusion.log"
            conda:  ""+FENV+".yaml"
            container: "oras://jfallmann/monsda:"+FENV+""
            threads: MAXTHREAD
            params: fpara = lambda wildcards: tool_params(wildcards.file, None, config, "FUSIONS", FENV)['OPTIONS'].get('FUSION', ""),
                    sf = FBIN
            shell:  "OUTDIR=$(dirname {output.preds}); {params.sf} --genome_lib_dir {input.lib} --left_fq {input.r1} --right_fq {input.r2} --output_dir $OUTDIR --CPU {threads} {params.fpara} &> {log}; touch {output.preds} {output.abridged}"
    else:
        rule starfusion:
            input:  r1 = expand("TRIMMED_FASTQ/{scombo}/{{file}}_trimmed.fastq.gz", scombo=scombo),
                    lib = ctat_lib_input
            output: preds = "FUSIONS/{combo}/{file}/star-fusion.fusion_predictions.tsv",
                    abridged = "FUSIONS/{combo}/{file}/star-fusion.fusion_predictions.abridged.tsv"
            log:    "LOGS/{combo}/{file}/FUSIONS/starfusion/starfusion.log"
            conda:  ""+FENV+".yaml"
            container: "oras://jfallmann/monsda:"+FENV+""
            threads: MAXTHREAD
            params: fpara = lambda wildcards: tool_params(wildcards.file, None, config, "FUSIONS", FENV)['OPTIONS'].get('FUSION', ""),
                    sf = FBIN
            shell:  "OUTDIR=$(dirname {output.preds}); {params.sf} --genome_lib_dir {input.lib} --left_fq {input.r1} --output_dir $OUTDIR --CPU {threads} {params.fpara} &> {log}; touch {output.preds} {output.abridged}"
