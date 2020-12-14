logid = 'drimseq_DTU.smk '
DTUBIN, DTUENV = env_bin_from_config3(config,'DTU')
log.info(logid+"DTUENV: "+str(DTUENV))
COUNTBIN, COUNTENV = ['salmon','salmon']#env_bin_from_config2(SAMPLES,config,'COUNTING') ##PINNING subreads package to version 1.6.4 due to changes in 2.0.1 gene_id length cutoff that interfers

outdir = "DTU/DRIMSEQ"
comparison = comparable_as_string2(config,'DTU')
compstr = [i.split(":")[0] for i in comparison.split(",")]
log.info(logid+"COMPARISON: "+str(comparison))

rule themall:
    input:  session = expand("{outdir}/DRIMSEQ_DTU_SESSION.gz", outdir=outdir),
            sig = expand("{outdir}/Sig_DTU_DRIMSEQ_{comparison}_results_genes.tsv.gz", outdir=outdir, comparison=compstr),
            sig_d = expand("{outdir}/SigDOWN_DTU_DRIMSEQ_{comparison}_results_genes.tsv.gz", outdir=outdir, comparison=compstr),
            sig_u = expand("{outdir}/SigUP_DTU_DRIMSEQ_{comparison}_results_genes.tsv.gz", outdir=outdir, comparison=compstr)

            # res_stager = expand("{outdir}/DTU_DRIMSEQ_{comparison}_results_stageR-filtered.tsv.gz", outdir=outdir, comparison=compstr),
            # res_posthoc = expand("{outdir}/DTU_DRIMSEQ_{comparison}_results_post-hoc-filtered-on-SD.tsv.gz", outdir=outdir, comparison=compstr)

rule salmon_index:
    input:  fa = REFERENCE
    output: idx = directory(expand("{refd}/INDICES/{mape}", refd=REFDIR, mape=COUNTENV))
    log:    expand("LOGS/{sets}/salmon.idx.log", sets=SETS)
    conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: mapp = COUNTBIN,
            ipara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None, config, 'DTU')['OPTIONS'][0].items()),
    shell:  "{params.mapp} index {params.ipara} -p {threads} -t {input.fa} -i {output.idx} 2>> {log}"

if paired == 'paired':
    rule mapping:
        input:  r1 = "FASTQ/{file}_R1.fastq.gz",
                r2 = "FASTQ/{file}_R2.fastq.gz",
                index = rules.salmon_index.output.idx
        output: ctsdir = directory("COUNTS/Salmon/{file}")
        log:    "LOGS/{file}/salmonquant.log"
        conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
        threads: MAXTHREAD
        params: cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'DTU')['OPTIONS'][1].items()),
                mapp=COUNTBIN,
                stranded = lambda x: '-l ISF' if (stranded == 'fr' or stranded == 'ISF') else '-l ISR' if (stranded == 'rf' or stranded == 'ISR') else ''
        shell: "{params.mapp} quant -p {threads} -i {input.index} {params.stranded} {params.cpara} -o {output.ctsdir} -1 {input.r1} -2 {input.r2} 2>> {log}"

else:
    rule mapping:
        input:  r1 = "FASTQ/{file}.fastq.gz",
                index = rules.salmon_index.output.idx
        output: ctsdir = directory("COUNTS/Salmon/{file}")
        log:    "LOGS/{file}/salmonquant.log"
        conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
        threads: MAXTHREAD
        params: cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'DTU')['OPTIONS'][1].items()),
                mapp=COUNTBIN,
                stranded = lambda x: '-l SF' if (stranded == 'fr' or stranded == 'SF') else '-l SR' if (stranded == 'rf' or stranded == 'SR') else ''
        shell: "{params.mapp} quant -p {threads} -i {input.index} {params.stranded} {params.cpara} -o {output.ctsdir} -1 {input.r1} 2>> {log} "

rule create_annotation_table:
    input:   dir  = expand(rules.mapping.output.ctsdir, file=samplecond(SAMPLES,config)),
    output:  anno = expand("{outdir}/Tables/ANNOTATION.gz",outdir=outdir)
    log:     expand("LOGS/{outdir}/create_DTU_table.log",outdir=outdir)
    conda:   "nextsnakes/envs/"+COUNTENV+".yaml"
    threads: 1
    params:  dereps = lambda wildcards, input: get_reps(input.dir,config,'DTU'),
             bins = BINS
             # tpara  = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(samplecond(SAMPLES,config)[0], None ,config, "DTU")['OPTIONS'][1].items())
    shell: "python3 {params.bins}/Analysis/build_DTU_table.py {params.dereps} --anno {output.anno} --loglevel DEBUG 2> {log}"

rule run_DTU:
    input:  anno = rules.create_annotation_table.output.anno,
    output: session = rules.themall.input.session,
            res_t = expand("{outdir}/DTU_DRIMSEQ_{comparison}_results_transcripts.tsv.gz", outdir=outdir, comparison=compstr),
            res_g = expand("{outdir}/DTU_DRIMSEQ_{comparison}_results_genes.tsv.gz", outdir=outdir, comparison=compstr),
            res_p = expand("{outdir}/DTU_DRIMSEQ_{comparison}_results_proportions.tsv.gz", outdir=outdir, comparison=compstr)
            # res_stager = rules.themall.input.res_stager,
            # res_posthoc = rules.themall.input.res_posthoc,
    log:    expand("LOGS/{outdir}run_DTU.log",outdir=outdir)
    conda:  "nextsnakes/envs/"+DTUENV+"_DTU.yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins   = str.join(os.sep,[BINS,DTUBIN]),
            compare = comparison,
            outdir = outdir,
            ref = ANNOTATION,
            cutts = get_cutoff_as_string(config, 'DTU')
            # pvcut = lambda wildcards: ' '.join(f"{val}" for (key,val) in tool_params(SAMPLES[0], None ,config, 'DTU')['OPTIONS'][2].items())
    shell: "Rscript --no-environ --no-restore --no-save {params.bins} {input.anno} {params.ref} {params.outdir} {params.compare} {threads} {params.cutts} 2> {log}"

rule filter_significant:
    input:  res_t = rules.run_DTU.output.res_t,
            res_g = rules.run_DTU.output.res_g
    output: sig  = rules.themall.input.sig,
            sig_d  = rules.themall.input.sig_d,
            sig_u  = rules.themall.input.sig_u,
    log:    expand("LOGS/{outdir}filter_drimseq.log",outdir=outdir)
    conda:  "nextsnakes/envs/"+DTUENV+"_DTU.yaml"
    threads: 1
    params: #pv_cut = get_cutoff_as_string(config, 'DTU')[0]['pval'] if get_cutoff_as_string(config, 'DTU')[0]['pval'] else 0.05,
            pv_cut = re.findall("\d+\.\d+", get_cutoff_as_string(config, 'DTU').split("-")[0]),
            lfc_cut = re.findall("\d+\.\d+", get_cutoff_as_string(config, 'DTU').split("-")[1])
    shell: "for i in {input};do fn=\"${{i##*/}}\"; if [[ -s \"$i\" ]]; then zcat $i| grep -v -w 'NA'|perl -F\'\\t\' -wlane ' next if (!$F[1] || !$F[2]);if ($F[1] =~ /adj_pvalue/ || $F[1] < 0.05 && ($F[2] <= -1.5 ||$F[2] >= 1.5)){{print}}' |gzip > {outdir}/Sig_$fn && zcat $i| grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[1] || !$F[2]);if ($F[1] =~ /adj_pvalue/ || $F[1] < 0.05 && ($F[2] >= 1.5)){{print}}' |gzip > {outdir}/SigUP_$fn && zcat $i| grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[1] || !$F[2]);if ($F[1] =~ /adj_pvalue/ || $F[1] < 0.05 && ($F[2] <= -1.5)){{print}}' |gzip > {outdir}/SigDOWN_$fn; else touch {outdir}/Sig_$fn {outdir}/SigUP_$fn {outdir}/SigDOWN_$fn; fi; done 2> {log}"
