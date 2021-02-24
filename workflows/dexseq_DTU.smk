logid = 'dexseq_DTU.smk '
DTUBIN, DTUENV = env_bin_from_config3(config,'DTU')
COUNTBIN, COUNTENV = ['salmon','salmon'] #env_bin_from_config3(config,'COUNTING') ##PINNING subreads package to version 1.6.4 due to changes in 2.0.1 gene_id length cutoff that interfers

comparison = comparable_as_string2(config,'DTU')
compstr = [i.split(":")[0] for i in comparison.split(",")]
log.info(logid+"COMPARISON: "+str(comparison))

rule themall:
    input:  session = expand("DTU/{combo}/DTU_DEXSEQ_{scombo}_SESSION.gz", combo=combo, scombo=scombo),
            res     = expand("DTU/{combo}/Tables/DTU_DEXSEQ_{scombo}_{comparison}_table_results.tsv.gz", combo=combo, scombo=scombo, comparison=compstr),
            # sig     = expand("DTU/{combo}/Tables/Sig_DTU_DEXSEQ_{scombo}_{comparison}_results.tsv.gz", combo=combo, scombo=scombo, comparison=compstr),
            # sig_d   = expand("DTU/{combo}/Tables/SigDOWN_DTU_DEXSEQ_{scombo}_{comparison}_results.tsv.gz", combo=combo, scombo=scombo, comparison=compstr),
            # sig_u   = expand("DTU/{combo}/Tables/SigUP_DTU_DEXSEQ_{scombo}_{comparison}_results.tsv.gz", combo=combo, scombo=scombo, comparison=compstr),
            Rmd     = expand("REPORTS/SUMMARY/RmdSnippets/{combo}.Rmd", combo=combo)

rule salmon_index:
    input:  fa = REFERENCE
    output: idx = directory(expand("{refd}/INDICES/{mape}", refd=REFDIR, mape=COUNTENV))
    log:    expand("LOGS/{sets}/salmon.idx.log", sets=SETS)
    conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: mapp = COUNTBIN,
            ipara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(SAMPLES[0], None, config, 'DTU', DTUENV.split("_")[0])['OPTIONS'][0].items()),
    shell:  "{params.mapp} index {params.ipara} -p {threads} -t {input.fa} -i {output.idx} 2>> {log}"

if paired == 'paired':
    rule mapping:
        input:  r1 = "FASTQ/{file}_R1.fastq.gz",
                r2 = "FASTQ/{file}_R2.fastq.gz",
                index = rules.salmon_index.output.idx
        output: ctsdir = directory("COUNTS/Salmon/{file}")
        log:    expand("LOGS/DTU/{combo}/salmonquant.log", combo=combo)
        conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
        threads: MAXTHREAD
        params: cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None , config, 'DTU', DTUENV.split("_")[0])['OPTIONS'][1].items()),
                mapp=COUNTBIN,
                stranded = lambda x: '-l ISF' if (stranded == 'fr' or stranded == 'ISF') else '-l ISR' if (stranded == 'rf' or stranded == 'ISR') else ''
        shell: "{params.mapp} quant -p {threads} -i {input.index} {params.stranded} {params.cpara} -o {output.ctsdir} -1 {input.r1} -2 {input.r2} 2>> {log}"

else:
    rule mapping:
        input:  r1 = "FASTQ/{file}.fastq.gz",
                index = rules.salmon_index.output.idx
        output: ctsdir = directory("COUNTS/Salmon/{file}")
        log:    expand("LOGS/DTU/{combo}/salmonquant.log", combo=combo)
        conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
        threads: MAXTHREAD
        params: cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None , config, 'DTU', DTUENV.split("_")[0])['OPTIONS'][1].items()),
                mapp=COUNTBIN,
                stranded = lambda x: '-l SF' if (stranded == 'fr' or stranded == 'SF') else '-l SR' if (stranded == 'rf' or stranded == 'SR') else ''
        shell: "{params.mapp} quant -p {threads} -i {input.index} {params.stranded} {params.cpara} -o {output.ctsdir} -1 {input.r1} 2>> {log} "

rule create_annotation_table:
    input:  dir  = expand(rules.mapping.output.ctsdir, file=samplecond(SAMPLES, config)),
    output: anno = expand("DTU/{combo}/Tables/{scombo}_ANNOTATION.gz", combo=combo, scombo=scombo)
    log:    expand("LOGS/DTU/{combo}/create_DTU_table.log", combo=combo)
    conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
    threads: 1
    params: dereps = lambda wildcards, input: get_reps(input.dir, config,'DTU'),
            bins = BINS
    shell:  "python3 {params.bins}/Analysis/build_DTU_table.py {params.dereps} --anno {output.anno} --loglevel DEBUG 2> {log}"

rule run_DTU:
    input:  anno = rules.create_annotation_table.output.anno,
    output: session = rules.themall.input.session,
            res = rules.themall.input.res
    log:    "LOGS/DTU/{combo}_{scombo}_{comparison}/run_DTU.log"
    conda:  "nextsnakes/envs/"+DTUENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins   = str.join(os.sep,[BINS, DTUBIN]),
            compare = comparison,
            scombo = scombo,
            outdir = 'DTU/'+combo,
            ref = ANNOTATION,
            # pv_cut = get_cutoff_as_string(config, 'DTU', 'pval'),
            # lfc_cut = get_cutoff_as_string(config, 'DTU', 'lfc')
    shell: "Rscript --no-environ --no-restore --no-save {params.bins} {input.anno} {params.ref} {params.outdir} {params.scombo} {params.compare} {threads} 2> {log}"

rule filter_significant:
    input:  res = rules.run_DTU.output.res
    output: sig   = rules.themall.input.sig,
            sig_d = rules.themall.input.sig_d,
            sig_u = rules.themall.input.sig_u
    log:    expand("LOGS/DTU/{combo}_{scombo}_{comparison}/filter_drimseq.log", combo=combo, scombo=scombo, comparison=compstr)
    conda:  "nextsnakes/envs/"+DTUENV+".yaml"
    threads: 1
    params: pv_cut = get_cutoff_as_string(config, 'DTU', 'pval'),
            lfc_cut = get_cutoff_as_string(config, 'DTU', 'lfc')
    shell: "set +o pipefail; for i in {input};do fn=\"${{i##*/}}\"; if [[ -s \"$i\" ]]; then zcat $i| grep -v -w 'NA'|perl -F\'\\t\' -wlane ' next if (!$F[1] || !$F[2]);if ($F[1] =~ /adj_pvalue/ || $F[1] < {params.pv_cut} && ($F[2] <= -{params.lfc_cut} ||$F[2] >= {params.lfc_cut})){{print}}' |gzip > DTU/{combo}/Sig_$fn && zcat $i| grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[1] || !$F[2]);if ($F[1] =~ /adj_pvalue/ || $F[1] < {params.pv_cut} && ($F[2] >= {params.lfc_cut})){{print}}' |gzip > DTU/{combo}/SigUP_$fn && zcat $i| grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[1] || !$F[2]);if ($F[1] =~ /adj_pvalue/ || $F[1] < {params.pv_cut} && ($F[2] <= -{params.lfc_cut})){{print}}' |gzip > DTU/{combo}/Sig_$fn; else touch DTU/{combo}/Sig_$fn DTU/{combo}/SigUP_$fn DTU/{combo}/SigDOWN_$fn; fi; done 2> {log}"

rule create_summary_snippet:
    input:  rules.run_DTU.output.res,
            # rules.filter_significant.output.sig,
            # rules.filter_significant.output.sig_d,
            # rules.filter_significant.output.sig_u,
    output: rules.themall.input.Rmd
    log:    expand("LOGS/DTU/{combo}create_summary_snippet.log", combo=combo)
    conda:  "nextsnakes/envs/"+DTUENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins = BINS
    shell:  "python3 {params.bins}/Analysis/RmdCreator.py --files {input} --output {output} --loglevel DEBUG 2> {log}"
