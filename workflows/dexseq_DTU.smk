logid = 'dexseq_DTU.smk '
DTUBIN, DTUENV = env_bin_from_config3(config,'DTU')
COUNTBIN, COUNTENV = ['salmon','salmon'] #env_bin_from_config3(config,'COUNTING') ##PINNING subreads package to version 1.6.4 due to changes in 2.0.1 gene_id length cutoff that interfers

comparison = comparable_as_string(config,'DTU')
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
    output: idx = directory(INDEX),
            uidx = directory(expand("{refd}/INDICES/{mape}_{unikey}", refd=REFDIR, mape=COUNTENV, unikey=get_dict_hash(subDict(tool_params(SAMPLES[0], None, config, 'DTU', DTUENV)['OPTIONS'], ['INDEX']))))
    log:    expand("LOGS/{sets}/{cape}.idx.log", sets=SETS, cape=COUNTENV)    
    conda:  ""+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: mapp = COUNTBIN,
            ipara = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'DTU', DTUENV)['OPTIONS'].get('INDEX', ""),
            linkidx = lambda wildcards, output: str(os.path.abspath(output.uidx[0]))
    shell:  "set +euo pipefail; {params.mapp} index {params.ipara} -p {threads} -t {input.fa} -i {output.uidx} &>> {log} && ln -fs {params.linkidx} {output.idx}"


if paired == 'paired':
    rule mapping:
        input:  r1 = expand("TRIMMED_FASTQ/{scombo}/{{file}}_R1_trimmed.fastq.gz", scombo=scombo),
                r2 = expand("TRIMMED_FASTQ/{scombo}/{{file}}_R2_trimmed.fastq.gz", scombo=scombo),
                index = rules.salmon_index.output.idx,
                uix = rules.salmon_index.output.uidx
        output: cnts = report("DTU/{combo}/salmon/{file}_counts.sf.gz", category="COUNTING"),
                ctsdir = report(directory("DTU/{combo}/salmon/{file}"), category="COUNTING")
        log:    "LOGS/{combo}/{file}/salmonquant.log"
        conda:  ""+COUNTENV+".yaml"
        threads: MAXTHREAD
        params: cpara = lambda wildcards: tool_params(wildcards.file, None, config, 'DTU', DTUENV)['OPTIONS'].get('QUANT', ""),
                mapp=COUNTBIN,
                stranded = lambda x: '-l ISF' if (stranded == 'fr' or stranded == 'ISF') else '-l ISR' if (stranded == 'rf' or stranded == 'ISR') else '-l IU',
                linksf = lambda wildcards, output: str(os.path.abspath(output.ctsdir))
        shell: "set +euo pipefail; {params.mapp} quant -p {threads} -i {input.index} {params.stranded} {params.cpara} -o {output.ctsdir} -1 {input.r1} -2 {input.r2} &>> {log} && gzip {output.ctsdir}/quant.sf && ln -fs {params.linksf}/quant.sf.gz {output.cnts} &>> {log}"

else:
    rule mapping:
        input:  r1 = expand("TRIMMED_FASTQ/{scombo}/{{file}}_trimmed.fastq.gz", scombo=scombo),
                index = rules.salmon_index.output.idx,
                uix = rules.salmon_index.output.uidx
        output: cnts = report("DTU/{combo}/salmon/{file}_counts.sf.gz", category="COUNTING"),
                ctsdir = report(directory("DTU/{combo}/salmon/{file}"), category="COUNTING")
        log:    "LOGS/{combo}/{file}/salmonquant.log"
        conda:  ""+COUNTENV+".yaml"
        threads: MAXTHREAD
        params: cpara = lambda wildcards: tool_params(wildcards.file, None, config, 'DTU', DTUENV)['OPTIONS'].get('QUANT', ""),
                mapp=COUNTBIN,
                stranded = lambda x: '-l SF' if (stranded == 'fr' or stranded == 'SF') else '-l SR' if (stranded == 'rf' or stranded == 'SR') else '-l U',
                linksf = lambda wildcards, output: str(os.path.abspath(output.ctsdir))
        shell: "set +euo pipefail; {params.mapp} quant -p {threads} -i {input.index} {params.stranded} {params.cpara} -o {output.ctsdir} -r {input.r1} &>> {log} && gzip {output.ctsdir}/quant.sf ; ln -fs {params.linksf}/quant.sf.gz {output.cnts} &>> {log}"


rule create_annotation_table:
    input:  dir  = expand(rules.mapping.output.ctsdir, combo=combo, file=samplecond(SAMPLES, config)),
    output: anno = expand("DTU/{combo}/Tables/{scombo}_ANNOTATION.gz", combo=combo, scombo=scombo)
    log:    expand("LOGS/DTU/{combo}/create_DTU_table.log", combo=combo)
    conda:  ""+COUNTENV+".yaml"
    threads: 1
    params: dereps = lambda wildcards, input: get_reps(input.dir, config,'DTU'),
            bins = BINS
    shell:  "python3 {params.bins}/Analysis/build_DTU_table.py {params.dereps} --anno {output.anno} --loglevel DEBUG 2> {log}"

rule run_DTU:
    input:  anno = expand(rules.create_annotation_table.output.anno, combo=combo, scombo=scombo)
    output: session = expand(rules.themall.input.session, combo=combo, scombo=scombo),
            res = rules.themall.input.res
    log:    expand("LOGS/DTU/{combo}_{scombo}_{comparison}/run_DTU.log", combo=combo, scombo=scombo, comparison=compstr)
    conda:  ""+DTUENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins   = str.join(os.sep,[BINS, DTUBIN]),
            compare = comparison,
            pcombo = scombo if scombo != '' else 'none',
            outdir = 'DTU/'+combo,
            ref = os.path.abspath(ANNOTATION)
            # pv_cut = get_cutoff_as_string(config, 'DTU', 'pval'),
            # lfc_cut = get_cutoff_as_string(config, 'DTU', 'lfc')
    shell: "Rscript --no-environ --no-restore --no-save {params.bins} {input.anno} {params.ref} {params.outdir} {params.pcombo} {params.compare} {threads} 2> {log}"

# rule filter_significant:
#     input:  res = rules.run_DTU.output.res
#     output: sig   = rules.themall.input.sig,
#             sig_d = rules.themall.input.sig_d,
#             sig_u = rules.themall.input.sig_u
#     log:    expand("LOGS/DTU/{combo}_{scombo}_{comparison}/filter_drimseq.log", combo=combo, scombo=scombo, comparison=compstr)
#     conda:  ""+DTUENV+".yaml"
#     threads: 1
#     params: pv_cut = get_cutoff_as_string(config, 'DTU', 'pval'),
#             lfc_cut = get_cutoff_as_string(config, 'DTU', 'lfc')
#     shell: "set +o pipefail; for i in {input};do fn=\"${{i##*/}}\"; if [[ -s \"$i\" ]]; then zcat $i| grep -v -w 'NA'|perl -F\'\\t\' -wlane ' next if (!$F[1] || !$F[2]);if ($F[1] =~ /adj_pvalue/ || $F[1] < {params.pv_cut} && ($F[2] <= -{params.lfc_cut} ||$F[2] >= {params.lfc_cut})){{print}}' |gzip > DTU/{combo}/Sig_$fn && zcat $i| grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[1] || !$F[2]);if ($F[1] =~ /adj_pvalue/ || $F[1] < {params.pv_cut} && ($F[2] >= {params.lfc_cut})){{print}}' |gzip > DTU/{combo}/SigUP_$fn && zcat $i| grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[1] || !$F[2]);if ($F[1] =~ /adj_pvalue/ || $F[1] < {params.pv_cut} && ($F[2] <= -{params.lfc_cut})){{print}}' |gzip > DTU/{combo}/Sig_$fn; else touch DTU/{combo}/Sig_$fn DTU/{combo}/SigUP_$fn DTU/{combo}/SigDOWN_$fn; fi; done 2> {log}"

rule create_summary_snippet:
    input:  rules.run_DTU.output.res,
            rules.themall.input.session
            # rules.filter_significant.output.sig,
            # rules.filter_significant.output.sig_d,
            # rules.filter_significant.output.sig_u,
    output: rules.themall.input.Rmd
    log:    expand("LOGS/DTU/{combo}create_summary_snippet.log", combo=combo)
    conda:  ""+DTUENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins = BINS,
            abspathfiles = lambda w, input: [os.path.abspath(x) for x in input]
    shell:  "python3 {params.bins}/Analysis/RmdCreator.py --files {params.abspathfiles} --output {output} --env {DTUENV} --loglevel DEBUG 2>> {log}"
