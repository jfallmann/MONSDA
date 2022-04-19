logid = 'drimseq_DTU.smk '
DTUBIN, DTUENV = env_bin_from_config3(config,'DTU')
log.info(logid+"DTUENV: "+str(DTUENV))
COUNTBIN, COUNTENV = ['salmon','salmon']#env_bin_from_config2(SAMPLES, config,'COUNTING') ##PINNING subreads package to version 1.6.4 due to changes in 2.0.1 gene_id length cutoff that interfers

comparison = comparable_as_string(config,'DTU')
compstr = [i.split(":")[0] for i in comparison.split(",")]
log.info(logid+"COMPARISON: "+str(comparison))

rule themall:
    input:  session = expand("DTU/{combo}/DTU_DRIMSEQ_{scombo}_SESSION.gz", combo=combo, scombo=scombo),
            res_t   = expand("DTU/{combo}/Tables/DTU_DRIMSEQ_{scombo}_{comparison}_table_transcripts.tsv.gz", combo=combo, comparison=compstr, scombo=scombo),
            res_g   = expand("DTU/{combo}/Tables/DTU_DRIMSEQ_{scombo}_{comparison}_table_genes.tsv.gz", combo=combo, comparison=compstr, scombo=scombo),
            res_p   = expand("DTU/{combo}/Tables/DTU_DRIMSEQ_{scombo}_{comparison}_table_proportions.tsv.gz", combo=combo, comparison=compstr, scombo=scombo),
            res_gwp = expand("DTU/{combo}/Tables/DTU_DRIMSEQ_{scombo}_{comparison}_table_genewise-precision.tsv.gz", combo=combo, comparison=compstr, scombo=scombo),
            fig_F   = expand("DTU/{combo}/Figures/DTU_DRIMSEQ_{scombo}_{comparison}_figure_FeatPerGene.png", combo=combo, comparison=compstr, scombo=scombo),
            fig_P   = expand("DTU/{combo}/Figures/DTU_DRIMSEQ_{scombo}_{comparison}_figure_Precision.png", combo=combo, comparison=compstr, scombo=scombo),
            fig_PV  = expand("DTU/{combo}/Figures/DTU_DRIMSEQ_{scombo}_{comparison}_figure_PValues.png", combo=combo, comparison=compstr, scombo=scombo),
            fig_files = expand("DTU/{combo}/Figures/DTU_DRIMSEQ_{scombo}_{comparison}_list_sigGenesFigures.tsv", combo=combo, comparison=compstr, scombo=scombo),
            sig_t    = expand("DTU/{combo}/Tables/Sig_DTU_DRIMSEQ_{scombo}_{comparison}_table_transcripts.tsv.gz", combo=combo, comparison=compstr, scombo=scombo),
            sig_dt   = expand("DTU/{combo}/Tables/SigDOWN_DTU_DRIMSEQ_{scombo}_{comparison}_table_transcripts.tsv.gz", combo=combo, comparison=compstr, scombo=scombo),
            sig_ut   = expand("DTU/{combo}/Tables/SigUP_DTU_DRIMSEQ_{scombo}_{comparison}_table_transcripts.tsv.gz", combo=combo, comparison=compstr, scombo=scombo),
            sig_g    = expand("DTU/{combo}/Tables/Sig_DTU_DRIMSEQ_{scombo}_{comparison}_table_genes.tsv.gz", combo=combo, comparison=compstr, scombo=scombo),
            sig_dg   = expand("DTU/{combo}/Tables/SigDOWN_DTU_DRIMSEQ_{scombo}_{comparison}_table_genes.tsv.gz", combo=combo, comparison=compstr, scombo=scombo),
            sig_ug   = expand("DTU/{combo}/Tables/SigUP_DTU_DRIMSEQ_{scombo}_{comparison}_table_genes.tsv.gz", combo=combo, comparison=compstr, scombo=scombo),
            # res_stager = expand("DTU/{combo}/DTU_DRIMSEQ_{comparison}_results_stageR-filtered.tsv.gz", combo=combo, comparison=compstr),
            # res_posthoc = expand("DTU/{combo}/DTU_DRIMSEQ_{comparison}_results_post-hoc-filtered-on-SD.tsv.gz", combo=combo, comparison=compstr),
            Rmd = expand("REPORTS/SUMMARY/RmdSnippets/{combo}.Rmd", combo=combo)

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
    params: dereps = lambda wildcards, input: get_reps(input.dir, config, 'DTU'),
            bins = BINS
    shell:  "python3 {params.bins}/Analysis/build_DTU_table.py {params.dereps} --anno {output.anno} --loglevel DEBUG 2> {log}"

rule run_DTU:
    input:  anno = expand(rules.create_annotation_table.output.anno, combo=combo, scombo=scombo)
    output: session = rules.themall.input.session,
            res_t   = rules.themall.input.res_t,
            res_g   = rules.themall.input.res_g,
            res_p   = rules.themall.input.res_p,
            res_gwp = rules.themall.input.res_gwp,
            fig_F   = rules.themall.input.fig_F,
            fig_P   = rules.themall.input.fig_P,
            fig_PV  = rules.themall.input.fig_PV,
            fig_files = rules.themall.input.fig_files
    log:    expand("LOGS/DTU/{combo}/run_DTU.log", combo=combo)
    conda:  ""+DTUENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins   = str.join(os.sep,[BINS, DTUBIN]),
            compare = comparison,
            pcombo = scombo if scombo != '' else 'none',
            outdir = 'DTU/'+combo,
            ref = os.path.abspath(ANNOTATION)
    shell: "Rscript --no-environ --no-restore --no-save {params.bins} {input.anno} {params.ref} {params.outdir} {params.pcombo} {params.compare} {threads} 2> {log}"

rule filter_significant_drimseq:
    input:  res_g = rules.themall.input.res_g,
            res_t = rules.themall.input.res_t
    output: sig_g   = rules.themall.input.sig_g,
            sig_dg  = rules.themall.input.sig_dg,
            sig_ug  = rules.themall.input.sig_ug,
            sig_t   = rules.themall.input.sig_t,
            sig_dt  = rules.themall.input.sig_dt,
            sig_ut  = rules.themall.input.sig_ut
    log:    "LOGS/DTU/filter_drimseqDTU.log"
    conda:  ""+DTUENV+".yaml"
    threads: 1
    params: pv_cut = get_cutoff_as_string(config, 'DTU', 'pvalue'),
            lfc_cut = get_cutoff_as_string(config, 'DTU', 'lfc')
    shell:  "set +o pipefail; arr=({input.res_g}); arrt=({input.res_t}); orr=({output.sig_g}); orrd=({output.sig_dg}); orru=({output.sig_ug}); orrt=({output.sig_t}); orrtd=({output.sig_dt}); orrtu=({output.sig_ut}); for i in \"${{!arr[@]}}\"; do a=\"${{arr[$i]}}\"; fn=\"${{a##*/}}\"; if [[ -s \"$a\" ]];then zcat $a| head -n1 |gzip > \"${{orr[$i]}}\"; cp \"${{orr[$i]}}\" \"${{orrd[$i]}}\"; cp \"${{orr[$i]}}\" \"${{orru[$i]}}\"; zcat $a| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[1]);if ($F[1] < {params.pv_cut} && ($F[2] <= -{params.lfc_cut} ||$F[2] >= {params.lfc_cut}) ){{print}}' |gzip >> \"${{orr[$i]}}\" && zcat $a| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[1]);if ($F[1] < {params.pv_cut} && ($F[2] >= {params.lfc_cut}) ){{print}}' |gzip >> \"${{orru[$i]}}\" && zcat $a| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[1]);if ($F[1] < {params.pv_cut} && ($F[2] <= -{params.lfc_cut}) ){{print}}' |gzip >> \"${{orrd[$i]}}\"; else touch \"${{orr[$i]}}\" \"${{orrd[$i]}}\" \"${{orru[$i]}}\"; fi;done; for i in \"${{!arrt[@]}}\"; do a=\"${{arrt[$i]}}\"; fn=\"${{a##*/}}\"; if [[ -s \"$a\" ]];then zcat $a| head -n1 |gzip > \"${{orrt[$i]}}\"; cp \"${{orrt[$i]}}\" \"${{orrtd[$i]}}\"; cp \"${{orr[$i]}}\" \"${{orrtu[$i]}}\"; zcat $a| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[1]);if ($F[1] < {params.pv_cut} && ($F[2] <= -{params.lfc_cut} ||$F[2] >= {params.lfc_cut}) ){{print}}' |gzip >> \"${{orrt[$i]}}\" && zcat $a| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[1]);if ($F[1] < {params.pv_cut} && ($F[2] >= {params.lfc_cut}) ){{print}}' |gzip >> \"${{orrtu[$i]}}\" && zcat $a| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[1]);if ($F[1] < {params.pv_cut} && ($F[2] <= -{params.lfc_cut}) ){{print}}' |gzip >> \"${{orrtd[$i]}}\"; else touch \"${{orrt[$i]}}\" \"${{orrtd[$i]}}\" \"${{orrtu[$i]}}\"; fi;done 2> {log}"

rule create_summary_snippet:
    input:  rules.run_DTU.output.res_t,
            rules.run_DTU.output.res_g,
            rules.run_DTU.output.res_p,
            rules.run_DTU.output.fig_F,
            rules.run_DTU.output.fig_P,
            rules.run_DTU.output.fig_PV,
            rules.run_DTU.output.fig_files,
            rules.themall.input.sig_g,
            rules.themall.input.sig_dg,
            rules.themall.input.sig_ug,
            rules.themall.input.sig_t,
            rules.themall.input.sig_dt,
            rules.themall.input.sig_ut,
            rules.themall.input.session
    output: rules.themall.input.Rmd
    log:    expand("LOGS/DTU/{combo}/create_summary_snippet.log", combo=combo)
    conda:  ""+DTUENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins = BINS,
            abspathfiles = lambda w, input: [os.path.abspath(x) for x in input]
    shell:  "python3 {params.bins}/Analysis/RmdCreator.py --files {params.abspathfiles} --output {output} --env {DTUENV} --loglevel DEBUG 2>> {log}"
