logid = 'dexseq_DTU.smk '
DTUBIN, DTUENV = env_bin_from_config(config,'DTU')
log.debug(logid+"DTUENV: "+str(DTUENV))
COUNTBIN, COUNTENV = ['salmon','salmon'] #env_bin_from_config(config,'COUNTING') ##PINNING subreads package to version 1.6.4 due to changes in 2.0.1 gene_id length cutoff that interfers

comparison = comparable_as_string(config,'DTU')
compstr = [i.split(":")[0] for i in comparison.split(",")]
usededup = config.get('RUNDEDUP', False)

TERMINUSENV = 'terminus'
termpara = tool_params(SAMPLES[0], None, config, 'DTU', DTUENV)['OPTIONS'].get('TERMINUS', None)
runterminus = termpara is not None

keydict = sub_dict(tool_params(SAMPLES[0], None, config, 'DTU', DTUENV)['OPTIONS'], ['INDEX'])
keydict["REF"] = REFERENCE
keydict["DECOY"] = DECOY
keydict["ENV"] = DTUENV
unik = get_dict_hash(keydict)

rule themall:
    input:  session = expand("DTU/{combo}/DTU_DEXSEQ_{scombo}_SESSION.gz", combo=combo, scombo=scombo),
            res     = expand("DTU/{combo}/Tables/DTU_DEXSEQ_{scombo}_{comparison}_table_results.tsv.gz", combo=combo, scombo=scombo, comparison=compstr),
            # sig     = expand("DTU/{combo}/Tables/Sig_DTU_DEXSEQ_{scombo}_{comparison}_results.tsv.gz", combo=combo, scombo=scombo, comparison=compstr),
            # sig_d   = expand("DTU/{combo}/Tables/SigDOWN_DTU_DEXSEQ_{scombo}_{comparison}_results.tsv.gz", combo=combo, scombo=scombo, comparison=compstr),
            # sig_u   = expand("DTU/{combo}/Tables/SigUP_DTU_DEXSEQ_{scombo}_{comparison}_results.tsv.gz", combo=combo, scombo=scombo, comparison=compstr),
            Rmd     = expand("REPORTS/SUMMARY/RmdSnippets/{combo}.Rmd", combo=combo)

include: "dtu_base.smk"

rule run_DTU:
    input:  anno = expand(rules.create_annotation_table.output.anno, combo=combo, scombo=scombo)
    output: session = expand(rules.themall.input.session, combo=combo, scombo=scombo),
            res = rules.themall.input.res
    log:    expand("LOGS/{combo}_{scombo}_{comparison}/DTU/dexseq/run_DTU.log", combo=combo, scombo=scombo, comparison=compstr)
    conda:  ""+DTUENV+".yaml"
    container: "oras://jfallmann/monsda:"+DTUENV+""
    threads: 1  # Due to BPPARAM errors, else int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins   = str.join(os.sep,[BINS, DTUBIN]),
            compare = comparison,
            pcombo = scombo if scombo != '' else 'none',
            outdir = 'DTU/'+combo,
            ref = os.path.abspath(ANNOTATION),
            dtuopt = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'DTU', DTUENV)['OPTIONS'].get('DTU', "")
            # pv_cut = get_cutoff_as_string(config, 'DTU', 'pval'),
            # lfc_cut = get_cutoff_as_string(config, 'DTU', 'lfc')
    shell: "Rscript --no-environ --no-restore --no-save {params.bins} {input.anno} {params.ref} {params.outdir} {params.compare} {params.pcombo} {threads} \'{params.dtuopt}\' 2> {log}"

# rule filter_significant:
#     input:  res = rules.run_DTU.output.res
#     output: sig   = rules.themall.input.sig,
#             sig_d = rules.themall.input.sig_d,
#             sig_u = rules.themall.input.sig_u
#     log:    expand("LOGS/{combo}_{scombo}_{comparison}/DTU/dexseq/filter_drimseq.log", combo=combo, scombo=scombo, comparison=compstr)
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
    log:    expand("LOGS/{combo}/DTU/dexseq/create_summary_snippet.log", combo=combo)
    conda:  ""+DTUENV+".yaml"
    container: "oras://jfallmann/monsda:"+DTUENV+""
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins = BINS,
            abspathfiles = lambda w, input: [os.path.abspath(x) for x in input]
    shell:  "python3 {params.bins}/Analysis/RmdCreator.py --files {params.abspathfiles} --output {output} --env {DTUENV} --loglevel DEBUG 2>> {log}"
