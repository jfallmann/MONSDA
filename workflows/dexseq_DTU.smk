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

rule salmon_index:
    input:  fa = REFERENCE
    output: idx = directory(INDEX),
            uidx = directory(expand("{refd}/INDICES/{mape}_{unikey}", refd=REFDIR, mape=COUNTENV, unikey=unik))
    log:    expand("LOGS/{sets}/{cape}.idx.log", sets=SETS, cape=COUNTENV)    
    conda:  ""+COUNTENV+".yaml"
    container: "oras://jfallmann/monsda:"+COUNTENV+""
    threads: MAXTHREAD
    params: mapp = COUNTBIN,
            ipara = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'DTU', DTUENV)['OPTIONS'].get('INDEX', ""),
            decoy = f"-d {os.path.abspath(DECOY)}" if DECOY else '',
            linkidx = lambda wildcards, output: str(os.path.abspath(output.uidx[0]))
    shell:  "set +euo pipefail; {params.mapp} index {params.ipara} {params.decoy} -p {threads} -t {input.fa} -i {output.uidx} &>> {log} && ln -fs {params.linkidx} {output.idx}"


if paired == 'paired':
    rule simulate_trim:
        input:  r1 = lambda wildcards: "FASTQ/{rawfile}_R1.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]) if not usededup else "DEDUP_FASTQ/{combo}/{file}_R1_dedup.fastq.gz",
                r2 = lambda wildcards: "FASTQ/{rawfile}_R2.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]) if not prededup else "DEDUP_FASTQ/{combo}/{file}_R2_dedup.fastq.gz"
        output: r1 = "TRIMMED_FASTQ/{scombo}/{file}_R1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{scombo}/{file}_R2_trimmed.fastq.gz"
        threads: 1
        params: filetolink = lambda w, input: "{r}".format(r=os.path.abspath(input.r1)),
                filetolink2 = lambda w, input: "{r}".format(r=os.path.abspath(input.r2))
        shell:  "ln -s {params.filetolink} {output.r1} && ln -s {params.filetolink2} {output.r2}"

else:
    rule simulate_trim:
        input:  r1 = lambda wildcards: "FASTQ/{rawfile}.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]) if not usededup else "DEDUP_FASTQ/{combo}/{file}_dedup.fastq.gz"
        output: r1 = "TRIMMED_FASTQ/{scombo}/{file}_trimmed.fastq.gz"
        threads: 1
        params: filetolink = lambda w, input: "{r}".format(r=os.path.abspath(input.r1))
        shell:  "ln -s {params.filetolink} {output.r1}"
        

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
        container: "oras://jfallmann/monsda:"+COUNTENV+""
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
        container: "oras://jfallmann/monsda:"+COUNTENV+""
        threads: MAXTHREAD
        params: cpara = lambda wildcards: tool_params(wildcards.file, None, config, 'DTU', DTUENV)['OPTIONS'].get('QUANT', ""),
                mapp=COUNTBIN,
                stranded = lambda x: '-l SF' if (stranded == 'fr' or stranded == 'SF') else '-l SR' if (stranded == 'rf' or stranded == 'SR') else '-l U',
                linksf = lambda wildcards, output: str(os.path.abspath(output.ctsdir))
        shell: "set +euo pipefail; {params.mapp} quant -p {threads} -i {input.index} {params.stranded} {params.cpara} -o {output.ctsdir} -r {input.r1} &>> {log} && gzip {output.ctsdir}/quant.sf ; ln -fs {params.linksf}/quant.sf.gz {output.cnts} &>> {log}"


if runterminus:
    rule terminus_group:
        input:  ctsdir = rules.mapping.output.ctsdir
        output: grp = "DTU/{combo}/terminus/{file}/groups.txt"
        log:    "LOGS/{combo}/{file}/terminus_group.log"
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
        log:    expand("LOGS/DTU/{combo}/terminus_collapse.log", combo=combo)
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
    log:    expand("LOGS/DTU/{combo}/create_DTU_table.log", combo=combo)
    conda:  ""+COUNTENV+".yaml"
    container: "oras://jfallmann/monsda:"+COUNTENV+""
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
    container: "oras://jfallmann/monsda:"+DTUENV+""
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins = BINS,
            abspathfiles = lambda w, input: [os.path.abspath(x) for x in input]
    shell:  "python3 {params.bins}/Analysis/RmdCreator.py --files {params.abspathfiles} --output {output} --env {DTUENV} --loglevel DEBUG 2>> {log}"
