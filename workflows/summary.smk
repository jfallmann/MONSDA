from datetime import datetime
logid = 'summary.smk '

outdir = "SUMMARY"

rule themall:
    # input:  summary_files = expand("{outdir}/rmarkdown_summary.{format}", outdir=outdir, format=config["SUMMARY"]["FORMAT"]),
    #         summary_Rmd = expand("{outdir}/rmarkdown_summary.Rmd", outdir=outdir, analyses="-".join(config["SUMMARY"].keys()))
    input:  summary_files = expand("{outdir}/rmarkdown_summary.txt", outdir=outdir)

rule make_rmd:
    input:  expand("{dir}", dir=get_summary_dirs(config))
    output: files = rules.themall.input.summary_files,
    log:    expand("LOGS/{outdir}/summary.log", outdir=outdir)
    conda:  "nextsnakes/envs/summary.yaml"
    params: bins = os.path.join(BINS,config["SUMMARY"]["BIN"]),
            formats =   '+'.join(config["SUMMARY"]["FORMAT"]),
            cutoff = get_summary_cutoff(config),
            dirs = '+'.join(get_summary_dirs(config))
    shell:  "Rscript --no-environ --no-restore --no-save {params.bins} {params.dirs} {params.formats} {outdir} {params.cutoff} 2> {log}"

rule knitr_rmd:
    input :
