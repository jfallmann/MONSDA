from datetime import datetime
logid = 'summary.smk '

outdir = "REPORTS"

rule themall:
    # input:  summary_files = expand("{outdir}/rmarkdown_summary.{format}", outdir=outdir, format=config["SUMMARY"]["FORMAT"]),
    #         summary_Rmd = expand("{outdir}/rmarkdown_summary.Rmd", outdir=outdir, analyses="-".join(config["SUMMARY"].keys()))
    input:  summary_allpost = expand("{outdir}/summary_postprocessing.pdf", outdir=outdir),
            summarys = expand("{dir}/summary.pdf", dir=get_summary_dirs(config))

rule make_rmd:
    input:  expand("{files}", files=get_summary_files(config))
    output: expand("{outdir}/summary_postprocessing.Rmd", outdir=outdir)
    log:    expand("LOGS/{outdir}/summary.log", outdir=outdir)
    conda:  "nextsnakes/envs/summary.yaml"
    shell:  "Rscript --no-environ --no-restore --no-save nextsnakes/scripts/Analysis/SUMMARY.Rmd {input} {outdir} 2> {log}"

rule knitr_rmd:
    input:  rules.make_rmd.output
    output: rules.themall.input.summary_allpost,
            rules.themall.input.summarys
    log:    expand("LOGS/{outdir}/summary.log", outdir=outdir)
    conda:  "nextsnakes/envs/summary.yaml"
    shell:  "R rmarkdown::render('{input}',output_file='{output}')"
