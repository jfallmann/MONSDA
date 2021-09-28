from datetime import datetime
logid = 'summary.smk '

outdir = "REPORTS/SUMMARY"

rule themall:
    input:  summary_all = expand("{outdir}/SUMMARY.html", outdir=outdir)
            # summarys = expand("{dir}/SUMMARY.pdf", dir=get_summary_dirs(config))

rule make_rmd:
    input:  os.path.join(outdir,'summary.Rmd')
    output: rules.themall.input.summary_all
            # rules.themall.input.summarys
    log:    expand("LOGS/{outdir}/make_rmd.log", outdir=outdir)
    conda:  "summary.yaml"
    params: outdir = outdir,
            currentpath = os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"..")
    shell:  "Rscript --vanilla -e \"rmarkdown::render('{input}',params=list(root='{params.currentpath}/'),output_file='{params.currentpath}/{params.outdir}/SUMMARY.html', quiet=TRUE)\" 2> {log}"
