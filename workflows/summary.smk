from datetime import datetime
logid = 'summary.smk '

outdir = "REPORTS/SUMMARY"

rule themall:
    input:  summary_all = expand("{outdir}/SUMMARY.pdf", outdir=outdir)
            # summarys = expand("{dir}/SUMMARY.pdf", dir=get_summary_dirs(config))

rule make_rmd:
    input:  expand("{files}", files=get_summary_files(config))
    output: rules.themall.input.summary_all
            # rules.themall.input.summarys
    log:    expand("LOGS/{outdir}/make_rmd.log", outdir=outdir)
    conda:  "nextsnakes/envs/summary.yaml"
    params: outdir = outdir,
            rmd = os.path.join(outdir,'summary.Rmd'),
            inputstring = lambda wildcards, input: '-'.join(input),
            currentpath = os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"..")
    shell:  "Rscript -e \"rmarkdown::render('{params.rmd}',params=list(files='{params.inputstring}',root='{params.currentpath}/'),output_file='{params.currentpath}/{params.outdir}/SUMMARY.pdf')\" 2> {log}"
