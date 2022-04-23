FETCHBIN, FETCHENV = env_bin_from_config3(config,'FETCH')

if paired == 'paired':
    log.info('Downloading paired fastq files from SRA')
    rule themall:
        input: expand("FASTQ/{srafile}_{read}.fastq.gz", srafile=SAMPLES, read=['R1','R2'])
        
    rule fetch_from_sra:
        output: prefetch = temp("TMP/{srafile}.sra")
        log:    "LOGS/FETCH/prefetch_{srafile}.log"
        conda:  ""+FETCHENV+".yaml"
        params: ids = lambda w: expand("{accession}", accession = [os.path.basename(x) for x in SAMPLES]),
                fpara = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'FETCH', FETCHENV)['OPTIONS'].get('PREFETCH', ""),
        shell: "arr=({params.ids}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do prefetch ${{arr[$i]}} -o {output.prefetch} {params.fpara} &> {log};done"

    rule get_from_sra:
        input: prefetch = rules.fetch_from_sra.output.prefetch
        output: fq = expand("FASTQ/{{srafile}}_{read}.fastq.gz", read=['R1','R2'])
        log:    "LOGS/FETCH/{srafile}.log"
        conda:  ""+FETCHENV+".yaml"
        threads: min(MAXTHREAD, 6)
        params: outdir = lambda w, output: expand("{cond}", cond=[os.path.dirname(x) for x in output.fq]),
                ids = lambda w, input: os.path.abspath(input.prefetch),
                spara = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'FETCH', FETCHENV)['OPTIONS'].get('DOWNLOAD', ""),
                wd = lambda wildcards: os.getcwd()+os.sep,
                sra = lambda wildcards: os.path.splitext(os.path.basename(wildcards.srafile))[0]
        shell: "set +euo pipefail; fasterq-dump -O {params.outdir[0]} -e {threads} -t TMP {params.spara} --split-files {params.ids} &> {params.wd}{log} ; cd {params.outdir[0]} ; rename 's/(.sra)*_([1|2])/_R$2/' {params.sra}*.fastq &>> {params.wd}{log} ; for i in {params.sra}*.fastq;do gzip $i &>> {params.wd}{log};done ; exit 0"

else:
    log.info('Downloading single-end fastq files from SRA')
    rule themall:
        input: expand("FASTQ/{srafile}.fastq.gz", srafile=SAMPLES)

    rule fetch_from_sra:
        output: prefetch = temp("TMP/{srafile}.sra")
        log:    "LOGS/FETCH/prefetch_{srafile}.log"
        conda:  ""+FETCHENV+".yaml"
        params: fpara = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'FETCH', FETCHENV)['OPTIONS'].get('PREFETCH', ""),
        shell: "arr=({params.ids}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do prefetch ${{arr[$i]}} -o {output.prefetch} {params.fpara} &> {log};done"

    rule get_from_sra:
        input: prefetch = rules.fetch_from_sra.output.prefetch
        output: fq = "FASTQ/{srafile}.fastq.gz"
        log:    "LOGS/FETCH/{srafile}.log"
        conda:  ""+FETCHENV+".yaml"
        threads: min(MAXTHREAD, 6)
        params: outdir = lambda w, output: expand("{cond}", cond=os.path.dirname(output.fq)),
                ids = lambda w, input: os.path.abspath(input.prefetch),
                spara = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'FETCH', FETCHENV)['OPTIONS'].get('DOWNLOAD', ""),
                wd = lambda wildcards: os.getcwd()+os.sep,
                sra = lambda wildcards:  os.path.splitext(os.path.basename(wildcards.srafile))[0]
        shell: "set +euo pipefail; fasterq-dump -O {params.outdir[0]} -e {threads} -t TMP {params.spara} {params.ids} &> {params.wd}{log} ; cd {params.outdir[0]} ; rename 's/(.sra)*_([1|2])/_R$2/' {params.sra}*.fastq &>> {params.wd}{log} ; for i in {params.sra}*.fastq; do gzip $i &>> {params.wd}{log};done ; exit 0"
