#SAMPLES = download_samples(config)
SRABIN, SRAENV = env_bin_from_config3(config,'SRA')

if paired == 'paired':
    log.info('Downloading paired fastq files from SRA')
    rule themall:
        input: expand("FASTQ/{srafile}_{read}.fastq.gz", srafile=SAMPLES, read=['R1','R2'])

    rule get_from_sra:
        output: fq = expand("FASTQ/{{srafile}}_{read}.fastq.gz", read=['R1','R2'])
        log:    "LOGS/SRA/{srafile}.log"
        conda:  ""+SRAENV+".yaml"
        threads: MAXTHREAD
        params: outdir = lambda w, output: expand("{cond}", cond=[os.path.dirname(x) for x in output.fq]),
                ids = lambda w: expand("{accession}", accession = [os.path.basename(x) for x in SAMPLES]),
                spara = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'SRA', SRAENV)['OPTIONS'].get('DOWNLOAD', ""),
        shell:  "arr=({params.ids}); orr=({params.outdir}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do fasterq-dump -O ${{orr[$i]}} -e {threads} -t TMP {params.spara} --split-files ${{arr[$i]}} 2> {log};done && cd ${{orr[$i]}} && rename 's/_1/_R1/' *.fastq && rename 's/_2/_R2/' *.fastq && pigz -p {threads} *.fastq"

else:
    log.info('Downloading single-end fastq files from SRA')
    rule themall:
        input: expand("FASTQ/{srafile}.fastq.gz", srafile=SAMPLES)

    rule get_from_sra:
        output: fq = "FASTQ/{srafile}.fastq.gz"
        log:    "LOGS/SRA/{srafile}.log"
        conda:  ""+SRAENV+".yaml"
        threads: MAXTHREAD
        params: outdir = lambda w, output: expand("{cond}", cond=os.path.dirname(output.fq)),
                ids = lambda w: expand("{accession}", accession = [os.path.basename(x) for x in SAMPLES]),
                spara = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'SRA', SRAENV)['OPTIONS'].get('DOWNLOAD', ""),
        shell: "arr=({params.ids}); orr=({params.outdir}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do fasterq-dump -O ${{orr[$i]}} -e {threads} -t TMP {params.spara} ${{arr[$i]}} 2> {log};done && cd ${{orr[$i]}} && pigz -p {threads} *.fastq"
