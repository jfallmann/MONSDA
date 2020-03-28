RAWBIN, RAWENV = env_bin_from_config2(SAMPLES,config,'RAW')

if paired == 'paired':
    log.info('Downloading paired fastq files from SRA')
    rule themall:
        input: expand("FASTQ/{rawfile}_{read}.fastq.gz", rawfile=SAMPLES, read=['R1','R2'])

    rule get_from_sra:
        output: fq = expand("FASTQ/{{rawfile}}_{read}.fastq.gz", read=['R1','R2'])
        log:    "LOGS/RAW/{rawfile}.log"
        conda:  "snakes/envs/"+RAWENV+".yaml"
        threads: MAXTHREAD
        params: outdir = lambda w, output: expand("{cond}",cond=[os.path.dirname(x) for x in output.fq]),
                ids = lambda w: expand("{accession}",accession = [os.path.basename(x) for x in download_samples(config)])
        shell:  "arr=({params.ids}); orr=({params.outdir}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do fasterq-dump -O ${{orr[$i]}} -e {threads} -t TMP --split-files ${{arr[$i]}} 2> {log};done && cd ${{orr[$i]}} && rename _1 _R1 *.fastq && rename _2 _R2 *.fastq && pigz -p {threads} *.fastq"

else:
    log.info('Downloading single-end fastq files from SRA')
    rule themall:
        input: expand("FASTQ/{rawfile}.fastq.gz", rawfile=SAMPLES)

    rule get_from_sra:
        output: fq = "FASTQ/{rawfile}.fastq.gz"
        log:    "LOGS/RAW/{rawfile}.log"
        conda:  "snakes/envs/"+RAWENV+".yaml"
        threads: MAXTHREAD
        params: outdir = lambda w, output: expand("{cond}",cond=os.path.dirname(output.fq)),
                ids = lambda w: expand("{accession}",accession = [os.path.basename(x) for x in download_samples(config)])
        shell: "arr=({params.ids}); orr=({params.outdir}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do fasterq-dump -O ${{orr[$i]}} -e {threads} -t TMP ${{arr[$i]}} 2> {log};done && cd ${{orr[$i]}} && pigz -p {threads} *.fastq"
