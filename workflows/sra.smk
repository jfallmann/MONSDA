RAWBIN, RAWENV = env_bin_from_config2(SAMPLES,config,'RAW')

if paired == 'paired':
    log.info('Downloading paired fastq files from SRA')
    rule themall:
        input: expand("FASTQ/{rawfile}_{read}.fastq.gz", rawfile=SAMPLES, read=['R1','R2'])

    rule get_from_sra:
        output: fq = "FASTQ/{rawfile}_{read}.fastq.gz"
        log:    "LOGS/RAW/{rawfile}_{read}.log"
        conda:  "snakes/envs/"+RAWENV+".yaml"
        threads: MAXTHREAD
        params: outdir = lambda w: expand("FASTQ/{cond}",cond=os.path.dirname(w.rawfile)),
                ids = lambda w: expand("{accession}",accession = [os.path.basename(x) for x in download_samples(config)])
        shell:  "for i in {params.ids};do fasterq-dump -O {params.outdir} -e {threads} -t TMP --split-files $i 2> {log};done && cd {params.outdir} && rename _1 _R1 *.fastq && rename _2 _R2 *.fastq && pigz -p {threads} *.fastq"

else:
    log.info('Downloading single-end fastq files from SRA')
    rule themall:
        input: expand("FASTQ/{rawfile}.fastq.gz", rawfile=SAMPLES)

    rule get_from_sra:
        output: fq = "FASTQ/{rawfile}.fastq.gz"
        log:    "LOGS/RAW/{rawfile}.log"
        conda:  "snakes/envs/"+RAWENV+".yaml"
        threads: MAXTHREAD
        params: outdir = lambda w: expand("FASTQ/{cond}",cond=os.path.dirname(w.rawfile)),
                ids = lambda w: expand("{accession}",accession = [os.path.basename(x) for x in download_samples(config)])
        shell:  "for i in {params.ids};do fasterq-dump -O {params.outdir} -e {threads} -t TMP --concatenate-reads {params.ids} 2> {log};done && cd {params.outdir} && pigz -p {threads} *.fastq"
