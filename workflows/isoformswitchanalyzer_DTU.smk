DTUBIN, DTUENV = env_bin_from_config3(config,'DTU')
COUNTBIN, COUNTENV = ['featureCounts','countreads_de']#env_bin_from_config3(config,'COUNTING') ##PINNING subreads package to version 1.6.4 due to changes in 2.0.1 gene_id length cutoff that interfers

outdir = "DTU/DEXSEQ"
comparison = comparable_as_string2(config,'DTU')
compstr = [i.split(":")[0] for i in comparison.split(",")]

rule themall:
    input:


rule salmon_index:
    input: "assembly/transcriptome.fasta"
    output:
        directory("salmon/transcriptome_index")
    log:
        "logs/salmon/transcriptome_index.log"
    threads: 2
    params:
        # optional parameters
        extra=""
    wrapper:
        "master/bio/salmon/index"

if paired == 'paired':
    rule salmon_quant_reads:
        input:
            # If you have multiple fastq files for a single sample (e.g. technical replicates)
            # use a list for r1 and r2.
            r1 = "reads/{sample}_1.fq.gz",
            r2 = "reads/{sample}_2.fq.gz",
            index = "salmon/transcriptome_index"
        output:
            quant = 'salmon/{sample}/quant.sf',
            lib = 'salmon/{sample}/lib_format_counts.json'
        log:
            'logs/salmon/{sample}.log'
        params:
            # optional parameters
            libtype ="A",
            #zip_ext = bz2 # req'd for bz2 files ('bz2'); optional for gz files('gz')
            extra=""
        threads: 2
        wrapper:
            "master/bio/salmon/quant"
else:
        rule salmon_quant_reads:
            input:
                # If you have multiple fastq files for a single sample (e.g. technical replicates)
                # use a list for r1 and r2.
                r1 = "reads/{sample}_1.fq.gz",
                r2 = "reads/{sample}_2.fq.gz",
                index = "salmon/transcriptome_index"
            output:
                quant = 'salmon/{sample}/quant.sf',
                lib = 'salmon/{sample}/lib_format_counts.json'
            log:
                'logs/salmon/{sample}.log'
            params:
                # optional parameters
                libtype ="A",
                #zip_ext = bz2 # req'd for bz2 files ('bz2'); optional for gz files('gz')
                extra=""
            threads: 2
            wrapper:
                "master/bio/salmon/quant"
