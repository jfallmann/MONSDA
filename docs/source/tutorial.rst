TUTORIAL
========

Here we will show three ``MONSDA`` runs from simple to most-complex, explaining the configuration and how to run this in real life.

First please create a working directory and access it, then install ``MONSDA`` following :ref:`_install` and activate the conda environment.

All examples are based on data from SRA_ , which you can either download beforehand or you use the "FETCH" workflow as demonstrated in the tutorial.

.. _SRA: https://www.ncbi.nlm.nih.gov/sra/SRX12601836[accn]

Genome and Annotation files need to be downloaded beforehand and can be found here: GENOME_, GFF_, GTF_

.. _GENOME: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
.. _GFF: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz
.. _GTF: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gtf.gz

Ideally you download those files into a sub-directory of your workspace named "GENOME/Ecoli" and rename them to something simpler, in the tutorial we will refer to "GCF_000005845.2_ASM584v2_genomic" as "ecoli". You also need to change the ".fna.gz" suffix of the genome file to ".fa.gz" and we should be ready to go.

All config files this tutorial is based on can be found in the "configs" directory of the ``MONSDA`` repository.

A simple run
#############

This run is based on a simple mapping only use-case. We have a single paired end FASTQ file and want to map that. In this case we have no information about strandedness and ignore this setting.

.. literalinclude:: ../../configs/tutorial_quick.json
    :language: json


A more complex run
###################

This slightly more complex use case involves multiple input files, two conditions (WT/KO) and a more or less standard DE analysis workflow.

.. literalinclude:: ../../configs/tutorial_de.json
    :language: json


Run it all
###########

.. literalinclude:: ../../configs/tutorial_exhaustive.json
    :language: json

