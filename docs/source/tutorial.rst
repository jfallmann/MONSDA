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

All config files this tutorial is based on can be found in ``envs/monsda/share/MONSDA/configs`` of your local ``conda`` installation or the "configs" directory of the ``MONSDA`` repository.

A simple run
#############

This run is based on a simple mapping only use-case. We have a single paired end FASTQ file and want to map that. In this case we have no information about strandedness and ignore this setting.

The corresponding config file looks as follows:

.. literalinclude:: ../../configs/tutorial_quick.json
    :language: json


This tells ``MONSDA`` version "1.0.0" to run the workflows "FETCH" and "MAPPING" on the sample "SRR16324019". 

"FETCH" will activate the environment "sra" and run the binary "fasterq-dump" with no extra options to fetch the sample from SRA.

"MAPPING" will then run the executable "STAR" in environment "star" on the downloaded FASTQ file. First it will check for existence of the needed INDEX file, and as no INDEX file is defined in the config it will try to generate one. For that it will use the REFERENCE "GENOMES/Ecoli/ecoli.fa.gz" as well as the ANNOTATION GTF "GENOMES/Ecoli/ecoli.gtf.gz". GTF is always selected first and GFF is used as a fallback if no GTF is found. 
After successfully building the *STAR* index, ``MONSDA`` will start mapping with *STAR* in "paired" mode and applying a set of options as can be seen under the "MAP" key in section "MAPPING"->"SIMPLE"->"star"->"OPTIONS". We do not define a pre- or suffix for output files, so "EXTENSION" can be skipped or left blank.

Starting the run with 4 cores (defining more will be capped by the config file as MAXTHREADS is set to 4)

.. code-block:: bash

    monsda -j 4 -c ${CONDA_PREFIX}/share/MONSDA/configs/tutorial_quick.json --directory ${PWD}

Will start the run in the current directory and generate a "FASTQ" sub-directory containing the downloaded sample, a "GENOME/INDICES" directory containing the built index and a "MAPPING" directory containing the mapped files. Furthermore, ``MONSDA`` will create a "LOG" directory containing it's own log, as well as logs of all executed jobs.

A more complex run
###################

This slightly more complex use case involves multiple input files, two conditions (WT/KO) and a more or less standard DE analysis workflow.

.. literalinclude:: ../../configs/tutorial_de.json
    :language: json


Run it all
###########

.. literalinclude:: ../../configs/tutorial_exhaustive.json
    :language: json

