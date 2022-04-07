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
A successful run will show the message 'Workflow finished, no error'


A more complex run
###################

This slightly more complex use case involves multiple input files, two conditions (WT/KO) and a more or less standard DE analysis workflow. We also include a "dummylevel" that is a placeholder for settings or other subdivisions of the WT level, to demonstrate that ``MONSDA`` can work on condition-trees of differing depth. 

Workflows include: 

    - FETCH: Download from SRA
    - QC: FASTQC of input and output
    - TRIMMING: Adaptor removal with cutadapt/trimgalore
    - MAPPING: Read mapping with STAR, hisat2, bwa, segemehl3 and minimap2 
    - DEDUP: Read deduplication with umi_tools and picard
    - DE: Differential Expression analysis with EdgeR and DESeq2

The more complex config for this analysis follows

.. literalinclude:: ../../configs/tutorial_de.json
    :language: json

Note that "SETTINGS" now also contain "GROUPS" which hold identifiers for the DE stage and can be used to define "COMPARABLES" that tell ``MONSDA`` which Samples should be treated as replicates and compared to which other samples. By default this will make an all-vs-all comparison, in our simple case ctrl-vs-ko. Keeping this "GROUPS" setting separated from the condition-tree makes it possible to mix samples for different conditions for all pre- and processing steps and separating them for postprocessing without further struggle. This means that in this simple case we could have mixed all samples under one condition as all other tool options stay the same and reference and annotation do not differ and would still be able to split by "GROUPS" for the DE step. However, sometimes we need to compare data that has to be processed differently, e.g. mixing paired and single-ended reads, so we can apply different settings as needed for all workflow steps by simply splitting our condition-tree accordingly.

Starting the run with 12 cores (defining more will be capped by the config file as MAXTHREADS is set to 12)

.. code-block:: bash

    monsda -j 21 -c ${CONDA_PREFIX}/share/MONSDA/configs/tutorial_de.json --directory ${PWD}

Will start the run in the current directory and generate a "FASTQ" sub-directory containing the downloaded sample, a "GENOME/INDICES" directory containing the built index, a "QC" directory containing all FASTQC reports and MULTIQC output, a "TRIMMED_FASTQ" directory for trimgalore and cutadapt output, a "DEDUP" directory for umi_tools (runs before trimming and after mapping) and picard (runs after mapping) output and a "MAPPING" directory containing the mapped files. Furthermore, a "DE" directory will be created which will hold output from counting with featurecounts and DE input and output from EDGER and DESeq2. Again, ``MONSDA`` will create a "LOG" directory containing it's own log, as well as logs of all executed jobs. 

A successful run will show the message 'Workflow finished, no error'

Run it all
###########

This exhaustive use case involves multiple input files, two conditions (WT/KO) and a set of postprocessing workflows. We also include a "dummylevel" that is a placeholder for settings or other subdivisions of the WT level, to demonstrate that ``MONSDA`` can work on condition-trees of differing depth. 

Workflows include: 

    - FETCH: Download from SRA
    - QC: FASTQC of input and output
    - TRIMMING: Adaptor removal with cutadapt/trimgalore
    - MAPPING: Read mapping with STAR, hisat2, bwa, segemehl3 and minimap2 
    - DEDUP: Read deduplication with umi_tools and picard
    - DE: Differential Expression analysis with EdgeR and DESeq2
    - DEU: Differential Exon Usage analysis with EdgeR and DEXSeq
    - DAS: Differential Alternative Splicing analysis with EdgeR and DIEGO
    - DTU: Differential Transcript Usage analysis with DRIMSeq and DEXSeq
    - COUNTING: Read counting with FeaturCounts ot quantification with Salmon
    - TRACKS: Generation of tracks for UCSC or other genome browsers
    - PEAKS: Analysis of ChIP-Seq or CLIP-Seq or cyPhyRNA-Seq Peaks

The more complex config for this analysis follows


.. literalinclude:: ../../configs/tutorial_exhaustive.json
    :language: json

Note that "SETTINGS" now also contain "GROUPS" which hold identifiers for the DE stage and can be used to define "COMPARABLES" that tell ``MONSDA`` which Samples should be treated as replicates and compared to which other samples. By default this will make an all-vs-all comparison, in our simple case ctrl-vs-ko. Keeping this "GROUPS" setting separated from the condition-tree makes it possible to mix samples for different conditions for all pre- and processing steps and separating them for postprocessing without further struggle. This means that in this simple case we could have mixed all samples under one condition as all other tool options stay the same and reference and annotation do not differ and would still be able to split by "GROUPS" for the DE step. However, sometimes we need to compare data that has to be processed differently, e.g. mixing paired and single-ended reads, so we can apply different settings as needed for all workflow steps by simply splitting our condition-tree accordingly.

Starting the run with 12 cores (defining more will be capped by the config file as MAXTHREADS is set to 12)

.. code-block:: bash

    monsda -j 21 -c ${CONDA_PREFIX}/share/MONSDA/configs/tutorial_de.json --directory ${PWD}

Will start the run in the current directory and generate a "FASTQ" sub-directory containing the downloaded sample, a "GENOME/INDICES" directory containing the built index, a "QC" directory containing all FASTQC reports and MULTIQC output, a "TRIMMED_FASTQ" directory for trimgalore and cutadapt output, a "DEDUP" directory for umi_tools (runs before trimming and after mapping) and picard (runs after mapping) output and a "MAPPING" directory containing the mapped files. Furthermore, a "DE" directory will be created which will hold output from counting with featurecounts and DE input and output from EDGER and DESeq2. Again, ``MONSDA`` will create a "LOG" directory containing it's own log, as well as logs of all executed jobs. 

A successful run will show the message 'Workflow finished, no error'