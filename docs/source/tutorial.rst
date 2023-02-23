.. _tutorials:

TUTORIAL
========

Here we will show three **MONSDA** runs from simple to most-complex, explaining the configuration and how to run this in real life.

First please create a working directory and access it, then install **MONSDA** following :ref:`install` and activate the conda environment.

All examples are based on data from SRA_ , which you can either download beforehand or you use the "FETCH" workflow as demonstrated in the tutorial.

.. _SRA: https://www.ncbi.nlm.nih.gov/sra/SRX12601836[accn]

Genome and Annotation files need to be downloaded beforehand and can be found here: GENOME_, GFF_, GTF_

.. _GENOME: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
.. _GFF: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz
.. _GTF: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gtf.gz

Ideally you download those files into a sub-directory of your workspace named "GENOME/Ecoli" and rename them to something simpler, in the tutorial we will refer to "GCF_000005845.2_ASM584v2_genomic" as "ecoli". You also need to change the ".fna.gz" suffix of the genome file to ".fa.gz" and we should be ready to go.

All config files this tutorial is based on can be found in **envs/monsda/share/MONSDA/configs** of your local **conda** installation or the "configs" directory of the **MONSDA** repository.

Please be aware that we are currently working on the implementation of postprocessing workflows for Nextflow, so at the moment you can only run the *quick* and *toolmix* tutorial in Nextflow mode, but this will change soon.

A quick run
#############

This run is based on a simple mapping only use-case. We have a single paired-end FASTQ file and want to map that. In this case we have no information about strandedness and ignore this setting.

The corresponding config file looks as follows:

.. literalinclude:: ../../configs/tutorial_quick.json
    :language: json


This tells **MONSDA** version "1.0.0" to run the workflows "FETCH" and "MAPPING" on the sample "SRR16324019". 

"FETCH" will activate the environment "sra" and run the binary "fasterq-dump" with no extra options to fetch the sample from SRA.

"MAPPING" will then run the executable "STAR" in environment "star" on the downloaded FASTQ file. First it will check for existence of the needed INDEX file, and as no INDEX file is defined in the config it will try to generate one. For that it will use the REFERENCE "GENOMES/Ecoli/ecoli.fa.gz" as well as the ANNOTATION GTF "GENOMES/Ecoli/ecoli.gtf.gz". GTF is always selected first and GFF is used as a fallback if no GTF is found. 
After successfully building the *STAR* index, **MONSDA** will start mapping with *STAR* in "paired" mode and applying a set of options as can be seen under the "MAP" key in section "MAPPING"->"SIMPLE"->"star"->"OPTIONS". We do not define a pre- or suffix for output files, so "EXTENSION" can be skipped or left blank.

Starting the run with 4 cores (defining more will be capped by the config file as MAXTHREADS is set to 4)

.. code-block:: bash

    monsda -j 4 -c ${CONDA_PREFIX}/share/MONSDA/configs/tutorial_quick.json --directory ${PWD}

Will start the run in the current directory and generate a "FASTQ" sub-directory containing the downloaded sample, a "GENOME/Ecoli/INDICES" directory containing the built index and a "MAPPING" directory containing the mapped files. Furthermore, **MONSDA** will create a "LOG" directory containing it's own log, as well as logs of all executed jobs and a "JOBS" directory containing the 'monsda.commands' file, which contains all command-line calls that **MONSDA** runs.
A successful run will show the message 'Workflow finished, no error' at the end.


A more complex toolmix run
###########################

This slightly more complex use case involves multiple input files, two conditions (WT/KO) and a more or less standard DE analysis workflow. We also include a "dummylevel" that is a placeholder for settings or other subdivisions of the WT level, to demonstrate that **MONSDA** can work on condition-trees of differing depth. 

Workflows include: 

    - FETCH: Download from SRA
    - QC: FASTQC of input and output
    - TRIMMING: Adaptor removal with cutadapt/trimgalore
    - MAPPING: Read mapping with STAR, hisat2, bwa, segemehl3 and minimap2 
    - DEDUP: Read deduplication with umi_tools and picard

The more complex config for this analysis follows

.. literalinclude:: ../../configs/tutorial_toolmix.json
    :language: json


Starting the run with 12 cores (defining more will be capped by the config file as MAXTHREADS is set to 12)

.. code-block:: bash

    monsda -j 12 -c ${CONDA_PREFIX}/share/MONSDA/configs/tutorial_toolmix.json --directory ${PWD}

Will start the run in the current directory and generate a "FASTQ" sub-directory containing the downloaded sample, a "GENOME/INDICES" directory containing the built index, a "QC" directory containing all FASTQC reports and MULTIQC output, a "TRIMMED_FASTQ" directory for trimgalore and cutadapt output, a "DEDUP" directory for umi_tools (runs before trimming and after mapping) and picard (runs after mapping) output and a "MAPPING" directory containing the mapped files. Furthermore, a "DE" directory will be created which will hold output from counting with featurecounts and DE input and output from EDGER and DESeq2. Again, **MONSDA** will create a "LOG" directory containing it's own log, as well as logs of all executed jobs and the "JOBS" directory with all command-line calls. 

A successful run will show the message 'Workflow finished, no error' at the end.

Run postprocessing
###################

This postprocessing use case involves multiple input files, two conditions (WT/KO) and almost all workflows available. We also include a "dummylevel" that is a placeholder for settings or other subdivisions of the WT level, to demonstrate that **MONSDA** can work on condition-trees of differing depth. 

Workflows include: 

    - FETCH: Download from SRA
    - QC: FASTQC of input and output
    - TRIMMING: Adaptor removal with cutadapt/trimgalore
    - MAPPING: Read mapping with STAR, hisat2, bwa, segemehl3 and minimap2 
    - DEDUP: Read deduplication with umi_tools and picard
    - DE: Differential Expression analysis with EdgeR and DESeq2
    - DEU: Differential Exon Usage analysis with EdgeR (DEXSeq skipped, runtime)
    - COUNTING: Read counting with FeatureCounts
    - TRACKS: Generation of tracks for UCSC or other genome browsers
    - PEAKS: Analysis of ChIP-Seq or CLIP-Seq or cyPhyRNA-Seq Peaks

The more complex config for this analysis follows

Note that "SETTINGS" now also contain "GROUPS" which hold identifiers for the DE stage and can be used to define "COMPARABLES" that tell **MONSDA** which samples should be treated as replicates and compared to which other samples. By default this will make an all-vs-all comparison, in our simple case ctrl-vs-ko. Keeping this "GROUPS" setting separated from the condition-tree makes it possible to mix samples for different conditions for all pre- and processing steps and separating them for postprocessing without further struggle. This means that in this simple case we could have mixed all samples under one condition as all other tool options stay the same and reference and annotation do not differ and would still be able to split by "GROUPS" for the DE step. However, sometimes we need to compare data that has to be processed differently, e.g. mixing paired and single-ended reads, so we can apply different settings as needed for all workflow steps by simply splitting our condition-tree accordingly.

Same as for the more complex run we can see that "SETTINGS" contains "GROUPS" which hold identifiers for the DE/DEU/DAS/DTU stages that will be used to define "COMPARABLES" that tell **MONSDA** which samples should be treated as replicates and compared to which other samples. By default this will make an all-vs-all comparison, in our case ctrl-vs-ko. Keeping this "GROUPS" setting separated from the condition-tree makes it possible to mix samples for different conditions for all pre- and processing steps and separating them for postprocessing without further struggle. This means that in this simple case we could have mixed all samples under one condition, as all other tool options stay the same and reference and annotation do not differ and would still be able to split by "GROUPS" for the DE/DEU/DAS steps. However, sometimes we need to compare data that has to be processed differently, e.g. mixing paired and single-ended reads, so we can apply different settings as needed for all workflow steps by simply splitting our condition-tree accordingly. Another addition to the "SETTINGS" key is the "IP" (ImmunoPrecipitation, aka enrichment) part. This is intended to be combined with the *PEAKS* workflow and tools that have certain preprocessing step requirements dependent on the way reads where enriched. For example our novel "scyphy" (Software for cyPhyRNA-Seq) workflow can work on 5-prime and 3-prime end enriched data, which can be indicated as "iCLIP" and "revCLIP" respectively. For this tutorial we will simulate standard "CLIP" and set the corresponding key.

Starting the run with 12 cores (defining more will be capped by the config file as MAXTHREADS is set to 12)

.. code-block:: bash

    monsda -j 12 -c ${CONDA_PREFIX}/share/MONSDA/configs/tutorial_postprocess.json --directory ${PWD}

Will start the run in the current directory and generate a "FASTQ" sub-directory containing the downloaded sample, a "GENOME/Ecoli/INDICES" directory containing the built indices, including the one built for salmon later on, a "QC" directory containing all FASTQC reports and MULTIQC output, a "TRIMMED_FASTQ" directory for trimgalore and cutadapt output, a "DEDUP" directory for umi_tools (runs before trimming and after mapping) and picard (runs after mapping) output and a "MAPPING" directory containing the mapped files. Furthermore, "DE/DEU" directories will be created which will hold output from counting with FeatureCounts and DE/DEU input and output from EDGER and DESeq2 respectively. Again, **MONSDA** will create a "LOG" directory containing it's own log, as well as logs of all executed jobs and again a "JOBS" directory for command-line calls. 

A successful run will show the message 'Workflow finished, no error'. Be aware that this is indeed an exhaustive workflow and will require a decent amount of disk-space, memory and compute-time, depending on the hardware at your disposal.


Coming soon
############

We are working on implementing this exhaustive tutorial, but running analysis steps like DTU or DAS requires more complex datasets and thus more time to test.

This exhaustive use case involves multiple input files, two conditions (WT/KO) and almost all workflows available. We also include a "dummylevel" that is a placeholder for settings or other subdivisions of the WT level, to demonstrate that **MONSDA** can work on condition-trees of differing depth. 

Workflows include: 

    - FETCH: Download from SRA
    - QC: FASTQC of input and output
    - TRIMMING: Adaptor removal with cutadapt/trimgalore
    - MAPPING: Read mapping with STAR, hisat2, bwa, segemehl3 and minimap2 
    - DEDUP: Read deduplication with umi_tools and picard
    - DE: Differential Expression analysis with EdgeR and DESeq2
    - DEU: Differential Exon Usage analysis with EdgeR (DEXSeq skipped, runtime)
    - DAS: Differential Alternative Splicing analysis with EdgeR and DIEGO 
    - DTU: Differential Transcript Usage analysis with DRIMSeq 
    - COUNTING: Read counting with FeaturCounts
    - TRACKS: Generation of tracks for UCSC or other genome browsers
    - PEAKS: Analysis of ChIP-Seq or CLIP-Seq or cyPhyRNA-Seq Peaks

The more complex config for this analysis follows

You will see here that for DTU workflows we will need an additional FASTA and GTF file covering transcripts. As our sample data is from E.coli, transcripts and genes can be seen as synonymous, so we can easily create the needed files from our already downloaded files.

Creating a transcript gtf file from the downloaded gtf by subsetting for gene features:

.. code-block:: bash

    zcat ecoli.gtf.gz |head -n+3 > bla; zcat ecoli.gtf.gz|grep -P '\tgene\t' >> bla;cat bla|gzip > Ecoli_trans.gtf.gz;rm -f bla

We now fix the empty transcript-id tag by replacing it with the gene-id followed by "_1":abbr:

.. code-block:: bash
    
    zcat Ecoli_trans.gtf.gz |perl -wlane 'BEGIN{%ids}{if($_=~/^#/){print}else{$n=$F[9];$ids{$n}++;$F[11]=substr($n,0,-2)."_".$ids{$n}."\";";print join("\t",@F[0..7],join(" ",@F[8..$#F]))}}'|perl -F'\t' -wlane 'print $_;if($_ !~/^#/){print join("\t",@F[0..1])."\ttranscript\t".join("\t",@F[3..$#F])}'|gzip > Ecoli_trans_fix.gtf.gz

From this gtf we can now create a FASTA file by writing a helper BED file and using *BEDtools* to extract sequences from our genome FASTA file in strand specific manner. We then only need to clean up the name tag as follows:

.. code-block:: bash
    
    zcat Ecoli_trans_fix.gtf.gz|grep -P '\ttranscript\t'|perl -wlane 'next if($_=~/^#/);($name=(split(/;/,$F[11]))[0])=~s/\"//g;print join("\t",$F[0],$F[3]-1,$F[4],$name,100,$F[6])' > Ecoli_trans.bed
    bedtools getfasta -fi ecoli.fa -bed Ecoli_trans.bed -fo Ecoli_trans.fa -nameOnly -s
    perl -wlane 'if($_=~/^>/){$_=~s/\(|\+|\-|\)//g;}print' Ecoli_trans.fa |gzip > tmp.gz;cat tmp.gz ecoli.fa.gz > Ecoli_trans.fa.gz && rm -f tmp.gz
    zcat ecoli.fa.gz|grep -P '^>'|cut -d " " -f1 > salmon_decoy
    sed -i.bak -e 's/>//g' salmon_decoy

With these "Ecoli_trans.fa.gz" and "Ecoli_trans_fix.gtf.gz" files and the corresponding salmon decoy file we can now also run the DTU workflow. The config file for the exhaustive run looks as follows:

.. literalinclude:: ../../configs/tutorial_exhaustive.json
    :language: json

Note that "SETTINGS" now contains "GROUPS" which hold identifiers for the DE/DEU/DAS/DTU stages that will be used to define "COMPARABLES" that tell **MONSDA** which samples should be treated as replicates and compared to which other samples. By default this will make an all-vs-all comparison, in our case ctrl-vs-ko. Keeping this "GROUPS" setting separated from the condition-tree makes it possible to mix samples for different conditions for all pre- and processing steps and separating them for postprocessing without further struggle. This means that in this simple case we could have mixed all samples under one condition, as all other tool options stay the same and reference and annotation do not differ and would still be able to split by "GROUPS" for the DE/DEU/DAS steps. However, sometimes we need to compare data that has to be processed differently, e.g. mixing paired and single-ended reads, so we can apply different settings as needed for all workflow steps by simply splitting our condition-tree accordingly. Another addition to the "SETTINGS" key is the "IP" part. This is intended to be combined with the *PEAKS* workflow and tools that have certain preprocessing step requirements. For example our novel "scyphy" (Software for cyPhyRNA-Seq) workflow can work on 5-prime and 3-prime end enriched data, which can be indicated as "iCLIP" and "revCLIP" respectively. For this tutorial we will simulate "iCLIP" and set the corresponding key.

Somewhat special is the DTU step. Here we use the pseudo-mapper salmon to quantify expression of transcripts. These transcripts are first defined in the preprocessing part of this tutorial chapter. As we can see the REFERENCE and ANNOTATION differ from what is listed in SETTINGS, so we make use of the ability of **MONSDA** to define those settings on a per-workflow basis. This will ALWAYS over-rule the settings made in SETTINGS. So we simply define the REFERENCE and ANNOTATION key for each tool within the regular tool setting part of the config and this will make sure that these files are used instead of the globally defined ones.

Starting the run with 12 cores (defining more will be capped by the config file as MAXTHREADS is set to 12)

.. code-block:: bash

    monsda -j 12 -c ${CONDA_PREFIX}/share/MONSDA/configs/tutorial_exhaustive.json --directory ${PWD}

Will start the run in the current directory and generate a "FASTQ" sub-directory containing the downloaded sample, a "GENOME/Ecoli/INDICES" directory containing the built indices, including the one built for salmon later on, a "QC" directory containing all FASTQC reports and MULTIQC output, a "TRIMMED_FASTQ" directory for trimgalore and cutadapt output, a "DEDUP" directory for umi_tools (runs before trimming and after mapping) and picard (runs after mapping) output and a "MAPPING" directory containing the mapped files. Furthermore, "DE/DEU/DAS/DTU" directories will be created which will hold output from counting with FeatureCounts (or salmon for DTU) and DE/DEU/DAS/DTU input and output from EDGER, DESeq2, Diego and DrimSeq respectively. Again, **MONSDA** will create a "LOG" directory containing it's own log, as well as logs of all executed jobs and again a "JOBS" directory for command-line calls. 

A successful run will show the message 'Workflow finished, no error'. Be aware that this is indeed an exhaustive workflow and will require a decent amount of disk-space, memory and compute-time, depending on the hardware at your disposal.