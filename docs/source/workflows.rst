PREPROCESSING
=============

FETCH
#####

.. list-table::
   :widths: 25 50 25 50
   :header-rows: 1

   * - Tool
     - Description
     - Env/Bin
     - Link
   * - SRAtools fasterq-dump
     - The fasterq-dump tool extracts data in FASTQ- or FASTA-format from SRA-accessions
     - sra/fasterq-dump
     - `LINK <https://github.com/ncbi/sra-tools>`_


BASECALL
########

.. list-table::
   :widths: 25 50 25 50
   :header-rows: 1

   * - Tool
     - Description
     - Env/Bin
     - Link
   * - Guppy
     - Data processing toolkit that contains Oxford Nanopore’s basecalling algorithms, and several bioinformatic post-processing features, such as barcoding/demultiplexing, adapter trimming, and alignment. Needs to be installed locally as no ``conda`` version is available
     - guppy/$PATH_to_local_installation
     - `LINK <https://nanoporetech.com/nanopore-sequencing-data-analysis>`_


Quality Control
################

This workflow step can be run as preprocessing step if none of the processing workflows is defined in the config.json.

.. list-table::
   :widths: 25 50 25 50
   :header-rows: 1

   * - Tool
     - Description
     - Env/Bin
     - Link
   * - FASTQC
     - A quality control tool for high throughput sequence data.
     - fastqc/fastqc
     - `LINK <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_


PROCESSING
==========

Quality Control
###############

If any of the below listed processing steps is defined in the config.json, quality control will run for all generated output files if enabled. 

.. list-table::
   :widths: 25 50 25 50
   :header-rows: 1

   * - Tool
     - Description
     - Env/Bin
     - Link
   * - FASTQC
     - A quality control tool for high throughput sequence data.
     - fastqc/fastqc
     - `LINK <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`_


Trimming
########

.. list-table::
   :widths: 25 50 25 50
   :header-rows: 1

   * - Tool
     - Description
     - Env/Bin
     - Link
   * - trim_galore
     - A wrapper tool around Cutadapt and FastQC to consistently apply quality and adapter trimming to FastQ files, with some extra functionality for MspI-digested RRBS-type (Reduced Representation Bisufite-Seq) libraries.
     - trimgalore/trim_galore
     - `LINK <https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/>`_
   * - Cutadapt
     - Cutadapt finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads.
     - cutadapt/cutadapt
     - `LINK <https://cutadapt.readthedocs.io/en/stable/>`_


Mapping
#######

.. list-table::
   :widths: 25 50 25 50
   :header-rows: 1

   * - Tool
     - Description
     - Env/Bin
     - Link
   * - HISAT2
     - HISAT2 is a fast and sensitive alignment program
     - hisat2/hisat2
     - `LINK <http://daehwankimlab.github.io/hisat2/manual/>`_
   * - STAR
     - Spliced Transcripts Alignment to a Reference
     - star/star
     - `LINK <https://github.com/alexdobin/STAR>`_
   * - Segemehl2|3
     - Segemehl is a software to map short sequencer reads to reference genomes.
     - segmehl2|3/segemehl.x
     - `LINK <https://www.bioinf.uni-leipzig.de/Software/segemehl/>`_
   * - BWA
     - BWA is a software package for mapping low-divergent sequences against a large reference genome
     - bwa/bwa mem
     - `LINK <http://bio-bwa.sourceforge.net/>`_
   * - Minimap2
     - Minimap2 is a versatile sequence alignment program that aligns DNA or mRNA sequences against a large reference database. 
     - minimap/minimap2
     - `LINK <https://github.com/lh3/minimap2>`_    


DEDUP
=============

.. list-table::
   :widths: 25 50 25 50
   :header-rows: 1

   * - Tool
     - Description
     - Env/Bin
     - Link
   * - tools
     - blabla
     - env/bin
     - `LINK <https://github.com/>`_


POSTPROCESSING
==============

Read-Counting and Quantification
################################

.. list-table::
   :widths: 25 50 25 50
   :header-rows: 1

   * - Tool
     - Description
     - Env/Bin
     - Link
   * - tools
     - blabla
     - env/bin
     - `LINK <https://github.com/>`_

Differential Analyses
#####################

+-----------+-------------------------------------+------------------+-----------------+----------------+---------------------------------+----------------+------------------------------------------------------+-----------------------------------------+-----------------------------------------+-------------------+-------------------------------------------------------------------+-------+
| Tool      | Analysis                            | Filtering        | Normalization   | Distribution   | Testing                         | Significance   | Results Table                                        | further                                 | SigTables                               | Clustering        | further                                                           | Rmd   |
+===========+=====================================+==================+=================+================+=================================+================+======================================================+=========================================+=========================================+===================+===================================================================+=======+
| edgeR     | Differential Gene Expression        | filterByExpr()   | TMM             | NB             | Fisher’s exact test             | pValue, LFC    | results, sorted-results                              | normalized                              | Sig, SigUP, SigDOWN                     | MDS-plot          | BCV, QLDisp, MD(per comparison)                                   | ✓     |
+-----------+-------------------------------------+------------------+-----------------+----------------+---------------------------------+----------------+------------------------------------------------------+-----------------------------------------+-----------------------------------------+-------------------+-------------------------------------------------------------------+-------+
| edgeR     | Differential Exon Usage             | filterByExpr()   | TMM             | NB             | Fisher’s exact test             | pValue, LFC    | results                                              | normalized                              |                                         | MDS-plot          | BCV, QLDisp, MD(per comparison)                                   | ✓     |
+-----------+-------------------------------------+------------------+-----------------+----------------+---------------------------------+----------------+------------------------------------------------------+-----------------------------------------+-----------------------------------------+-------------------+-------------------------------------------------------------------+-------+
| edgeR     | Differential Alternative Splicing   | filterByExpr()   | TMM             | NB             | Simes, gene-level, exon-level   | pValue, LFC    | results(diffSpliceExonTest, Simes-Test, Gene-Test)   |                                         | Sig, SigUP, SigDOWN                     | MDS-plot          | BCV, QLDisp, MD(per comparison), topSpliceSimes-plots(per Gene)   | ✓     |
+-----------+-------------------------------------+------------------+-----------------+----------------+---------------------------------+----------------+------------------------------------------------------+-----------------------------------------+-----------------------------------------+-------------------+-------------------------------------------------------------------+-------+
| DESeq2    | Differential Gene Expression        | RowSums >= 10    | RLE             | NB             | Wald test                       | pValue, LFC    | results                                              | rld, vsd, results(per comparison)       | Sig, SigUP, SigDOWN                     | PCA               | Heatmaps, MA(per comparison), VST-and-log2                        | ✓     |
+-----------+-------------------------------------+------------------+-----------------+----------------+---------------------------------+----------------+------------------------------------------------------+-----------------------------------------+-----------------------------------------+-------------------+-------------------------------------------------------------------+-------+
| DEXSeq    | Differential Exon Usage             | RowSums >= 10    | RLE             | Cox-Reid       | likelihood ratio test           |                |                                                      |                                         |                                         |                   |                                                                   |       |
+-----------+-------------------------------------+------------------+-----------------+----------------+---------------------------------+----------------+------------------------------------------------------+-----------------------------------------+-----------------------------------------+-------------------+-------------------------------------------------------------------+-------+
| DEXSeq    | Differential Transcript Usage       | dmFilter()       | RLE             | Cox-Reid       | likelihood ratio test           | pValue         | results                                              |                                         |                                         |                   |                                                                   | ✓     |
+-----------+-------------------------------------+------------------+-----------------+----------------+---------------------------------+----------------+------------------------------------------------------+-----------------------------------------+-----------------------------------------+-------------------+-------------------------------------------------------------------+-------+
| DIEGO     | Differential Alternative Splicing   |                  |                 |                | Mann-Whitney U test             | pValue         | results                                              |                                         | Sig                                     | Dendrogram-plot   |                                                                   | ✓     |
+-----------+-------------------------------------+------------------+-----------------+----------------+---------------------------------+----------------+------------------------------------------------------+-----------------------------------------+-----------------------------------------+-------------------+-------------------------------------------------------------------+-------+
| DRIMSeq   | Differential Transcript Usage       | dmFilter()       |                 | DM             |                                 | pValue, LFC    | results(transcript, genes)                           | Proportions-table, genewise precision   | Sig, SigUP, SigDOWN (transcipt, gene)   |                   | FeatPerGene, precision, Pvalues (per comparison)                  | ✓     |
+-----------+-------------------------------------+------------------+-----------------+----------------+---------------------------------+----------------+------------------------------------------------------+-----------------------------------------+-----------------------------------------+-------------------+-------------------------------------------------------------------+-------+

TRACKS
###############

.. list-table::
   :widths: 25 50 25 50
   :header-rows: 1

   * - Tool
     - Description
     - Env/Bin
     - Link
   * - tools
     - blabla
     - env/bin
     - `LINK <https://github.com/>`_

PEAKS
#####

.. list-table::
   :widths: 25 50 25 50
   :header-rows: 1

   * - Tool
     - Description
     - Env/Bin
     - Link
   * - tools
     - blabla
     - env/bin
     - `LINK <https://github.com/>`_
