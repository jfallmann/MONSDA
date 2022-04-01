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
   * - tools
     - blabla
     - env/bin
     - `LINK <https://github.com/>`_


Mapping
=======

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

Deduplication
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
================================

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
=====================

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

Track Generator
=======================

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

Peak-calling
============

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
