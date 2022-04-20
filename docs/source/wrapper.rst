Wrapping Workflows
==================

In general **MONSDA** is *Python3* software, that wraps workflows by assembling subworkflows or single tools from `.smk` or `.nf` templates. The idea here is that tools and subworkflows for similar tasks are designed in a way that starts from the same input and results in the same output. This is not only true for single workflow steps which can be performed by multiple tools, but also for the wrapped workflow management systems (WMS). In principal output generated in `Nextflow` mode should be suitable as input for `Snakemake` and vice versa. This means that, for example, mapping output generated in `Nextflow` mode can be used as input for *DE* analysis in `Snakemake` mode, while both work from the same **config.json**.

As Snakemake is also written in *Python*, wrapping workflows is similar to the built-in way of submodule assembly, although we take care that submodules for the same task remain interchangeable. Wrapping Nextflow is slightly different, as `MONSDA` has to assemble *Groovy* text blocks, which does not make any difference to the end user, but requires to translate configuration from the config file to Nextflow parsable command lines. However, the idea of creating interchangeable subworkflows or tool specific code fragments stays the same.

Independently of the wrapped WMS, workflows are split internally in three independent stages. *PREPROCESSING* includes all workflow steps that generate or manipulate FASTQ files, to make them available to the *PROCESSING* stage. This includes download of files from SRA, basecalling with Guppy and pre-quality-control, so all workflow steps that do not require identical input formats, but lead to similar output.

*PROCESSING* starts from FASTQ files and includes trimming, deduplication, mapping and quality control for all subprocesses.

*POSTPROCESSING* builds upon *PROCESSING* output and includes quantification, differential expression analysis on gene, transcript and exon level, generation of tracks for UCSC or other genome browsers, peak finding and circular RNA identification. In contrast to *PROCESSING* workflows, these steps do not require output to be of similar format but are able to work from the same input.

In case dedicated workflows need to be established, as is the case for example for cyPhyRNA-Seq, the main idea is to split all preprocessing and processing steps in units that remain interchangeable, and deliver dedicated post-processing subworkflows which can work on their output. In the mentioned example we have quality control, trimming, mapping and deduplication embedded in standard workflows and dedicated, specific postprocessing of mapped reads wrapped in the PEAKS postprocessing step.

For new workflows, we aim to split those into as small subunits as possible to make subworkflows available for other pipelines and add dedicated parts as postprocessing to currently established categories. In case new categories need to be defined, please contact us to discuss how this can be embedded.

