=========================
The config.json explained
=========================

It consists of multiple sections which will be explained in detail.  For the template please refer to the
``template.json`` file in the ``config`` directory of the repo.  The ``template.json`` is the default
blueprint for a config file used by ``Configurator.py``.

The config file contains all the information needed to run stated workflows and to find the sample/genome
files.  It starts with a key/value pair defining which workflows and pre-/postprocessing steps to run. Be
aware that every workflow and postproccessing value has to correspond to a key later in the config.json that
defines parameters specific for the job:

::

   "WORKFLOWS": "SRA,QC,MAPPING,TRIMMING,COUNTING,UCSC,DE,DEU,DAS", # Here you define which main workflow steps should be run,
   "REFERENCE": "GENOMES", #where to find the reference genomes
   "BINS": "nextsnakes/scripts", #where to find the scripts used in the workflow, if you soft-link or clone the snake git to your working directory use this path
   "MAXTHREADS": "20", #maximum number of cores to use, make sure your cluster/machine can handle the load

The next part defines where the path to the genome files and its main name plus an extension in case you use
specific genomes for different runs.  The directory structure that ``RunSnakemake.py|RunNextflow.py`` will
follow is in this example *GENOME\Dm6* and it would look for a genome file named *dm6.fa.gz* without further
extension.  Here already the split by condition is applied, so you can define different genomes to use for
different conditions/settings, e.g. when running on standard RNA-Seq and Bisulfite-Seq in parallel or on
different genomes at once.  In case you use different genomes just add them as key/value pairs to the *GENOME*
key of the config file, the ``.fa.gz`` postfix is *ALWAYS* assumed.  In the *SOURCE* section you then define
which condition/setting should use which genome by stating the *GENOME* key as innermost value of the *SOURCE*
key.

::

    "GENOME": { #which genomes to use and how the reference fasta is named, key is subdir of REFERENCE and value is name of fasta
                "Dm6": "dm6"
              },
    "NAME": { #extension for genome file name in case different types of reference are used, e.g. genomic, artificial, organelle ...
              "Dm6": {#first key should be same than for GENOME
                      "untreated": { #Whatever identifier you want to use for this dataset e.g. Date, Labname, treatment, ...
                                     "std": ""  #Key is identifier for this set of options, e.g. mapper name, standard, fancynewmethod, ... Value is either empty or can be a string that is part of the genome name, e.g. extended, artificial ...
                                   }
                     }
            },
    "SOURCE": {  #which organisms are to be parsed for which samples
                 "Dm6": { #key is subdir of REFERENCE, see GENOME, you can specify which genomes to use for which dataset identifier, e.g. untreated with setting std will use dm6 here
                          "untreated": {
                              "std": "Dm6"
                          }
                        }
              }


The next part defines the samples to run the analysis on, just add a list of sample names as innermost value
to the *SAMPLES* key for each condition.  In case of single-end sequencing make sure to include the _R1 _R2
tag, in case of paired end skip those as the pipeline will automatically look for _R1 and _R2 tags to find
read pairs.  *Make sure the naming of you samples follows this _R1 _R2 convention when running paired-end
analysis!* The *SEQUENCING* key allows you to define *unpaired* or *paired* as values to enable analysis of a
mix of single/paired end sequences at once, defined by condition/setting.  You can also specify strandedness
of the protocol used, if unstranded leave empty, else add strandedness according to
http://rseqc.sourceforge.net/#infer-experiment-py as comma separated value (rf Assumes a stranded library
fr-firststrand [1+-,1-+,2++,2--], fr Assumes a stranded library fr-secondstrand [1++,1--,2+-,2-+])

::

    "SAMPLES": {  #which samples to analyze
                  "Dm6": { #key for source and genome
                           "untreated": {      # sample id
                                               "std": ["GSM461177_untreat_paired_subset_r1","GSM461177_untreat_paired_subset_r2"] # setup and list of samples you whish to analyze
                                        }
                         }
               },
    "SEQUENCING" : {
        "Dm6": { #key for source and genome
                 "untreated": {      # sample id
                                     "std": "unpaired" # setup and sequencing type, either paired or unpaires, stranded or unstranded, if unstranded leave empty, if stranded see below
                                     #"std": "paired,fr" # if stranded add strandedness according to http://rseqc.sourceforge.net/#infer-experiment-py as comma separated value (rf Assumes a stranded library fr-firststrand [1+-,1-+,2++,2--], fr Assumes a stranded library fr-secondstrand [1++,1--,2+-,2-+])
                              }
               }
    }


Now the actual workflow section begins, where you can define for each combinatio of processing/postprocessing
step and condition/setting which environments and tool to use and which settings to apply to the run.  These
follow the same scheme for each step, optionally define *RUN* ON/OFF or simply skip the key in the
*WORKFLOW*/*POSTPROCESSING* section and here if not needed.  The *ENV* key defines the conda environment to
load from the *env* directory of this repository, feel free to add you own environment.yaml files there.  The
*BIN* key defines the name of the executable, this is needed in case the env and the bin differ as e.g. for
the mapping tool ``segemehl/segemehl.x``.  The next key is the *OPTIONS* key which is where you can define
additional parameters for each tool. It is not needed to define anything related to *unpaired/paired* end
sequencing, this is done automatically.  To add parameters simply add the *OPTION* key which holds as value a
list of hashes. Parameters are defined in this hashes again as key/value pairs corresponding to the parameter
name and the setting.  This should become clear having a look at the different processing steps.  If there are
no options just do not add the *OPTION* key

::

    "QC": {
        "RUN": "ON", #set to 'OFF' to skip QC
        "Dm6": { #key for source and genome
                 "untreated": {      # sample id
                                     "std": {
                                         "ENV" : "fastqc",  # name of conda env for QC
                                         "BIN" : "fastqc" # binary for trimming
                                     }
                              }
               }
    },
    "TRIMMING": { #options for trimming for each sample/condition
                  "RUN": ON", # set to 'OFF' if no trimming wanted
        "Dm6": {
            "untreated": {
                "std": { # See above
                    "ENV": "trimgalore", # name of conda env for trimming
                    "BIN": "trim_galore", # name of binary for trimming
                    "OPTIONS":
                    [
                        {  # trimming options here, --paired is not required, will be resolved by rules
                            "-q": "15",
                            "--length": "8", #READ_MINLEN discard reads shorter than that
                            "-e": "0.15"
                        }
                    ]
                }
            }
        }
    },
    "MAPPING": { #options for mapping for each sample/condition
        "Dm6": {
            "untreated": {
                "std": {# first entry in list is a dict of options for indexing, second for mapping, third can be e.g. appendix to index name, useful especially with minimap if using different kmer sizes
                    "ENV": "minimap", # which conda env to use for mapping
                    "BIN": "minimap2", #how the mapper binary is called
                    "OPTIONS":
                    [
                        {
                            "-k": "14"#option for setting kmer size while indexing
                        },
                        {
                            "-ax": "map-ont",
                            "-ub": "",
                            "-Y": "",
                            "-L": "",
                            "--MD": "",
                            "-d": ""
                        },
                        "k14" #name the index that is generated, if this is left empty the index will have the extention 'std'
                    ]
                }
            }
        }
    },
    "COUNTING": { #options for read counting for each sample/condition
        "FEATURES": { #which features to count (KEY) and which group they belong to (VALUE)
            "exon": "Parent",
            "gene": "ID"
        },
         "Dm6": {
            "untreated": {
                "std": {# See above
                    "ENV": "countreads", #see QC
                    "BIN": "featurecounts",
                    "OPTIONS":
                    [
                        {  # counting options here, --paired is not required, will be resolved by rules, annotation is resolved from ANNOTATION option, feature and group is resolved by the FEATURES key
                           "-f": "",
                           "--fraction": "",
                           "-p": "",
                           "-O": "",
                           "-M": "",
                           "-T": "5"
                        }
                    ]
               }
           }
       }
    },
    "ANNOTATE" : {
         "Dm6": {
            "untreated": {
                "std": { # See above
                    "ENV" : "annotatebed",
                    "BIN" : "annotate", #dummy as ucsc has no direct bin but we need the key
                    "ANNOFEATURE" : "", #You can specify a set of certain features to annotate here, e.g. 'exon' will only annotate exon overlaps, disable specific feature annotation by adding empty string ("") as value
                    "ANNOTATIONFILE": "dm6.gff.gz",
                    "OPTIONS":
                    [
                        {
                            "-w": "ON" #-w ON enables one line per feature annotation, including start/end of the feature, output can become quite large, disable by adding empty string ("") as value                        }
                    ]
                }
            }
         }
     },
    "UCSC" : {
         "Dm6": {
            "untreated": {
                "std": { # See above
                    "ENV" : "ucsc",
                    "BIN" : "ucsc", #dummy as ucsc has no direct bin but we need the key
                    "ANNOTATION": "dm6.gff.gz",
                    "OPTIONS":
                    [
                        {
                          "-n": "DM6 Standard Mapping", #name of the hub
                          "-s" : "dm6_st", #short name for hub
                          "-l" : "UCSC DM6 Standard Mapping", #long name for track
                          "-b" : "UCSC dm6 std" #short name for track
                       }
                    ]
                }
            }
         }
     }


Nextsnakes further supports DE/DEU/DAS-Analysis as postprocessing steps for a defined set of samples. The config for this step looks as follows:

::

    #DE/DEU/DAS options
	"DAS" : { # this can be DE, DEU or DAS
	    "TOOLS" : #in contrast to other analysis types you can already define a set of tools at this stage that will be run sequentially
        {
            "edger"  : "Analysis/DAS/EDGER.R",
            "diego"  : "diego.py"
        },
        "COMPARABLE" : #Here you can set the actual comparisons you are interested in, leace empty for ALLvsALL pairwise comparisons
        {
            "contrast_WTvsKOs": [["WT"],["KO1","KO2"]]
        },
        "id": {
            "condition": {
                "setting": {
                    "ANNOTATION": "genome_or_other.gtf.gz", #gtf file for featurecount and dexseq/edger
                    "GROUPS":  ["WT","KO1","KO2"], #Conditions of samples can be different than the condition setting
                    "REPLICATES": ["SAMPLE1_r1","SAMPLE2_r2","SAMPLE2_r3"], #replicates that belong to condition, one entry here for one entry in GROUPS
                    "TYPES": ["standard","standard","standard"], #sequencing type or additional condition to compare to, can be empty
                    "OPTIONS":
                    [
                        {# this options are used for the featurecount rule, there is no need to run COUNTING prior to DE/DEU/DAS as specific processing of count tables is needed anyway
                            "-t": "exon",
                            "-g": "gene_id",
                            "-f": "",
                            "--fraction": "",
                            "-O": ""
                        }
                    ]
                }
            }
        }
    }


Keep in mind that every workflow/postprocessing step needs a corresponding entry in the config file or
``RunSnakemake.py|RunNextflow.py`` will throw an error.
