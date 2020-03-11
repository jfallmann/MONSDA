THE config.json EXPLAINED
===============================

It consists of multiple sections which will be explained in detail.
For examples please refer to the ```config``` directory of the repo.
The ```skeleton.json``` is the default blueprint for a config file used by ```Configurator.py```.

The config file contains all the information needed to run stated workflows and to find the sample/genome files.
It starts with a key/value pair defining which workflows and postprocessing steps to run. Be aware that every worklow and postproccessing value has to correspond to a key later in the config.json that defines parameters specific for the job:

```
{
    "WORKFLOWS": "MAPPING,TRIMMING,QC", # Here you define which main workflow steps should be run,
    "POSTPROCESSING" : "COUNTING,UCSC,ANNOTATE", # no specific order needed
    "REFERENCE": "GENOMES", #where to find the reference genomes
    "BINS": "snakes/scripts", #where to find the scripts used in the workflow, if you soft-link the snake git to Workflow, use this path
    "MAXTHREADS": "20", #maximum number of cores to use
```

The next part defines where the path to the genome files and its main name plus an extension in case you use specific genomes for different runs.
The directory structure that ```RunSnakemake.py``` will follow is in this example *GENOME\Dm6* and it would look for a genome file named *dm6* without extension.
Here already the split by condition is applied, so you can define different genomes to use for different conditions/settings, e.g. when running on standard RNA-Seq and Bisulfite-Seq in parallel or on different genomes at once.
In case you use different genomes just add them as key/value pairs to the *GENOME* key of the config file.
In the *SOURCE* section you then define which condition/setting should use which genome by stating the *GENOME* key as innermost value of the *SOURCE* key.

```
    "GENOME": { #which genomes to use and how the reference fasta is named, key is subdir of REFERENCE and value is name of fasta
                "Dm6": "dm6"
              },
    "NAME": { #extension for genome file name in case different typed of reference are used, e.g. genomic, artificial, organelle ...
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
              },
```

The next part defines the samples to run the analysis on, just add a list of sample names as innermost value to the *SAMPLES* key for each condition.
In case of single-end sequencing make sure to include the _r1 _r2 tag, in case of paired end skip those as the pipeline will look for _r1 and _r2 tags to find read pairs.
*Make sure the naming of you samples follows this _r1 _r2 convention when running paired-end analysis!*
The *SEQUENCING* key allows you to define *unpaired* or *paired* as values to enable analysis of a mix of single/paired end sequences at once, defined by condition/setting.
You can also specify strandedness of the protocol used, if unstranded leave empty, else add strandedness according to http://rseqc.sourceforge.net/#infer-experiment-py as comma separated value (rf Assumes a stranded library fr-firststrand [1+-,1-+,2++,2--], fr Assumes a stranded library fr-secondstrand [1++,1--,2+-,2-+])

```
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
    },
```

Now the actual workflow section begins, where you can define for each combinatio of processing/postprocessing step and condition/setting which environments and tool to use and which settings to apply to the run.
This follow the same scheme for each step, optionally define *RUN* ON/OFF or simply skip the key in the *WORKFLOW*/*POSTPROCESSING* section and here if not needed.
The *ENV* key defines the conda environment to load from the *env* directory of this repository, feel free to add you own environment.yaml files there.
The *BIN* key defines the name of the executable, this is needed in case the env and the bin differ as e.g. for the mapping tool ```segemehl/segemehl.x```.
The next key is the *OPTIONS* key which is where you can define additional parameters for each tool. It is not needed to define anything related to *unpaired/paired* end sequencing, this is done automatically.
To add parameters simply add the *OPTION* key which holds as value a list of hashes. Parameters are defined in this hashes again as key/value pairs corresponding to the parameter name and the setting.
This should become clear having a look at the different processing steps.
If there are no options just do not add the *OPTION*

```
#QC options
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
#Trimming options
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
    #Mapping software options
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
    #Count options
    "COUNTING": { #options for trimming for each sample/condition
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
    #Annotation options
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
                          "-b" : "UCSC dm6 std", #short name for track
                       }
                    ]
                }
            }
         }
     }
}
```

The pipeline now also supports DE-Analysis as postprocessing step for a defined set of samples. The config for this step looks as follows:

```
    #DE options
    "DE": { #options for each sample/condition
            "Dm6": {
                "unpaired": {
                    "star": {
                             "ENV": "deseq2", #choose which environment to use (currently only deseq2 is supported, others will follow)
                             "BIN": "deseq2", #choose which binary to use (currently only deseq2 is supported, others will follow)
                             "GROUP":  ["WT","WT"], # define the conditions per sample as list of same lenght as replicates
                             "REPLICATES": ["GSM461177_untreat_paired_subset_r1","GSM461177_untreat_paired_subset_r2"], #this are the replicates for each condition, length has to match with CONDITION
                             "OPTIONS": #This is not needed currently, other DE pipelines may need this later on
                             [
                                 {  #Should not be needed
                                 }
                             ]
                            }
                },
                "paired": {#see above, in principle conditions can also be mixed between subsets, this example tests all unpaired reads as WT versus all paired reads as KO
                    "star": {
                             "ENV": "deseq2", #see QC
                             "BIN": "deseq2",
                             "GROUP":  ["KO","KO"],
                             "REPLICATES": ["GSM461177_untreat_paired_subset","GSM461177_untreat_paired_subset"],
                             "OPTIONS":
                             [
                                 {  #Should not be needed
                                 }
                             ]
                            }
                }
            }
          }

```

Keep in mind that every workflow/postprocessing step needs a corresponding entry in the config file or ```RunSnakemake.py``` will throw an error.
