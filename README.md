# snakes

Collection of snakemake workflows for NGS analysis from mapping to featurecount and track generation
Contains sumodule so clone with ```git clone --recursive``` or if already cloned with pull submodules with 
```
git submodule init
git submodule update 
```

See [stackoverflow](https://stackoverflow.com/questions/25200231/cloning-a-git-repo-with-all-submodules) and 
[git docs](https://git-scm.com/book/en/v2/Git-Tools-Submodules) for details.

For details on ```snakemake``` and it's features please refer to the [snakemake documentation](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html).

In general it is necessary to write a configuration file containing information on paths, files to process and settings beyond default for mapping tools and others.
For examples on this please have a look into the ```config``` directory.

For ```snakemake``` to be fully FAIR, one needs to use ```conda``` or similar environment management systems. For details on ```conda``` please refer to the [conda manual](https://docs.conda.io/en/latest/).

This workflow collection makes heavy use of ```conda``` and especially the [bioconda](https://bioconda.github.io) channel.

For distribution of jobs one can either rely on local hardware, use scheduling software like the [SGE](https://docs.oracle.com/cd/E19957-01/820-0699/chp1-1/index.html) or [Slurm](https://slurm.schedmd.com/documentation.html).

This manual will only show examples on local and SGE usage, but more information on how to use other scheduling software is available elsewhere.

## What is happening

This repository hosts the executable ```RunSnakemake.py``` which acts a wrapper around ```snakemake``` and the ```config.json``` file. 
The ```config.json``` holds all the information that is needed to run the jobs and will be parsed by ```RunSnakemake.py``` and split into sub-configs that can later be found in the directory ```SubSnakes```.

To successfully run an analysis pipeline, a few steps have to be followed:
  * Clone this repository to you working directory, or symlink it there. *DO NOT CHANGE THE NAME OF THE REPO!!!* This is needed for subworkflow generation.
  * Directory structure: The structure for the directories is dictated by the config file
  * Config file: This is the central part of the analysis. Depending on this file ```RunSnakemake.py``` will determine processing steps and generate according config and ```snakemake``` workflow files to run each subworkflow until all processing steps are done.

## Create config.json

The configuration file follows standard ```json``` format and contains all the information needed to run the analysis.
It consists of multiple sections which will be explained in detail.
For examples please refer to the ```config``` directory of the repo.

### The config.json explained

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
                                     "std": ""  #Identifier for this set of options, e.g. mapper name, standard, fancynewmethod, ...
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
The *SEQUENCING* key allows you to define *single* or *paired* as values to enable analysis of a mix of single/paired end sequences at once, defined by condition/setting.

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
                                     "std": "single" # setup and sequencing type, either paired or single
                              }
               }
    },
```

Now the actual workflow section begins, where you can define for each combinatio of processing/postprocessing step and condition/setting which environments and tool to use and which settings to apply to the run.
This follow the same scheme for each step, optionally define *RUN* ON/OFF or simply skip the key in the *WORKFLOW*/*POSTPROCESSING* section and here if not needed.
The *ENV* key defines the conda environment to load from the *env* directory of this repository, feel free to add you own environment.yaml files there.
The *BIN* key defines the name of the executable, this is needed in case the env and the bin differ as e.g. for the mapping tool ```segemehl/segemehl.x```.
The enxt key is the *OPTIONS* key which is where you can define additional parameters for each tool. It is not needed to define anything related to *single/paired* end sequencing, this is done automatically.
To add parameters simply add the *OPTION* key which holds as value a list of hashes. Parameters are defined in this hashes again as key/value pairs corresponding to the parameter name and the setting.
This should become clear having a look at the different processing steps.

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
                        "k14"
                    ]
                }
            }
        }
    },
    #Count options
    "COUNTING": { #options for trimming for each sample/condition
         "Dm6": {
            "untreated": {
                "std": {# See above
                    "ENV": "countreads", #see QC
                    "BIN": "featurecounts",
                    "OPTIONS":
                    [
                        {  # counting options here, --paired is not required, will be resolved by rules, annotation is resolved from Annotation options
                           "-f": "",
                           "--fraction": "",
                           "-p": "",
                           "-O": "",
                           "-M": "",
                           "-T": "5",
                           "-t": "exon"
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
                    "OPTIONS":
                    [
                        {
                            "-w": "ON", #-w ON enables one line per feature annotation, including start/end of the feature, output can become quite large, disable by adding empty string ("") as value
                            "ANNOFEATURE" : "", #You can specify a set of certain features to annotate here, e.g. 'exon' will only annotate exon overlaps, disable specific feature annotation by adding empty string ("") as value
                            "ANNOTATIONFILE": "dm6.gff.gz"
                        }
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
                    "OPTIONS":
                    [
                        {
                          "-n": "DM6 Standard Mapping", #name of the hub
                          "-s" : "dm6_st", #short name for hub
                          "-l" : "UCSC DM6 Standard Mapping", #long name for track
                          "-b" : "UCSC dm6 std", #short name for track
                          "ANNOTATION": "dm6.gff.gz"
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
                             "CONDITION":  ["WT","WT"], # define the conditions per sample as list of same lenght as replicates
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
                             "CONDITION":  ["KO","KO"],
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


### Cluster config
This is separate from the main configuration, for details please follow the explanations in the [snakemake documentation](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html).
For an example have a look at the ```cluster``` directory.

## Run the pipeline
Simply run

```
python snakes/RunSnakemake.py
```

to see the help and available options that will be passed through to ```snakemake```.

To start a simple run call
```
python snakes/RunSnakemake.py -j NUMBER_OF_CORES --configfile YOUR_CONFIG.json --directory ${PWD}
```
or add additional arguments for ```snakemake``` as you see fit.

### Run on cluster

Define the cluster config file and for SGE support simply append ```--cluster "qsub -q ${QUEUENAME} -e ${PWD}/sgeerror -o ${PWD}/sgeout -N ${JOBNAME}" --jobscript snakes/cluster/sge.sh```

## TODO
Implementation of Slurm support, Isoform analysis, you name it
