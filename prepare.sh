#!/usr/bin/env bash

echo "Creating skeleton dirs"

mkdir -p FASTQ TRIMMED_FASTQ RAW QC LOGS

#now link scripts and bins to the snakemake directory, e.g.
#ln -s snakes/scripts .
