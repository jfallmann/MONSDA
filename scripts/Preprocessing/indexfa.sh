#!/usr/bin/env bash

inf=$1

if [[ -s $inf ]]
then
	if [[ "$inf" == *.gz* ]]	   
	then
		filename="${inf%.*}"
		zcat $inf > $filename && samtools faidx $filename && rm -f $filename
	else
		samtools faidx $inf
	fi
else
	touch $inf.fai
fi
