#!/bin/bash

in=$1
out=$2
threads=$3

if [[ "$1" == *.gz* ]]
then
    samtools view -H <(zcat $in) | grep '@HD' | pigz -p $threads -f > $out
    samtools view -H <(zcat $in) | grep '@SQ' | sort -t$'\t' -k1,1 -k2,2V | pigz -p $threads -f >> $out
    samtools view -H <(zcat $in) | grep '@RG' | pigz -p $threads -f >> $out
    samtools view -H <(zcat $in) | grep '@PG' | pigz -p $threads -f >> $out
else
    samtools view -H <(cat $in)|grep '@HD' | pigz -p $threads -f > $out
    samtools view -H <(cat $in)|grep '@SQ' | sort -t$'\t' -k1,1 -k2,2V |pigz -p $threads -f >> $out
    samtools view -H <(cat $in)|grep '@RG' | pigz -p $threads -f >> $out
    samtools view -H <(cat $in)|grep '@PG' | pigz -p $threads -f >> $out
fi

if [[ "$1" == *-bwa* ]]
then
    if [[ "$1" == *.gz* ]]
    then
        zcat $in | grep -v "^@"| grep -v -e $'\t''XA:Z:' -e $'\t''SA:Z:' | pigz -p $threads -f >> $out
    else
        cat $in | grep -v "^@"| grep -v -e $'\t''XA:Z:' -e $'\t''SA:Z:' | pigz -p $threads -f >> $out
    fi
else
    if [[ "$1" == *.gz* ]]
    then
        zcat $in | grep -v "^@" | grep -w -P "NH:i:1|tp:A:P" | pigz -p $threads -f >> $out
    else
        cat $in | grep -v "^@" | grep -w -P "NH:i:1|tp:A:P" | pigz -p $threads -f >> $out
    fi
fi
