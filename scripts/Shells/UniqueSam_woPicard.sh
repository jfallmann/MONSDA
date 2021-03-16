#!/bin/bash

#/home/mescalin/fall/bin/java/bin/java -jar /home/mescalin/fall/bin/picard-tools-1.66/ViewSam.jar INPUT=$1 |grep "^@" > $1\_unique.sam
#samtools view -HS $1 |pigz > $3
#samtools view -H $1|grep '@HD' |pigz -f > $3 && samtools view -H $1|grep '@SQ'|sort -t$'\t' -k1,1 -k2,2V |pigz -f >> {output[1]} && samtools view -H {input[0]}|grep '@PG'|pigz -p {threads} -f >> {output[1]}
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
        zcat $in | grep -v "^@"| grep -v -e $'\t''XA:Z:' -e $'\t''SA:Z:' | pigz >> $out
    else
        cat $in | grep -v "^@"| grep -v -e $'\t''XA:Z:' -e $'\t''SA:Z:' | pigz >> $out
    fi
else
    if [[ "$1" == *.gz* ]]
    then
        zcat $in | grep -v "^@" | grep -w -P "\tNH:i:1|\ttp:A:P" | pigz >> $out
    else
        cat $in | grep -v "^@" | grep -w -P "\tNH:i:1|\ttp:A:P" | pigz >> $out
    fi
fi
