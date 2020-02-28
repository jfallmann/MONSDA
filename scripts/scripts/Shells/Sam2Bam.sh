#!/bin/bash

file=$1 ### Name of the sam file you want to convert
ref=$2  ### Path to reference genome fasta file
bins=$3
out=$4
threads=$5

echo "running samtools view -bT $ref -o $out --threads $threads $file"

if [ ! -f $out ];then
        echo "$out not found, creating new"
        zcat $file|samtools view -bT $ref -o $out --threads $threads -
fi
if [ ! -f $out".bai" ];then
       echo "$out.bai not found, creating new"
        samtools index $out
fi

