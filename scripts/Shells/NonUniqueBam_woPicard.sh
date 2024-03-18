#!/bin/bash

in=$1
out=$2
threads=$3
bwa="${4:-}"

samtools view -H ${in} | grep '@HD' > nhead
samtools view -H ${in} | grep '@SQ' | sort -t$'\t' -k1,1 -k2,2V >> nhead
samtools view -H ${in} | grep '@RG' >> nhead
samtools view -H ${in} | grep '@PG' >> nhead


if [[ "${in}" == *bwa* ]] || [[ -n "${bwa}" ]]
then
    cat nhead <(samtools view --threads ${threads} -F 4 ${in} | grep -v "^@"| grep -e $'\t''XA:Z:' -e $'\t''SA:Z:') | samtools view --threads ${threads} -hb - > ${out}
else
    cat nhead <(samtools view --threads ${threads} -F 4 ${in} | grep -v "^@" | grep -v -w -P "NH:i:0|NH:i:1|tp:A:P") | samtools view --threads ${threads} -hb - > ${out}
fi

rm -f nhead
