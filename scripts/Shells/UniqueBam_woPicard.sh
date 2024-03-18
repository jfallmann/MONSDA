#!/bin/bash

in=$1
out=$2
threads=$3
bwa="${4:-}"

samtools view -H ${in} | grep '@HD' > head
samtools view -H ${in} | grep '@SQ' | sort -t$'\t' -k1,1 -k2,2V >> head
samtools view -H ${in} | grep '@RG' >> head
samtools view -H ${in} | grep '@PG' >> head


if [[ "${in}" == *bwa* ]] || [[ -n "${bwa}" ]]
then
    cat head <(samtools view --threads ${threads} -F 4 ${in} | grep -v "^@"| grep -v -e $'\t''XA:Z:' -e $'\t''SA:Z:') | samtools view --threads ${threads} -hb - > ${out}
else
    cat head <(samtools view --threads ${threads} -F 4 ${in} | grep -v "^@" | grep -w -P "NH:i:1|tp:A:P") | samtools view --threads ${threads} -hb - > ${out}
fi

rm -f head
