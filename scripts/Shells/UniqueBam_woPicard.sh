#!/bin/bash

in=$1
out=$2
threads=$3
specialmappers="${4:-}"

samtools view -H ${in} | grep '@HD' > ${in}_head
samtools view -H ${in} | grep '@SQ' | sort -t$'\t' -k1,1 -k2,2V >> ${in}_head
samtools view -H ${in} | grep '@RG' >> ${in}_head
samtools view -H ${in} | grep '@PG' >> ${in}_head


if [[ "${in}" == *bwa* ]] || [[ "$specialmappers" == *bwa* ]]
then
    cat ${in}_head <(samtools view --threads ${threads} -F 4 ${in} | grep -v "^@"| grep -v -e $'\t''XA:Z:' -e $'\t''SA:Z:') | samtools view --threads ${threads} -hb - > ${out}
elif [[ "$1" == *minimap* ]] || [[ "$specialmappers" == *minimap* ]]
then
    cat ${in}_head <(samtools view --threads ${threads} -F 4 ${in} | grep -v "^@" | perl -wlane 'print if $F[4] >=60' | samtools view --threads ${threads} -hb - > ${out}
else
    cat ${in}_head <(samtools view --threads ${threads} -F 4 ${in} | grep -v "^@" | grep -w -P "NH:i:1|tp:A:P") | samtools view --threads ${threads} -hb - > ${out}
fi

rm -f ${in}_head
