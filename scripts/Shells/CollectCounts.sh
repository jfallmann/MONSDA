#!/bin/env bash
PATTERN=$1

for i in *$PATTERN;do
	echo $i
    grep -v '^#' $i|grep -v "^Gene"|cut -d$'\t' -f1,7 > ToJoin_$i
done

cut -d$'\t' -f1 ToJoin_L12270_09_WT_d4_1_2_IN_R1_mapped_sorted.counts > tmp_counts

for i in ToJoin_*$PATTERN;do
	echo $i
	join -1 1 -2 1 <(sort -k1,1d tmp_counts) <(sort -k1,1d $i) > COUNTS_$PATTERN && cp -f COUNTS_$PATTERN tmp_counts
done
rm -f tmp_counts

awk 'BEGIN{FS=" ";OFS=FS}{t=0;for(i=2;i<=NF;i++){t+=$i} print $1,t;t=0}' COUNTS_$PATTERN |sort -k1,1d > Gene_sum_$PATTERN

#join -1 1 -2 1 <(sort -k1,1d Gene_sum_bfx656.hg38.e81) <(sort -k1,1d COUNTS/Featurecounter/bfx656/Gene_sum_mapped_sorted_unique.counts) > Compare_D_L_uni
#join -1 1 -2 1 <(sort -k1,1d Gene_sum_bfx656.hg38.e81) <(sort -k1,1d COUNTS/Featurecounter/bfx656/Gene_sum_mapped_sortedall) > Compare_D_L

