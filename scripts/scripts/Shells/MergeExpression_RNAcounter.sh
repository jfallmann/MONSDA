#!/bin/env bash
COND=$1
SOURCE=$2
PATTERN=$3
SAMPLES=$4

cut -d$'\t' -f1 L12270_09_WT_d4_1_2_IN_R1_mapped_sorted.transcript_counts > tmp_counts

for i in *$COND\_*$SOURCE\_*$PATTERN\.transcript_counts;do
	echo $i
	join -1 1 -2 1 <(sort -k1,1d tmp_counts) <(sort -k1,1d $i|cut -d$'\t' -f1,2,3) > COUNTS_$COND\_$SOURCE && cp -f COUNTS_$COND\_$SOURCE tmp_counts
done
rm -f tmp_counts

awk 'BEGIN{FS=" ";OFS=FS}{s=1;t=0;for(i=2;i<=NF;i+=2){if($i >=1){t+=$i;s++}} print $1,t/s;t=0;s=1}' COUNTS_$COND\_$SOURCE |sort -k1,1d > $SOURCE\_$COND\_$PATTERN\.transcript.counts

#join -1 1 -2 1 <(sort -k1,1d Gene_sum_bfx656.hg38.e81) <(sort -k1,1d COUNTS/Featurecounter/bfx656/Gene_sum_mapped_sorted_unique.counts) > Compare_D_L_uni
#join -1 1 -2 1 <(sort -k1,1d Gene_sum_bfx656.hg38.e81) <(sort -k1,1d COUNTS/Featurecounter/bfx656/Gene_sum_mapped_sortedall) > Compare_D_L

#bash ../../../Workflows/scripts/Shells/MergeExpression_RNAcounter.sh KO_d4 IP mapped_sorted 3
#bash ../../../Workflows/scripts/Shells/MergeExpression_RNAcounter.sh KO_d0 IP mapped_sorted 3
#bash ../../../Workflows/scripts/Shells/MergeExpression_RNAcounter.sh WT_d4 IP mapped_sorted 3
#bash ../../../Workflows/scripts/Shells/MergeExpression_RNAcounter.sh WT_d0 IP mapped_sorted 3
#bash ../../../Workflows/scripts/Shells/MergeExpression_RNAcounter.sh WT_d0 IN mapped_sorted 3
#bash ../../../Workflows/scripts/Shells/MergeExpression_RNAcounter.sh WT_d4 IN mapped_sorted 3
#bash ../../../Workflows/scripts/Shells/MergeExpression_RNAcounter.sh KO_d4 IN mapped_sorted 3
#bash ../../../Workflows/scripts/Shells/MergeExpression_RNAcounter.sh KO_d0 IN mapped_sorted 3
