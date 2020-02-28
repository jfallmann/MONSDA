#!/bin/env bash
COND=$1
SOURCE=$2
PATTERN=$3
SAMPLES=$4

#cut -d$'\t' -f1 L12270_09_WT_d4_1_2_IN_R1_mapped_sorted.transcript_counts > tmp_counts

for i in *$COND\_*$SOURCE\_*$PATTERN\_unique/genes.fpkm_tracking;do
	echo $i
	tail -n+2 $i|cut -d$'\t' -f1,5,10  >>  GENECOUNTS_$COND\_$SOURCE\_tmp
done

#env sa=$SAMPLES perl -lan -F'\t' -e 'BEGIN{%exp=()}; $exp{$F[0]}+=$F[1]; END{foreach $key(keys %exp){print $key,"\t",$exp{$key}/$ENV{sa}}}' COUNTS_$COND\_$SOURCE\_tmp |sort -k1,1d > COUNTS_$COND\_$SOURCE.transcript.fpkm

#rm -f COUNTS_$COND\_$SOURCE\_tmp

#join -1 1 -2 1 <(sort -k1,1d Gene_sum_bfx656.hg38.e81) <(sort -k1,1d COUNTS/Featurecounter/bfx656/Gene_sum_mapped_sorted_unique.counts) > Compare_D_L_uni
#join -1 1 -2 1 <(sort -k1,1d Gene_sum_bfx656.hg38.e81) <(sort -k1,1d COUNTS/Featurecounter/bfx656/Gene_sum_mapped_sortedall) > Compare_D_L

#bash ../../../Workflows/scripts/Shells/MergeExpression_Cufflinks.sh KO_d4 IP R1 3
#bash ../../../Workflows/scripts/Shells/MergeExpression_Cufflinks.sh KO_d0 IP R1 3
#bash ../../../Workflows/scripts/Shells/MergeExpression_Cufflinks.sh WT_d4 IP R1 3
#bash ../../../Workflows/scripts/Shells/MergeExpression_Cufflinks.sh WT_d0 IP R1 3
#bash ../../../Workflows/scripts/Shells/MergeExpression_Cufflinks.sh WT_d0 IN R1 3
#bash ../../../Workflows/scripts/Shells/MergeExpression_Cufflinks.sh WT_d4 IN R1 3
#bash ../../../Workflows/scripts/Shells/MergeExpression_Cufflinks.sh KO_d4 IN R1 3
#bash ../../../Workflows/scripts/Shells/MergeExpression_Cufflinks.sh KO_d0 IN R1 3
