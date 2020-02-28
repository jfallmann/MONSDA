#!/bin/bash
#iabulk_rnaplex args: path_fasta_directory fst_fasta_index snd_fasta_index path_accessibility_directory path_output_directory
echo "iabulk_rnaplex $1 $2 $3 $4 $5";
#set l to length of shorter seq
#RNAplfold only up to u=13, so l > 13 does not make any sense
#SEQLENGTH1=$((`tail -n 1 ${3}SplitFasta${1}.fa | wc -m`-1)); 
#SEQLENGTH2=$((`tail -n 1 ${3}SplitFasta${2}.fa | wc -m`-1));
#LENGTH=$((SEQLENGTH1 < SEQLENGTH2 ? SEQLENGTH1 : SEQLENGTH2));
LENGTH=13;
#if [ $LENGTH -lt 40 ];
#then
  gunzip -c ${4}${1}_openen.gz > ${4}${1}_openen && RNAplex -l $LENGTH -q  ${3}SplitFasta${1}.fa -t ${3}SplitFasta${2}.fa -a $4 > $5$2.plex && rm -f ${4}${1}_openen
#else
#  RNAplex -l 40 -q ${3}SplitFasta${1}.fa -t ${3}SplitFasta${2}.fa -a $4 > $5$2.plex
#fi
