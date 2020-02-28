#!/usr/bin/env bash

in=$1
out=$2

if [ ! -f $out ]; then
	echo -e "SAMPLE\ttRNA\tCount" > $out
fi

fn=${in#*/}

zgrep -w CCA $in | perl -sae 'BEGIN{$c={}}{if($F[1] eq "seq"){@trnas=split(",",$F[6]);foreach $rna (@trnas){$c->{$rna}+=$F[3]}}}END{foreach $rna (keys %{$c}){print "$e\t$rna\t$c->{$rna}\n"}}' -- -e=$fn >> $out

#awk -v e="$a" 'BEGIN{FS="\t";OFS="";c=0}{c+=$3}END{if(e != "A"){print c,"\t"}else{print c}}' >> $out
