#!/usr/bin/env bash

in=$1
out=$2

if [ ! -f $out ]; then
	echo -ne "SAMPLE\tCCACCA\tCCACC\tCCAC\tCCA\tCC\tCA\tC\tA" > $out
fi

fn=${in#*/}

echo -ne "\n$fn\t" >> $out
#perl -se '{$o = (split("\/",$f))[-1];print "\n".$o."\t"}' -- -f=$in >> $out

for a in CCACCA CCACC CCAC CCA CC CA C A;do
	zgrep -w $a $in | perl -sae 'BEGIN{$c=0}{if($F[1] eq "seq"){$c+=$F[3]}}END{if($e ne "A"){print $c."\t"}else{print $c}}' -- -e=$a >> $out
done

#awk -v e="$a" 'BEGIN{FS="\t";OFS="";c=0}{c+=$3}END{if(e != "A"){print c,"\t"}else{print c}}' >> $out
