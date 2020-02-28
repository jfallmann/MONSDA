#!/usr/bin/env bash

in=$1
out=$2

if [ ! -f $out ]; then
	echo -ne "SAMPLE\tCCA" > $out
fi

fn=${in#*/}

echo -ne "\n$fn\t" >> $out

for a in CCA;do
	zcat $in|perl -sae 'BEGIN{$c=0}{if($F[0] eq $e){$c+=$F[1]}}END{{print $c}}' -- -e=$a >> $out
done

#zcat $fn|perl -sae 'print $fn."\t".join("\t",@F)' -- -e=$fn >> $out
#echo -ne "\n" >> $out
