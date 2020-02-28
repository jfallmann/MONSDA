#!/usr/bin/env bash

inf=$1
name=${inf%*.fa}
out=$2
out2=$3

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ${inf} |perl -lane '{if ($_ =~ /^>/){@heads=split(/\s/,$_);@head=();for $h (@heads){$h=~s/\W$//g;push @head, $h if ($h =~ />|chrM|>M|>MT|mito|plastid|plasmid|organel/);}print join("_",@head)}elsif($_ =~ /^$/){next;}else{print $_}}' > ${name}.genomic.fa #get rid of lines with spaces and too much text in header
#sed -i '/^\s*$/d' ${name}.fa #remove empty lines
touch ${out}
grep -A1 -E '>chrM|>M|>MT|mito|plastid|plasmid|organel' ${name}.genomic.fa > ${out}
if [[ -s ${out} ]];then
	echo "Found ${out}";
else
	echo -e ">Decoy_MT\naaaaaaaaaaaaaaaaaaaaaaccccccccccccccccccccccccccccttttttttttttttttttttttttttggggggggggggggggggggggggggg" > ${out}
fi
touch ${out2}
perl -ne 'chomp;if($_ =~ /^>/){print $_."\t"} else{print $_."\n"}' ${name}.genomic.fa | grep -v -E '>chrM|>M|>MT|mitochond|plastid|plasmid|organel' | sed 's/\t/\n/g' > ${out2}
