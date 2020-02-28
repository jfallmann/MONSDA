#!/bin/bash
FILE=$1
OUT=$(basename $FILE)

if [[ "$FILE" == *.gz* ]]
then
    zcat $FILE|grep "HWI-ST"|awk '{FS="\t";OFS="\t"}{if(and($2,16)){print $3,$4-1,$4+length($10),$1,$2,"-"} else {print $3,$4-1,$4+length($10),$1,$2,"+"}}' - > $OUT"_Unique.bed";
else
    cat $FILE|grep "HWI-ST"|awk '{FS="\t";OFS="\t"}{if(and($2,16)){print $3,$4-1,$4+length($10),$1,$2,"-"} else {print $3,$4-1,$4+length($10),$1,$2,"+"}}' - > $OUT"_Unique.bed";
fi

