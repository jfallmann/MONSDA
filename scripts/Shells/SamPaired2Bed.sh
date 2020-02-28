#!/bin/bash
FILE=$1
OUT=$(basename $FILE)

if [[ "$FILE" == *.gz* ]]
then
    zcat $FILE|grep "HWI-ST"|awk '{FS="\t";OFS="\t"}{if(and($2,0x2)){if(and($2,0x40)>0 && and($2,0x10)>0){print $3,$4-1,$4+length($10),$1,$2,"-"} else if(and($2,0x40)>0 && and($2,0x10)==0){print $3,$4-1,$4+length($10),$1,$2,"+"} else if(and($2,0x80)>0 && and($2,0x20)>0){print $3,$4-1,$4+length($10),$1,$2,"-"} else if(and($2,0x80)>0 && and($2,0x20)==0){print $3,$4-1,$4+length($10),$1,$2,"+"}}}' - > $OUT"_Unique.bed";
else
    cat $FILE|grep "HWI-ST"|awk '{FS="\t";OFS="\t"}{if(and($2,0x2)){if(and($2,0x40)>0 && and($2,0x10)>0){print $3,$4-1,$4+length($10),$1,$2,"-"} else if(and($2,0x40)>0 && and($2,0x10)==0){print $3,$4-1,$4+length($10),$1,$2,"+"} else if(and($2,0x80)>0 && and($2,0x20)>0){print $3,$4-1,$4+length($10),$1,$2,"-"} else if(and($2,0x80)>0 && and($2,0x20)==0){print $3,$4-1,$4+length($10),$1,$2,"+"}}}' - > $OUT"_Unique.bed";
fi
#awk '{FS="\t";OFS="\t"}{if(and($2,0x1)){if(and($2,0x40)>0 && and($2,0x10)>0){print $3,$4-1,$4+length($10)-1,$1,$2,"-"} else if(and($2,0x40)>0 && and($2,0x10)==0){print $3,$4-1,$4+length($10)-1,$1,$2,"+"} else if(and($2,0x80)>0 && and($2,0x10)>0){print $3,$4-1,$4+length($10)-1,$1,$2,"-"} else if(and($2,0x80)>0 && and($2,0x10)==0){print $3,$4-1,$4+length($10)-1,$1,$2,"+"}}}' - > $OUT"_Unique.bed";
#awk '{FS="\t";OFS="\t"}{if(and($2,1)){if(and($2,64)>0 && and($2,16)>0){print $3,$4-1,$4+length($10)-1,$1,$2,"-"} else if(and($2,64)>0 && and($2,16)==0){print $3,$4-1,$4+length($10)-1,$1,$2,"+"} else if(and($2,128)>0 && and($2,32)>0){print $3,$4-1,$4+length($10)-1,$1,$2,"+"} else if(and($2,128)>0 && and($2,32)==0){print $3,$4-1,$4+length($10)-1,$1,$2,"-"}}}' - > $OUT"_Unique.bed";
#awk '{FS="\t";OFS="\t"}{if(and($2,1)){if(and($2,16)){print $3,$4-1,$4+length($10)-1,$1,$2,"-"} else if($1 ~ /\#0\/1/){print $3,$4-1,$4+length($10)-1,$1,$2,"+"} else if(and($2,32)){print $3,$4-1,$4+length($10)-1,$1,$2,"+"} else if($1 ~ /\#0\/2/){print $3,$4-1,$4+length($10)-1,$1,$2,"-"}}}' - > $OUT"_Unique.bed";
