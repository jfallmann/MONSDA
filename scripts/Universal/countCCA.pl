#!/usr/bin/perl
use strict;
use warnings;
use PerlIO::gzip;

my $in = shift; 	#bam file        
my $in2 = shift;	#cluster sequences


open BAM, "samtools view $in |";
open SEQ, "<:gzip(autopop)", "$in2" or die "can t open $in2\n";


my %hash =();
my $c1 =0;
my $c2 =0;
my $c3 =0;
my $c4 =0;
my $c5 = 0;


while(<SEQ>){
	chomp $_;
	
	if($_ =~m /^>(.*)/){
		my $id = $1;
		my $seq = <SEQ>;
		chomp $seq;
		$hash{$id} = $seq;
	}
}

my $all;
while(<BAM>){
        $all++;
	chomp $_;
	
	my @line = split(/\t/,$_);
	
	my $cluster = $line[2];			#Reference sequence NAME
	my $start = $line[3];			#1-based leftmost mapping POSition	
	my $seq = $line[9]; 			#segment SEQuence (read)
	my $seqL = length($seq);
	
	my $cigar = $line[5];
	my @array = split(/(M|I|D|N|S|H|P|X)/,$cigar);
	
	#sum number of deletions or sub number of insertions from the read length
	my $del = 0;
	my $ins = 0;
	for(my $i=0; $i < scalar @array ; $i++){
		if($array[$i] eq "D"){
			$del += $array[$i-1];
		}
		if($array[$i] eq "I"){
                	$ins += $array[$i-1];
		}
	}
	
	$seqL += $del;
	$seqL -= $ins;
	
	my $len = length($hash{$cluster});	#length of the cluster seq with CCACCA
	my $tailCCA = substr $seq, -3;
        my $tailCC = substr $seq, -2;
        my $tailC = substr $seq, -1;
	my $tailCCACCA = substr $seq, -6;

	### only reads with CCA and mapping position at 3' end -3 (CCA) position
        if(($start -1 + $seqL) == ($len - 3) && $tailCCA eq "CCA"){
		$c1++;
        }
	### only reads with CCACCA and mapping position at 3' end position
        elsif(($start -1 + $seqL) == ($len) && $tailCCACCA eq "CCACCA"){
		$c2++;
        }
	### only reads with no CCA and mapping position at 3' end -6 (CCACCA) position
        elsif(($start -1 + $seqL) == ($len - 6)){
		$c3++;
	}
	### only reads with C and mapping position at 3' end -5 (CACCA) position
	elsif(($start -1 + $seqL) == ($len - 5) && $tailC eq "C"){
		$c4++;
        }
	### only reads with CC and mapping position at 3' end -5 (ACCA) position
	elsif(($start -1 + $seqL) == ($len - 4) && $tailCC eq "CC"){
		$c5++;
      	}
}


print join ("\t", "all", "CCACCA", "CCA", "CC", "C", "no");
print "\n";
print join ("\t", $all, $c2, $c1, $c5, $c4, $c3);
print "\n";

