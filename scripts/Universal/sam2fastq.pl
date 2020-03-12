#!/usr/bin/perl
use strict;
use warnings;
use autodie;

my $in1 = shift; 	#bam file
my $in2 = shift;	#fastq

open IN1, "samtools view -h $in1 |";
#open IN2, "< $in2" or die "can t open $in2\n";
open(IN2, "gunzip -c $in2 |") or die "gunzip $in2: $!";


my %hash =();
while(<IN1>){
	next if($_ =~ /^@/);
	chomp $_;

	my @line = split(/\t/,$_);
	my $ID = "@".$line[0];
	$hash{$ID} = 0;
}


my @entry = ();

while(<IN2>){
	 chomp;
         push @entry, $_;

         if (scalar(@entry) == 4) {

                my ($id, $seq, $plusLine, $qual) = @entry;
                @entry = ();

		if(exists $hash{$id}){
			print join ("\n", $id, $seq, $plusLine, $qual)."\n";
		}
	}
}
