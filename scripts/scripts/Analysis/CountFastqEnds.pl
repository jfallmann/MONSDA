#!/usr/bin/perl
use strict;
use warnings;
use PerlIO::gzip;

my $in = shift; 	#fastq

open (SEQ, "<:gzip(autopop)", $in) or die "$!";

my $tails = ();
my $all = 0;
my $i = 1;

while(<SEQ>){

	chomp(my $line = $_);

	if ($i == 2){
		my $tail = substr $line,-3;
		$tails->{$tail}++;
		$i++;
	}
	elsif($i==4){
		$i=1;
		next;
	}
	else{
		$i++;
		next;
	}
}

foreach my $tail (sort{$a eq $b} keys %{$tails}){
	print join("\t", $tail, $tails->{$tail})."\n";
}
