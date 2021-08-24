#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use Getopt::Long;
use vars qw ($help $inlist $outfile $annotation);

$outfile="junction_table_dexgo";
GetOptions (
    "i=s"       => \$inlist,
    "h"         => \$help,
    "o=s"       => \$outfile,
    "a=s"       => \$annotation
    );
usage() if ($help || !$inlist);

open(my $IN, "<", $inlist);
my @files;
my @names;
my $count=0;

while(<$IN>) { 
	chomp(my $line = $_);
	next if ($line =~ /^#/);
    my @F = split("\t",$line);
    $names[$count]=$F[0] if $F[0];
    $files[$count]=$F[1] if $F[1];
    $count++;
}

die "$#names not $#files!" if ($#names != $#files or $#names <1 ) ;

#my $filelist=join(" ",@files);
my $headerlist=join("\t",@names);
close($IN);

#my %annotation;
#if($annotation) {
#    open($IN, "<",$annotation);
#}

open(my $OUT, ">>", $outfile);
print $OUT "junction\ttype\t$headerlist\tgeneID\tgeneName\n";

my $oldid="";
my $zaehl=0;
$count=0;

my @countfiles = @files;#split(" ",$filelist);

while (my $f = shift(@countfiles)){
	open($IN,'<',$f);
	while(<$IN>) {
		next if ($_ =~ /^#/);
		my @F=split;
		$F[0]=~/(\S+):\d+$/;
		my $geneid=$1;
		if($geneid ne $oldid) {
			$count++;
		}
		$oldid=$geneid;
		$zaehl+=100;
		my $z2=$zaehl+50;
		print $OUT "chrfoo:$zaehl"."-$z2\tN_w";
		for (my $i=1; $i<=$#F; $i+=2) {
			print $OUT "\t$F[$i]";
		}
		print $OUT "\t$geneid\tbar$count\n";
	}
	close($OUT);
}

printf STDERR "You now should hav a file $outfile to play with\nThank you for travelling with us, Good bye!\n";


sub usage {
  print STDERR "\nHTseq2DIEGO.pl\n";
  print STDERR "usage: HTseq2DIEGO.pl -i <file>  [OPTIONS]\n";
  print STDERR "\n";
  print STDERR "[INPUT]\n";
  print STDERR " -i <file>    file containing input files and ids\n \t\tid [tab] path.to/file\n";
  print STDERR " -o <file>    output file name (default:junction_table_dexdas )\n";
  print STDERR " -h <file>    this (usefull) help message\n";
  print STDERR "[VERSION]\n";
  print STDERR " 06-25-2012\n";
  print STDERR "[BUGS]\n";
  print STDERR " Please report bugs to salzamt\@bioinf.uni-leipzig.de\n";
  print STDERR "\n";
  exit(-1);
}
