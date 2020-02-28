#!/usr/bin/env perl
#Script NormalizePeaks.pl;
#Last changed Time-stamp: <2018-05-02 14:22:04 fall> by joerg

#### use things ###
use strict;
use warnings;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Pod::Usage;
use PerlIO::gzip;
use Cwd;
use Data::Dumper;
use File::Slurp;
use Math::Round;
### use own modules
use FindBin::Real qw(Bin); # locate this script
#use lib Bin() . "/../lib";
#use Collection;

###############################################################################
# start reading command line options
###############################################################################

my ( $dir, $odir, $VERBOSE, $peak, $expression, $program, $cutoff, $unit, $alpha, $libpath );

BEGIN{
	$VERBOSE=0;
	pod2usage(-verbose => 0)
		unless GetOptions(
			"dir|d=s"        => \$dir,
			"odir|o=s"       => \$odir,
			"peakfile|b=s"   => \$peak,
			"expression|e=s" => \$expression,
			"program|p=s"    => \$program,
			"cutoff|c=s"     => \$cutoff,
			"unit|u=s"       => \$unit,
			"alpha|a=s"      => \$alpha,
			"libpath|l=s"    => \$libpath,
			#			"withcoords|w=s" => \$coords,	
			"help|h"         => sub{pod2usage(-verbose => 1)},
			"man|m"          => sub{pod2usage(-verbose => 2)},      
			"verbose"        => sub{ $VERBOSE++ }
		);
	$libpath = Bin() . "/../lib" unless(defined $libpath);
}

use lib "$libpath";
use Collection;	

##################
# Main
##################

if(!$dir){$dir	 = cwd();}
if(!$odir){$odir = $dir."\/Normalized";}

$program = 'cufflinks' unless (defined $program);
$cutoff = 1 unless (defined $cutoff);
warn "You chose an expression cutoff below 1, this means that normalization will increase the signal on transcripts with very low expression, are you sure you want this?!\n" if ($cutoff <1) ;
$expression or die "No expression file given!\n";
$peak or die "No peakfile for processing chosen!\n";
$alpha = 1 unless (defined $alpha);
###process call
my $pid = $$;
(my $job = `cat /proc/$pid/cmdline`)=~ s/\0/ /g;
print STDERR $job,"\n";
print STDERR "Program: $program, Cutoff: $cutoff\n";

chdir ($dir) or die "Directory $dir could not be found!\n";

#Read peaks;
push my @beds,split(',',$peak);

##generate hash of unique bed entries
#parse bedfile
my ($unique) = Collection::parse_bed(\@beds);

#Parse expression file, define what we see as expressed (cutoff)

my @exp = split(/ /,$expression);

my $expressed = Collection::parse_expression($expression, $cutoff, $program, $unit);
#print Dumper($expressed);
### Normalize peaks

foreach my $peak ( keys %{$unique} ){
	my @tempuni = split(/\_/,$peak);
	push @tempuni , split(/\_/,$unique->{$peak});
	my $chromosome = $tempuni[0];		
	my $start      = $tempuni[1];
	my $end        = $tempuni[2];
	my $strand     = $tempuni[3];
	(my $name       = $tempuni[4])=~s/=/_/g;
	my $score      = $tempuni[5];
	my $summit     = $tempuni[6];
	(my $rest       = $tempuni[7])=~s/=/_/g;

	$strand = '.' if ($strand =~ /u/ || ($strand =~ /\-/ && $strand =~ /\+/) || $strand =~ /undef/);
	
	my @annotation = split(/\|/,join("\|",split("\t",$rest)));
	my $EX=0;
	foreach my $id ( @annotation ){
#		print STDERR $id,"\t",$expressed->{$id},"\n" if ($id =~ /^ENS[a-zA-Z]*T/ && defined $expressed->{$id} && $expressed->{$id} >= $cutoff);
		$EX += $expressed->{$id} if ($id =~ /^ENS[a-zA-Z]*T/ && defined $expressed->{$id} && $expressed->{$id} >= $cutoff);
	}
	next if ($EX == 0);
	my $normscore = nearest(.0001,$score/($EX+$alpha));

	#normalize profile
	my @restall = split("\t",$rest);
	my @profile = @{Collection::parse_peakprofile([split(/\|/,$restall[0])],$start)};
	my $poscov = {};
	foreach my $pos (@profile){
		my @bla = split (/\:/,$pos);
		$bla[1] = nearest(.0001,$bla[1]/($EX+$alpha));
		$poscov->{$bla[0]} = $bla[1];
	}
	my ($normprof, $normarea) = Collection::generate_peakprofile($start, $end, $poscov);
	my $normprofile = substr($normprof,0,-1);
	$restall[0] = $normprofile;
	$restall[2] = $normarea;
	$rest = join("\t",@restall);
	
	print STDOUT join("\t","chr".$chromosome, $start, $end, $name, $normscore, $strand, $rest."\n");
	
}

#############END####################

###############################################################################
# POD
###############################################################################

=pod

=head1 NAME

foldem.pl - Retrieve subsequence of Genbank File and fold it with RNAfold

=head1 SYNOPSIS

foldem.pl [-dir I<STRING>] [-positions I<ARRAY>] [-transcript I<STRING>] [-man] [-help]

=head1 OPTIONS

=over 4

=item B<--help|-h>

Show synopsis.

=item B<--man|-m>

Show man page.

=item B<--dir|-d> I<STRING>

Select working directory, e.g. foldem.pl --dir Human_ATTTA/ I<STRING>

=item B<--positions|-p> I<STRING>

Enter positions (start end) of ARE motifs in sequence, e.g. foldem.pl --positions 1200 1204 I<STRING>

=item B<--gene|-g> I<STRING>

Enter the name of the gene to work with, e.g. foldem.pl --gene ENSG00000000003 I<STRING>

=back
    
=head1 DESCRIPTION

This program will take positions of ARE motifs return the sequence with 100 flanking nucleotides on both sides.

=head1 AUTHOR

Joerg Fallmann
    
=cut
