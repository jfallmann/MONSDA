#!/usr/bin/env perl
#Script AnotateBed.pl;
#Last changed Time-stamp: <2020-03-12 16:53:46 fall> by joerg

#### use things ###
use strict;
#use warnings;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Pod::Usage;
use PerlIO::gzip;
use Set::IntervalTree;
use Cwd;
use Data::Dumper;
### use own modules
use FindBin::Real qw(Bin); # locate this script
use autodie; 
#use lib Bin() . "/../lib";
#use Collection;

###############################################################################
# start reading command line options
###############################################################################

my ( $species, $dir, $odir, $region, $spec, $bed, $VERBOSE, $eversion, $pversion, $dversion, $mysql, $dbversion, $file, $coords, $ignore, $enslibpath, $libpath, $specific, $usename );
my ( @dirs, @beds );

BEGIN{
	$VERBOSE=0;
	$eversion = "V79";
	$pversion = "1.2.3";
	$dversion = "79";

	pod2usage(-verbose => 0)
		unless GetOptions(
			"species|s=s"      => \$species,
			"region|r=s"       => \$region,
			"dir|d=s"          => \$dir,
			"odir|o=s"         => \$odir,
			"bedfile|b=s"      => \$bed,
			"mysql|y=s"        => \$mysql,
			"dbversion|v=s"	   => \$dbversion,
			"eversion|e=s"	   => \$eversion,
			"pversion|p=s"	   => \$pversion,
			"dversion|n=s"	   => \$dversion,
			"annotation|a=s"   => \$file,
			"enslibpath|x=s"   => \$enslibpath,
			"libpath|l=s"	   => \$libpath,
			"withcoords|w=s"   => \$coords,
			"ignorestrand|i=s" => \$ignore,
			"specific=s"       => \$specific,
			"usename|u"        => \$usename,
			"help|h"           => sub{pod2usage(-verbose => 1)},
			"man|m"            => sub{pod2usage(-verbose => 2)},
			"verbose"          => sub{ $VERBOSE++ }
		);

	$libpath = Bin() . "/../lib" unless(defined $libpath);
	print STDERR "Libpath: $libpath\n";
}

use lib "$libpath";
use Collection;
if (defined $enslibpath && $enslibpath ne ''){
		### use own modules relative to current location
		$enslibpath = Collection::parse_libpath($ENV{PWD},$enslibpath);
}
else{
	$enslibpath = Bin() . "/../lib";
}
print STDERR "ENSLIBPATH:".$enslibpath."\n";

#### use Bio things, check version of git repos to fit mysql if using ensembl-mysql for annotation ###
use lib $enslibpath . "/ensembl/modules";
use lib $enslibpath . "/ensembl-compara/modules";
use lib $enslibpath . "/ensembl-variation/modules";
use lib $enslibpath . "/bioperl-live";
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::LocatableSeq;
use Bio::EnsEMBL::Registry;
use Bio::SeqFeature::Generic;

##################
# Main
##################

if(!$dir){$dir	 = cwd();}
if(!$odir){$odir = $dir."\/Annotated";}

$file or die "No ensembl annotation gff given!\n";
$bed or die "No .bed file for processing chosen!\n";

###mysql stuff if wanted
die 'If you want to use a mysql-db please specify species! Optional also specify db version and mysql-host!\n' if ($file eq 'mysql' and !defined $species);

if(!$species){$species = "Human";}
if(!$region){$region   = "Core";}

###process call
my $pid = $$;
(my $job = `cat /proc/$pid/cmdline`)=~ s/\0/ /g;
print STDERR $job,"\n";
print STDERR "Species:$species, Dir:$dir, DB: $species\_$dversion\n" if ($file eq 'mysql');

#if (!-d $odir){
#    Collection::mkdircheck($odir);
#}
chdir ($dir) or die "Directory $dir could not be found!\n";

#Read beds;
push @beds,split(',',$bed);

##generate hash of unique bed entries
#my $unique={};
##Make chromosomes list
#my $chl = [];

#parse bedfile
print STDERR "Parsing Bed\n";
$ignore = '' unless (defined $ignore);
my ($unique,$chl) = Collection::parse_bed(\@beds,'','',$ignore);#, $unique, $chl);

#finish chromlist
my @chromlist =  sort { $a <=> $b } keys %{$chl};
#map{ $_->[0] } sort{ $a->[1] cmp $b->[1] } map{ [$_, -s $_] } @{$chl};

#Parse gff, define what we want to use for annotation here

my @checks;
push @checks, $specific if (defined $specific);
@checks = qw (exon intron gene transcript snoRNA snoRNA_gene snRNA snRNA_gene rRNA rRNA_gene CDS five_prime_UTR three_prime_UTR lincRNA lincRNA_gene miRNA miRNA_gene mRNA mt_gene NMD_transcript_variant processed_pseudogene processed_transcript pseudogene pseudogenic_transcript) unless (@checks and $specific);
my $check = {map { lc($_) => 1; } (@checks)};
@checks=();

print STDERR "Parsing Annotation\n";
my $intervals = Collection::parse_annotation($file, $check, $chl, $species, $mysql, $dversion, $usename);

$chl={};
#print Dumper($intervals);
#print STDERR Collection::get_features_interval($intervals->{1}->{"+"}->fetch(14404,29570));

print STDERR "Annotating\n";

foreach my $chroms (@chromlist){
	my $pat = quotemeta("\^$chroms\_");
	foreach my $uni ( sort{$a cmp $b} grep {${pat}} keys %{$unique} ){
		my @tempuni = split(/\_/,$uni);
		push @tempuni , split(/\_/,$unique->{$uni});

		(my $chromosome = $tempuni[0]) =~ s/=/_/g;
		my $start      = $tempuni[1];
		my $end        = $tempuni[2];
		my $strand     = $tempuni[3];
		(my $name      = $tempuni[4]) =~ s/=/_/g;
		(my $score      = $tempuni[5]) =~ s/u/0/g;
		my $summit     = $tempuni[6];
		(my $rest      = $tempuni[7]) =~ s/=/_/g;

		my @intersects;
		my $res;
		$res=$intervals->{"$chromosome"}->{"$strand"}->fetch($start,$end) if (defined $intervals->{"$chromosome"}->{"$strand"});
		push @intersects , $res if (defined $res);

		if (($strand eq "u" and !defined $intervals->{"$chromosome"}->{"u"}) or ($strand eq "." and !defined $intervals->{"$chromosome"}->{"."})){
			$res=$intervals->{"$chromosome"}->{"+"}->fetch($start,$end) if (defined $intervals->{"$chromosome"}->{"+"});
			push @intersects , @$res if (defined $res);
			$res=$intervals->{"$chromosome"}->{"-"}->fetch($start,$end) if (defined $intervals->{"$chromosome"}->{"-"});
			push @intersects , @$res if (defined $res);
		}

		my $tempbed = "$chromosome\t".$start."\t".$end;

		if (@intersects){
			my $fetch_coords = $coords || 0;
			my ($fstart,$fend,$ftype,$fid,$fstrand) = Collection::get_features_interval( \@intersects, $fetch_coords );
			if ($rest eq 'undef'){
				$rest = $name;
			}
			else{
				($rest = $name."\t".$rest) =~ s/=/_/g;
			}
			push @$fid, $name unless (defined @{$fid}[0]);
			if (defined $coords){
				for my $id (0..$#{$fid}){
					my $fstr = join(/,/,@$fstrand);
					$fstr = '.' if ($fstr =~ /u/ || ($fstr =~ /\-/ && $fstr =~ /\+/) || $fstr =~ /undef/ || $fstr =~ /\|/ || $fstr eq '');
					(my $featurename = @{$fid}[$id]) =~ s/=/_/g;
					my $featurestrand = @{$fstrand}[$id];
					$featurestrand = '.' if ($featurestrand =~ /u/ || ($featurestrand =~ /\-/ && $featurestrand =~ /\+/) || $featurestrand =~ /undef/ || $featurestrand =~ /\|/ || $featurestrand eq '');
					$tempbed.="\t$featurename\t$score\t$featurestrand";
					$tempbed .= "\t$rest";
					$tempbed .= "\t".join("\t",$featurename,@{$fstart}[$id],@{$fend}[$id],@{$ftype}[$id],@{$fstrand}[$id]) if defined ($fstart);
					print STDOUT $tempbed."\n" if (defined $tempbed);
					$tempbed = "$chromosome\t".$start."\t".$end;
				}
			}
			else{
				my $fstr = join(/,/,@$fstrand);
				$fstr = '.' if ($fstr =~ /u/ || ($fstr =~ /\-/ && $fstr =~ /\+/) || $fstr =~ /undef/ || $fstr =~ /\|/ || $fstr eq '');
				(my $featurename = join("|",grep(/ENS*G/,@{$fid}))) =~ s/=/_/g;
				$featurename = 'intergenic' if ($featurename eq '');
				$tempbed.="\t$featurename\t$score\t$fstr";
				$tempbed .= "\t$rest";
				$tempbed .= "\t".join("\t",join("|",@{$fid}),join("|",@{$ftype}));
				print STDOUT $tempbed."\n" if (defined $tempbed);
				$tempbed = "$chromosome\t".$start."\t".$end;
			}
		}
		else{
			$strand = '.' if ($strand =~ /u/ || ($strand =~ /\-/ && $strand =~ /\+/) || $strand =~ /undef/ || $strand eq '');
			$tempbed.="\t$name\t$score\t$strand";
			$tempbed .= "\t$rest" unless ($rest eq 'undef');
			print STDOUT $tempbed."\n" if (defined $tempbed);
		}
	}
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
