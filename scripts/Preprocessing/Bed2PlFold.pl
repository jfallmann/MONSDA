#!/usr/bin/env perl
#Script Bed2PlFold.pl;
#Last changed Time-stamp: <2017-05-04 16:42:46 fall> by Joerg Fallmann <fallmann.joerg@gmail.com>

##################
# Libs
##################
use strict;
use warnings;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Pod::Usage;
use Cwd;
use PerlIO::gzip;
use Data::Dumper;
### use own modules
use FindBin::Real qw(Bin);
use lib Bin() . "/../lib";
use Collection;

#### use Bio things, check version of git repos to fit mysql if using ensembl-mysql for annotation ###
use lib Bin() . "/../lib/ensembl/modules";
use lib Bin() . "/../lib/ensembl-compara/modules";
use lib Bin() . "/../lib/bioperl-live";
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::LocatableSeq;
use Bio::EnsEMBL::Registry;
use Bio::SeqFeature::Generic;

##################
# Options 
##################

my ($VERBOSE, $dir, $odir, $pdir, $bed, $ustart, $uend, $block, $rel, $coords, $species, $region, $eversion, $pversion, $dversion, $mysql, $dbversion, $highmem);#, $wobble, $minlength);

BEGIN{
	$VERBOSE=0;
	$eversion = "V79";
	$pversion = "1.2.3";
	$dversion = "79";

	pod2usage(-verbose => 0)
		unless GetOptions(
			"dir|d=s"				=> \$dir,     #Dir for Bed
			"odir|o=s"				=> \$odir,    #Dir for Output
			"plfold|f=s"			=> \$pdir,    #Dir for plfold output
			"ustart|u=s"			=> \$ustart,  #U range from
			"uend|t=s"				=> \$uend,    #U range to	    
			"bed|b=s"				=> \$bed,     #Bed filename
			"block|k=s"				=> \$block,   #Parse Bed12
			#	    "wobble|w=s"    => \$wobble,  #Embedding region ### NOT FUNCTIONAL PEAK MERGE FIRST OTHERWISE OVERLAP
			#	    "minlength|l=s" => \$minlength,  #Minimum length of peak ### NOT YET NECESSARY
			"rel|r=s"				=> \$rel,     #Peak already in relative coordinates
			"coords|c=s"			=> \$coords,  #Gene coordinates for relative peaks
			"highmem|x=s"			=> \$highmem,  #Enable high mem usage for faster processing
			"species|s=s"           => \$species,
			"region|r=s"            => \$region,
			"mysql|y=s"				=> \$mysql,
			"dbversion|v=s"			=> \$dbversion,
			"eversion|e=s"			=> \$eversion,
			"pversion|p=s"			=> \$pversion,
			"dversion|n=s"			=> \$dversion,
			"help|h"				=> sub{pod2usage(-verbose => 1)},
			"man|m"					=> sub{pod2usage(-verbose => 2)},      
			"verbose"				=> sub{ $VERBOSE++ }
		);
}

##################
# Main
##################

$bed or die "No .bed file for processing chosen!\n";

###process call
my $pid = $$;
(my $job = `cat /proc/$pid/cmdline`)=~ s/\0/ /g;
print STDERR $job,"\n";

#$wobble = 19 unless ($wobble);

$dir  =	 cwd() unless ($dir);
chdir ("$dir") || die ("Could not find directory $dir!\n");

my $gencoords;

unless (defined $rel){
	### Parse Coordinate File, if none is existent, create new one

	if (defined $coords && -e $coords){
		open(REL,"<:gzip(autopop)",$coords) || die ("Could not open $coords!\n");
		while(<REL>){
			chomp (my $raw = $_);
			my ($index, $gene, $coords) = split (/\t/,$raw);
			my ($chrom, $pos, $strand ) = split (/\:/,$coords);
			my ($start, $end) = split(/\-/,$pos);
			$strand =~ s/\"//g;
			$gencoords->{$gene} = join(":",$chrom,$start-1,$end,$strand);
		}
		close(REL);
	}
	else{
		$gencoords = Collection::fetch_gene_coords($species, $mysql, $dversion, $coords);
	}
}

my $genes = {};
$genes = Collection::parse_bed_plfold($bed, $block, $rel, $gencoords);

### Parse RNAplfold Output
my $RT = (-1.9872041*10**(-3))*(37+273.15);
die ("Could not find directory $pdir!\n") unless (-d $pdir || -l $pdir);
my $folded = Collection::parse_plfold($genes, $RT, $pdir, $ustart, $uend, $highmem);

### Generate Output

$odir = "Parsed_$pdir" unless (defined $odir);

if (!-d $odir){
    print STDERR "Will create dir $odir!\n";
    Collection::mkdircheck($odir);
}

foreach my $name (keys %{$folded}){
	print STDERR "Printing output for $name\n";
	open (OL,">:gzip",$odir."\/".$name."_plfold_parsed.bed.gz") or die "Could not open ".$odir."\/".$name."_plfold_parsed.bed.gz";
#	print OL join("\t","Chromosome","Start","End","Peak","Position","Strand",($ustart .. $uend))."\n";
	foreach my $pos (sort {$a<=>$b} keys %{$folded->{$name}->{"Peak"}}){		
		if (defined $folded->{$name}->{"Fold"}->{$pos}){
			print OL $folded->{$name}->{"Peak"}->{$pos}."\t".$folded->{$name}->{"Fold"}->{$pos}."\n";
		}
		else{
			warn "Could not find information for Position $pos in $name in plfold output!\n" unless (defined $folded->{$name}->{"Fold"}->{$pos});			
		}
	}
	close(OL);
}
