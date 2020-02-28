#!/usr/bin/env perl
#Script Link_RNAplfold.pl;
#Last changed Time-stamp: <2017-05-04 12:15:17 fall> by joerg

##################
# Libs
##################
use strict;
use warnings;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Pod::Usage;
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
my ($VERBOSE, $species, $dir, $pdir, $odir, $region, $spec, $id, $table, $filename, $eversion, $pversion, $dversion, $mysqlhost, $genelist, $setting);

BEGIN{
	$VERBOSE=0;
	$eversion = "V79";
	$pversion = "1.2.3";
	$dversion = "79";
	
	pod2usage(-verbose => 0)
		unless GetOptions(
			"species|s=s"	=> \$species,
			"region|r=s"	=> \$region,
			"dir|d=s"		=> \$dir,
			"plfold|f=s"	=> \$pdir,    #Dir for plfold output
			"odir|o=s"		=> \$odir,
			"table|t=s"		=> \$table,
			"id|i=s"        => \$id,
			"eversion|e=s"	=> \$eversion,
			"pversion|p=s"	=> \$pversion,
			"dversion|b=s"	=> \$dversion,
			"mysqlhost|y=s" => \$mysqlhost,
			"help|h"		=> sub{pod2usage(-verbose => 1)},
			"man|m"			=> sub{pod2usage(-verbose => 2)},      
			"verbose"		=> sub{ $VERBOSE++ }
		);
}



die 'Please provide a list of identifiers to link, either as file or comma separated commandline option!\n' unless ($id);
my $pid = $$;
(my $job = `cat /proc/$pid/cmdline`)=~ s/\0/ /g;
print STDERR $job,"Xref:$id\n";

##################
# Main
##################

die 'Please provide a lookup-table for the gene ids, e.g. Homo_sapiens_lookup_table.gz!\n' unless $table;

### Parsing hashtable for id
open (T, "<:gzip(autopop)",$table) || die ("Could not open $table!");

my %generef;
while (<T>){
    chomp $_;
    my ($index, $ensid, $pos) = split("\t",$_);
    $generef{$ensid} = $index;
}
close (T);

### Get list of genes
my @genes;
if (-e $id){
	open(L,"<:gzip(autopop)",$id) or die "$!\n";
	if ($id =~ /.bed[.gz]*$/){
		@genes = @{Collection::parse_bed_linking($id)};
	}
	else{
		while (<L>){
			chomp $_;
			push @genes, $_;
		}
	}
	close(L);
}
else{
	@genes = split(",",$id);
}

### Link plfold output

$odir = "PlFold_Linked" unless ($odir);

if (!-d $odir){
    print STDERR "Will create dir $odir!\n";
    Collection::mkdircheck($odir);
}

$setting = "240_160" unless $setting;

if ($species eq "Human" or $species eq "human" or $species eq "homo_sapiens"){$spec							= "Homo_sapiens";}
elsif ($species eq "Mouse" or $species eq "mouse" or $species eq "mus_musculus"){$spec						= "Mus_musculus";}
elsif ($species eq "Zebrafish" or $species eq "zebrafish" or $species eq "danio_rerio"){$spec				= "Danio_rerio";}
elsif ($species eq "Drosophila" or $species eq "drosophila" or $species eq "drosophila_melanogaster"){$spec = "Drosophila_melanogaster";}
elsif ($species eq "Worm" or $species eq "worm" or $species eq "caenorhabditis_elegans"){$spec				= "Caenorhabditis_elegans";}
else{$spec																									= $species;}

foreach my $gene (@genes){
	`ln -s $pdir\/RNAplfold_$spec\/$spec\_RNAplfold_$setting\/$generef{$gene}\_lunp.gz $odir\/$gene\_lunp.gz`;
}

############################## SUBS ##############################


############################## POD ##############################

=pod

=head1 NAME

FromXtoG.pl - Gives Gene information for Xref

=head1 SYNOPSIS

FromXtoG.pl [-file I<STRING>] [-species I<STRING>] [-man] [-help]

=head1 OPTIONS

=over 4

=item B<-help>

Show synopsis.

=item B<-man>

Show man page.

=item B<-file> I<STRING>

Give Xref containing file like: perl FromXtoG.pl -f Xref.csv

=item B<-species> I<STRING>

Search for species related Xrefs in given file,e.g. perl FromXtoF.pl -f Xref.csv -s human

=back
    
=head1 DESCRIPTION

This Program will search for the matching gene for a given Ensembl Xref
and print its stable_ID to make it easier finding your gene of interest.

=head1 AUTHOR

Joerg Fallmann
    
=cut

__END__
