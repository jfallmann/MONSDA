#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Path qw(make_path remove_tree);
use Cwd;
use PerlIO::gzip;
use Statistics::Basic qw(:all);
use Math::Round;

my $VERBOSE = 0;

# process command-line
pod2usage(-verbose => 0)
    unless GetOptions(
	    "dir|d=s"	   => \my $dir,
	    "file|f=s"	   => \my $file,
	    "output|o=s"   => \my $output,  
	    "nrg2prob|e=s" => \my $nrg,
	    "prob2nrg|p=s" => \my $prob,
	    "help"	   => sub{pod2usage(-verbose => 1)},
	    "man"	   => sub{pod2usage(-verbose => 2)},      
	    "verbose"	   => sub{ $VERBOSE++ }
    );

$dir  =	 cwd() unless ($dir);
($file && $output) || die ("$!, Give in and out file name!\n");
chdir ("$dir");

my $line;
my @data;

open (F, "<:gzip(autopop)",$file) || die ("Could not open $file!");
while ($line = <F>){
	chomp $line;
	push @data, $line;
}
close (F);

if(-e $output){
	rename ("$output", "$output\.bak");
	print STDERR "$output exists, renamed to $output.bak\n";
}

my $RT = (-1.9872041*10**(-3))*(37+273.15);

open (OUT,">:gzip","$output");
print OUT "#unpaired probabilities\n" if ($nrg);
print OUT "#opening energies\n" if ($prob);

for (1..$#data){
	if ($data[$_] !~ /\#/){
		my @value = (split('\t',$data[$_]));
		my @converted;
		for (1..$#value){
			if ($nrg){ ## if $nrg flag is set, we print opening energies, else we print prob of being unpaired
				if ($value[$_] eq 'NA'){
					push @converted,	'NA';
				}
				else{
					push @converted,	nearest(.00000001,($RT*log($value[$_])));
				}
			}
			elsif ($prob){
				if ($value[$_] eq 'NA'){
					push @converted,	'NA';
				}
				else{
					push @converted,	nearest(.00000001,(exp(-$RT/($value[$_]+0000000000.1))));
				}
			}
		}
		print OUT "$value[0]\t";
		print OUT join("\t",@converted),"\n";
	}
	else{
		print OUT $data[$_],"\n";
	}
} 
close (OUT);

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
