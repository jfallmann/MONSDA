#!/usr/bin/env perl

### use things ###
use strict;
use warnings;
use autodie;
use PerlIO::gzip;
use Data::Dumper;
use List::Util qw (max);
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Cwd;
use File::Slurp;
use Math::CDF;
### use own modules
use FindBin::Real qw(Bin); # locate this script
use lib Bin() . "/../lib";
use Collection;

my ( $dir, $odir, $peakout, $peakfile, $ratio, $limitratio, $abslim, $dist, $width, $cutoff );
my $VERBOSE=0;
pod2usage(-verbose => 0)
    unless GetOptions(
	    "dir|d=s"			=> \$dir,
	    "odir|o=s"			=> \$odir,
	    "peakout|x=s"       => \$peakout,
	    "peak|p=s"			=> \$peakfile,
	    "ratio|r=s"			=> \$ratio,
	    "limitratio|l=s"	=> \$limitratio,
	    "absolutelimit|a=s" => \$abslim,
	    "dist|t=s"			=> \$dist,
	    "width|w=s"			=> \$width,
	    "cutoff|c=s"		=> \$cutoff,
	    "help|h"			=> sub{pod2usage(-verbose => 1)},
	    "man|m"				=> sub{pod2usage(-verbose => 2)},
	    "verbose"			=> sub{ $VERBOSE++ }
    );

$ratio = 0.5 unless ($ratio);
$limitratio = 0.05 unless ($limitratio);
$dist  = 1 unless ($dist);
$width = 0 unless ($width);
$cutoff = 1 unless ($cutoff);
$dir = cwd() unless ($dir);
$odir =~ s/$odir/\Q$odir\E/g if($odir);
$odir = "$dir"."\/Subpeaks" unless ($odir);

my $pid = $$;
(my $job = `cat /proc/$pid/cmdline`)=~ s/\0/ /g;
print STDERR $job,"\n";

#print STDERR "Processing PrePeak File\n";
my @peakfiles = split(/,/,$peakfile);
my ($pks,$maxheight) = Collection::parse_Prepeaks(\@peakfiles);
### Define Limit and remove pks if necessary
my $limit = $abslim || $maxheight*$limitratio;
#print STDERR "Discarding peaks with height below: $limit\n";

foreach my $pk (keys %{$pks}){
	my $height = (split(/_/,$pks->{$pk}))[1];
#	print STDERR $pks->{$pk},"\t",$height,"\t",$limit,"\n";
	delete $pks->{$pk} if ($height < $limit);
}

#print STDERR "Merging Peaks\n";
my $merged = Collection::merge_peak_regions($pks,$dist);

#write_file 'merged.log', Dumper($merged);
#write_file 'pks.log', Dumper($pks);

#print STDERR "Merging Areas\n";
my $peaks = Collection::merge_peak_areas($pks,$merged);
$merged = {};
$pks = {};

#write_file 'peaks.log', Dumper($peaks);

#print STDERR "Splitting Peaks\n";
my $finalpks = Collection::split_peaks($peaks,$ratio,$limit,$width,$cutoff);
$peaks = {};

#write_file 'finalpks.log', Dumper($finalpks);

#print STDERR "Re-merging peaks in distance $dist\n";
#$dist = 2; ### At this moment peaks are half open!
my $finalmerged = Collection::merge_peak_regions($finalpks,$dist);

#write_file 'finalmerged.log', Dumper($finalmerged);

#print STDERR "Re-merging Areas\n";
my $finalpeaks = Collection::merge_peak_areas($finalpks,$finalmerged);

#write_file 'finalpeaks.log', Dumper($finalpeaks);

#print STDERR "Printing Results\n";

foreach my $newpeak ( keys %{$finalpeaks}){
    my @tmp     = split(/\_/,$newpeak);
    push @tmp,  split(/\_/,$finalpeaks->{$newpeak});
    (my $chrom   = $tmp[0]) =~ s/=/\_/g;
    my $strand  = $tmp[1];
    $strand     = '.' if ($strand eq 'u');
    my $start   = $tmp[2];
    my $end     = $tmp[3];
    my $profile = $tmp[4];
    my $height  = $tmp[5];
    my $summit  = $tmp[6];
    my $area    = $tmp[7];
	my $pval    = $tmp[8];

	my @pro = split(/\|/,$profile);

	if ((split("\:",$pro[-1]))[-1] == 0){
		pop(@pro);
		$profile = join("|",@pro);
		$end -= 1;
	}

    print STDOUT "$chrom\t$start\t$end\t$profile\t$height\t$strand\t$summit\t$area\t$pval\n";
}
