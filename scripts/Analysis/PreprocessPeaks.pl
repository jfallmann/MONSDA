#!/usr/bin/env perl

#### use things ###
use strict;
use warnings;
use autodie;
use PerlIO::gzip;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Cwd;
### use own modules
use FindBin::Real qw(Bin); # locate this script
use lib Bin() . "/../lib";
use Collection;

my ( $dir, $odir, $peakfile, $filter);
my $VERBOSE=0;
pod2usage(-verbose => 0)
    unless GetOptions(
		"dir|d=s"	 => \$dir,
		"odir|o=s"	 => \$odir,
		"peak|p=s"	 => \$peakfile,
		"filter|f=s" => \$filter,
		"help|h"	 => sub{pod2usage(-verbose => 1)},
		"man|m"		 => sub{pod2usage(-verbose => 2)},
		"verbose"	 => sub{ $VERBOSE++ }
    );

#my $pwd = cwd();
$dir = cwd() unless ($dir);
$odir =~ s/$odir/\Q$odir\E/g if($odir);
$odir = "$dir"."\/Subpeaks" unless ($odir);
$filter = qr(\Q$filter\E) if($filter);

my $pid = $$;
(my $job = `cat /proc/$pid/cmdline`)=~ s/\0/ /g;
print STDERR $job,"\n";

chdir ("$dir");

#Read beds;
push my @beds,split(',',$peakfile);

print STDERR "Processing Bed File and annotating profile\n";

#parse bedfile
my ($unique,$chl) = Collection::parse_bedgraph(\@beds);#, $unique, $chl);

print STDERR "Printing bed with profile\n";

foreach my $pk (keys %{$unique}){
	my @tempuni = split(/\_/,$pk);
	push @tempuni , split(/\_/,$unique->{$pk});
	(my $chromosome = $tempuni[0])=~ s/=/\_/g;
	$chromosome =~ s/:/\s/g;
	my $start      = $tempuni[1];
	my $end        = $tempuni[2];
	my $strand     = $tempuni[3];
	my $name       = $tempuni[4];
	my $score      = $tempuni[5];
	my $summit     = $tempuni[6];
	my $rest       = $tempuni[7];

	if ($filter){
		next if $name !~ /$filter/;
	}

	$strand = '.' if $strand eq 'u';

	my $profile;
	for ($start..$end-1){
		$profile->{$_}=$score;
	}
	my $area;
	my @tmp;
    for my $loci (sort{$a <=> $b} keys %{$profile}){
#	    push @tmp, join(':',$loci,$profile->{$loci});
	    $summit = $profile->{$loci} if ( $summit < $profile->{$loci} );
		$area+=$profile->{$loci};
    }
    if ($rest eq 'undef'){
		$rest = $area;
	}
	else{
		$rest=join('\t',$area,$rest);
	}
	#	my $peakprofile = join("|",@tmp);
	my $peakprofile = ($end-$start).':'.$score; #Changing to more sparse peak profile
    print STDOUT "$chromosome\t$start\t$end\t$peakprofile\t$summit\t$strand\t$rest\n";
}
