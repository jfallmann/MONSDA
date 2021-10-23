#!/usr/bin/env perl

#Script ExtendBed.pl;
#Last changed Time-stamp: <2020-11-16 13:22:19 fall> by Joerg Fallmann <joerg.fallmann@univie.ac.at>
###############
###Use stuff
###############
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Pod::Usage;
use PerlIO::gzip;
use Cwd;
use File::Path qw(make_path remove_tree);
use File::Basename qw(fileparse);
use Math::Round;
use autodie;

### use own modules
use FindBin::Real qw(Bin); # locate this script
use lib Bin() . "/../lib";
use Collection;

###############
###Variables
###############

my $VERBOSE = 0;
my ( $g, $b, $o, $l, $r, $e, $u, $d, $m, $trim, $keep, $introns );

$l=-inf;
$r=-inf;
$d=-inf;
$u=-inf;
$e=-inf;
$trim=-inf;

###############
###Command Line Options
###############
Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 0) unless GetOptions(
    "genome|g=s"   => \$g,
    "bedfile|b=s"  => \$b,
    "outfile|o=s"  => \$o,
    "left|l=i"     => \$l,
    "right|r=i"    => \$r,
    "extend|e=i"   => \$e,
    "up|u=i"       => \$u,
    "down|d=i"     => \$d,
    "trim|t=i"     => \$trim,
    "min|f=i"      => \$m,
    "original|k=s" => \$keep,
    "bedtwelve|i=s"=> \$introns,
    "help|h"       => sub{pod2usage(-verbose => 1)},
    "man|m"        => sub{pod2usage(-verbose => 2)},
    "verbose"      => sub{ $VERBOSE++ }
    );

###############
###MAIN
###############

my $pid = $$;
(my $job = `cat /proc/$pid/cmdline`)=~ s/\0/ /g;
print STDERR $job,"\n";

unless ( -f $b ) {
    warn "Could not open bed file for processing, please provide via the -b option!\n";
    pod2usage(-verbose => 0);
}

unless ( $e > -inf || $l > -inf || $r > -inf || $u > -inf || $d > -inf){
    warn "No number of flanking nucleotides chosen, output would be input! Please provide them via the -l, -r, -e, -u or -d shortoption\n";
    pod2usage(-verbose => 0);
}

if ($e > -inf){
    $l=nearest(1,$e/2);
    $r=nearest(1,$e/2);
	print STDERR "Extending both sides equally if possible, upstream $l and downstream $r for setting $e\n";
}

print STDERR "Extending left $l, right $r, downstream $d and upstream $u\n";

my $sizes;
if (-f $g){
    open (my $Genome, "<:gzip(autopop)", $g);
    while (<$Genome>){
        chomp $_;
        my ($chr, $size)=split (/\t/, $_);
        $sizes->{$chr}=$size;
    }
}
else{
    $sizes=Collection::fetch_chrom_sizes($g, $g.".chrom.sizes");
}

my @ks = keys %$sizes;
my $tag = $ks[0] =~ /^chr/ ? 1 : 0;

my $filextension;
if ($d > -inf){
    $filextension .= "_".$d."_fromEnd";
}
if ($u > -inf){
    $filextension .= "_".$u."_fromStart";
}
if ($l > -inf){
    $filextension .= "_".$l."_upstream";
}
if ($r > -inf){
    $filextension .= "_".$r."_downstream";
}
if ($e > -inf){
    $filextension = "_".$e."_equal";
}
if ($trim > -inf){
    $filextension = "_".$e."_trimmedto";
}

my($filename, $dirs, $suffix) = fileparse($b, '.bed');
if ($b =~ '.bed.gz$'){
    ($filename, $dirs, $suffix) = fileparse($b, '.bed.gz');
}

$o = $filename.$filextension.$suffix unless ($o);

if (-f $o){
    warn "Outputfile $o does exist, will be overwritten!\n";
}

open (my $Out, ">:gzip",$o) or die "$!";

open (my $Bed, "<:gzip(autopop)", $b) or die "$!";
my @featurelist; ## This will become a FeatureChain
while(<$Bed>){
    chomp $_;
    my ($chrom, $start, $end, $id, $score, $strand, @rest) = split (/\t/, $_);    
    $chrom =~ s/^chr//g if !$tag;
    my $right  = 0;
    my $left   = 0;
    my $width  = nearest(1, ($end-$start)/2);
    if ($keep){
        my $original = join("\t", $start, $end);
        push @rest, $original;
    }
    $width = 0 if ($d  > -inf || $u  > -inf || $e > -inf || $l  > -inf || $r  > -inf);
    $end+=1 if (($end-$start)%2) && $width > 0;

    if (defined $m && (($end - $start) >= $m) ){
        if (@rest){
            print $Out "$chrom\t$start\t$end\t$id\t$score\t$strand\t", join("\t", @rest), "\n";
        }
        else{
            print $Out "$chrom\t$start\t$end\t$id\t$score\t$strand\n";
        }
        next;
    }
    else{
        my $tstart = $start;
        my $tend = $end;
        
        if ($strand eq "+" || $strand eq '.' || $strand eq 'u'){
            if ($d >= 0){
                $start = $tend;
                if ($d == 0){
                    $start-=1;
                }
                $r = $d;
                $l = 0;
            }
            if ($u >= 0){
                $end = $tstart;
                if ($u == 0){
                    $end += 1;
                }
                $end += 1 if $end == 0;
                $l = $u;
                $r = 0;
            }
            $right=$r;
            $left=$l;
        }
        elsif ($strand eq "-"){
            if ($u >= 0){
                $start = $tend;
                if ($u == 0){
                    $start-=1;
                }
                $l = $u;
                $r = 0;
            }
            if ($d >= 0){
                $end = $tstart;
                if ($d == 0){
                    $end += 1;
                }
                $end += 1 if $end == 0;
                $r = $d;
                $l = 0;
            }
            $right=$l;
            $left=$r;
        }
        
        if (($right-$width) <= 0){
            $right = 0;
        }
        else{
            $right-=$width;
        }
        if (($left-$width) <= 0 ){
            $left = 0;
        }
        else{
            $left-=$width;
        }
        if ($end+$right >= $sizes->{$chrom}){
            $end = $sizes->{$chrom};
        }
        else{
            $end += $right;
        }
        if ( $start-$left >= 0 ){
            $start -= $left;
        }
        else{
            $start = 0;
        }
        if ($trim > -inf){
            if (($end-$start) > ($e+1)){
                my $finalsize = nearest(1,(($end-$start)-$e)/2);
                $start+=$finalsize;
                $end-=$finalsize;
            }
        }
        if ($start >= $end){
            if (($d >= 0 && ($strand eq "+" || $strand eq '.' || $strand eq 'u')) || ($u >= 0 && $strand eq "-")){
                $end = $start + 1
            }
            elsif (($d >= 0 && $strand eq "-") || ($u >= 0 && ($strand eq "+" || $strand eq '.' || $strand eq 'u'))){
                $start = $end - 1;
            }
        }
        print $Out "$chrom\t$start\t$end\t$id\t$score\t$strand";

        if (@rest){
            if ($introns){
                my @blocks = split(',',$rest[5]);
                my @blocksizes = split(',',$rest[4]);

                for (1..$#blocks){
                    $blocks[$_] += $tstart - $start;
                }
                $blocksizes[0] += $tstart - $start;
                $blocksizes[-1] += $end - $tend;

                $rest[4] = join(',',@blocksizes);
                $rest[5] = join(',',@blocks);
            }
            print $Out "\t";
            print $Out join ( "\t", @rest );
        }
        print $Out "\n";
    }
}


###############
###POD
###############

__END__

    =head1 NAME

    ExtendBed.pl - Extends bed entries strand specific one- or two-sided or equal.

    =head1 SYNOPSIS
    ExtendBed.pl [-g I<FILE>] [-b I<FILE>] [-o I<FILE>] [-e I<Interger>] [-l I<Interger>] [-r I<Interger>]
    [options]

    =head1 OPTIONS

    =over

    =item B<-g>

    File containing chromosome sizes

    =item B<-b>

    Bed file for extension

    =item B<-o>

    Output file name

    =item B<-e>

    Extension to total length from both sides

    =item B<-l>

    Extension to total length from 5' end

    =item B<-r>

    Extension to total length from 3' end

    =item B<-f>

    Extend if length smaller than -f (e.g. for MEME)

    =item B<--help -h>

    Print short help

    =item B<--man>

    Prints the manual page and exits

    =back

    =head1 DESCRIPTION

    This program extends Bed files to a total length of at least -e, -r or -l nucleotides, or at least to begin or end of the corresponding chromosome

    =head1 AUTHOR

    Joerg Fallmann E<lt>joerg.fallmann@univie.ac.atE<gt>

    =cut


    ##################################END################################
