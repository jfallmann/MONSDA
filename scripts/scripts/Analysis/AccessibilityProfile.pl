#!/usr/bin/perl
#Script AccessibilityProfile.pl;
#Last changed Time-stamp: <2016-09-12 16:57:33 fall> by joerg

##################
# Libs
##################
use strict;
use warnings;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Pod::Usage;
use Cwd;
use DateTime qw();
use PerlIO::gzip;
use All::Misc qw( unique_array mkdircheck rmdircheck );
use Data::Dumper;
use FindBin::Real qw(Bin); 

##################
# Options 
##################

my ($VERBOSE, $dir, $odir, $pdir, $bed, $ustart, $uend);#, $wobble, $minlength);
$VERBOSE=0;

pod2usage(-verbose => 0)
    unless GetOptions(
	    "dir|d=s"	    => \$dir,     #Dir for Bed
	    "odir|o=s"	    => \$odir,    #Dir for Output
	    "plfold|p=s"    => \$pdir,  #Dir for plfold output
	    "ustart|u=s"    => \$ustart,  #U range from
	    "uend|e=s"      => \$uend,    #U range to	    
	    "bed|b=s"	    => \$bed,     #Bed filename
#	    "wobble|w=s"    => \$wobble,  #Embedding region ### NOT FUNCTIONAL PEAK MERGE FIRST OTHERWISE OVERLAP
#	    "minlength|l=s" => \$minlength,  #Minimum length of peak ### NOT YET NECESSARY
	    "rel|r=s"	    => \my $rel,     #Peak already in relative coordinates
	    "coords|c=s"    => \my $coords,  #Gene coordinates for relative peaks
	    "help|h"	    => sub{pod2usage(-verbose => 1)},
	    "man|m"	    => sub{pod2usage(-verbose => 2)},      
	    "verbose"	    => sub{ $VERBOSE++ }
    );

##################
# Main
##################

my $pid = $$;
(my $job = `cat /proc/$pid/cmdline`)=~ s/\0/ /g;
#print STDERR $job,"\n";

#$wobble = 19 unless ($wobble);
my $appdir = Bin() . "/..";
#print STDERR "Running on dir $appdir, with context $wobble\n";

$dir  =	 cwd() unless ($dir);
chdir ("$dir") || die ("Could not find directory $dir!\n");

my %gencoords;

unless ($rel){
	### Parse Coordinate File
	open(REL,"<:gzip(autopop)",$coords) || die ("Could not open $coords!\n");
	while(<REL>){
		chomp (my $raw = $_);
		my ( $chrom, $start, $end, $gene, $nr, $strand ) = split (/\t/,$raw);
		$strand =~ s/\"//g;
#		print STDERR "$raw\t,$nr, $gene, $chrom, $strand, $start, $end\n";
		$gencoords{$gene} = join(":",$chrom,$start,$end,$strand);
	}
	close(REL);
}

#print Dumper (\%gencoords);

### Parse Bed File
open(BED,"<:gzip(autopop)",$bed) || die ("Could not open $bed in $dir!\n");
my %bed=();

while(<BED>){
        chomp (my $raw = $_);
        push my @line , split (/\t/,$raw);
        push @line, "\." if ( !$line[5] ); 
        (my $chromosome = $line[0])=~ s/chr//g;
        my $start	= $line[1]+1;
        my $end		= $line[2];
        my $name	= uc($line[3]);
        my $score	= $line[4];
        my $strand	= $line[5];
        my $rest	= join("|",@line[6 .. $#line]);
        $rest		= '' if ($rest eq ":");
        
	#        if ($file	 =~ /\.pk/){
	#                $score	 =  $line[5];
	#                $strand =  $line[4];
	#        }
	
#        if (($end-$start+1) < $minlength){
#                $wobble = nearest(1,($minlength-($end-$start+1))/2);
#        }
        
	my $wstart = $start;
	my $wend   = $end;
 
	unless ($rel){
		my ($gchrom,$gstart,$gend,$gstrand) = split(/:/,$gencoords{$name});
#		print STDERR "$gencoords{$name}\n";
		die ("Strand of Peak @line, $strand and Gene $name are not the same, please check!!!\n") if ($gstrand ne $strand);
		if ($gstrand eq "+"){
			$wstart -= $gstart;
			$wend   -= $gstart;
		}
		elsif ($gstrand eq "-"){
			my $tmpstart = $gend - $wstart;
			my $tmpend   = $gend - $wend;	
			$wend	     = $tmpstart;
			$wstart	     = $tmpend;
		}
	}
       
        push @{$bed{$name}}, "$chromosome\:$start\:$end\:$strand\:$score\:$rest\:$wstart\:$wend";
}

#print Dumper (\%bed{"CD551"});

### Parse RNAplfold Output
chdir ("$pdir") || die ("Could not find directory $pdir!\n");
my %folded;
my $RT = (-1.9872041*10**(-3))*(37+273.15);

foreach my $name (keys %bed){
	print STDERR "Processing $name\n";
	open (M, "<:gzip(autopop)","$name\_plfold_med") || die ("Could not open $!");
	open (L, "<:gzip(autopop)","$name\_plfold_long") || die ("Could not open $!");
	my %raw;
	my %us;
	while(<M>){
		chomp($_);
		if ($_ =~ /^Position/){
			my @u = split("\t",$_);
			foreach my $col (1 .. $#u){
				$us{$u[$col]}=$col;
			}
		}
		my @tmp = split('\t',$_);
		die ("Uend not in range of plfold prediction, please choose lower uend or rerun plfold with new params!\n") if ($uend > $#tmp+$ustart); ##Not tmp+1
		my $i = shift @tmp;
		$raw{"Med"}{$i} = join("\t",@tmp[$us{$ustart}-1 .. $us{$uend}-1]);
	}

	%us = ();
	while(<L>){
		chomp($_);
		if ($_ =~ /^Position/){
			my @u = split("\t",$_);
			foreach my $col (1 .. $#u){
				$us{$u[$col]}=$col;
			}
		}
		my @tmp = split('\t',$_);
		die ("Uend not in range of plfold prediction, please choose lower uend or rerun plfold with new params!\n") if ($uend > $#tmp+$ustart);
		my $i = shift @tmp;
		$raw{"Long"}{$i} = join("\t",@tmp[$us{$ustart}-1 .. $us{$uend}-1]);
	}		
	close(M); close(L);
	
	foreach my $peak (@{$bed{$name}}){
		my ($chromosome,$start,$end,$strand,$score,$rest,$wstart,$wend) = split(/\:/,$peak);
#		print STDERR "$peak\t$name\t$wstart\t$wend\n";		
		foreach my $pos ($wstart .. $wend){
			push @{$folded{$name}{"Med"}{$pos-$wstart+1}}  , $raw{"Med"}{$pos};
			push @{$folded{$name}{"Long"}{$pos-$wstart+1}} , $raw{"Long"}{$pos};
		}
	}
}

#print Dumper (\%folded);

### Generate Output
$odir = $dir."/../AccProfile" unless ($odir);

if (!-d" $dir\/$odir"){
	print STDERR "Will create dir $odir!\n";
	mkdircheck("$dir\/$odir");
}

chdir ("$dir\/$odir") or die "Not possible to change to dir $dir\/$odir!\n";

foreach my $name (keys %folded){
#	print STDERR "$name\n";
	open (OM,">",$name."_plfold_med_profile");
	print OM join("\t","Position",($ustart .. $uend))."\n";
	open (OL,">",$name."_plfold_long_profile");
	print OL join("\t","Position",($ustart .. $uend))."\n";
	foreach my $pos (sort {$a<=>$b} keys %{$folded{$name}{"Med"}}){
		foreach my $pk (@{$folded{$name}{"Med"}{$pos}}){
			print OM "$pos\t".$pk."\n";
		}
	}
	foreach my $pos (sort {$a<=>$b} keys %{$folded{$name}{"Long"}}){
		foreach my $pk (@{$folded{$name}{"Long"}{$pos}}){
			print OL "$pos\t".$pk."\n";
		}
	}
	close(OM);
	close(OL);
}
