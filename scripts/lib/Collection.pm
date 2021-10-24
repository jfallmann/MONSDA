package Collection;

#Last changed Time-stamp: <2019-10-18 16:12:54 fall> by joerg
use strict;
use Exporter qw(import);
use Tie::Hash::Indexed; ### Keeping the order
use List::Util qw'max sum shuffle';
use Data::Dumper;
use Set::IntervalTree;
use IPC::Cmd qw(can_run run);
use Path::Class qw(dir file);
use File::Path qw(make_path remove_tree);
use Math::Round;
use Carp;
use File::Slurp;
use FindBin::Real qw(Bin); # locate this script
use Math::CDF;

#use PadWalker qw/peek_my/; ### Check variables used in script

### use own modules
#use lib Bin() . "/../lib";
### use own modules relative to current location
#$libpath = Bin() . "/../lib";$ENV{PWD} unless ($libpath);
#### use Bio things, check version of git repos to fit mysql if using ensembl-mysql for annotation ###
#	use lib $libpath . "/ensembl/modules";
#	use lib $libpath . "/ensembl-compara/modules";
#	use lib $libpath . "/bioperl-live";
#	use Bio::Seq;
#	use Bio::SeqIO;
#	use Bio::AlignIO;
#	use Bio::LocatableSeq;
#	use Bio::EnsEMBL::Registry;
#	use Bio::SeqFeature::Generic;

our @ISA = qw(Exporter);
our @EXPORT_OK = qw( uniquearray mkdircheck rmdircheck fetch_chrom_sizes);
our %EXPORT_TAGS = ();
our $VERSION = 1.00;

sub parse_libpath{
	my $current = shift;
	my $lib = shift || '';
	#return
	if ($lib && $lib =~ /^\//){
		return($lib);
	}
	else{
		return($current."/".$lib);
	}
}

sub fetch_chrom_sizes{
  my $species = shift;
  my $file    = shift;
  my %sizes;
  my @chromsize;
  my $chrtag = 0;
  my $this_function = (caller(0))[3];

  my $genome  = 'hg38';
  $genome     = 'mm10' if ($species eq "Mus_musculus");
  $genome     = 'danRer10' if ($species eq "Danio_rerio");
  $genome     = 'dm6' if ($species eq "Drosophila_melanogaster");
  $genome     = 'ce245' if ($species eq "Caenorhabditis_elegans");

  if(!-e $file){
	  my $test_fetchChromSizes = can_run('fetchChromSizes') or
		  croak "ERROR [$this_function] fetchChromSizes utility not found";

	  my $cmd = "fetchChromSizes $genome";
	  my( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
		  run(command => $cmd, verbose => 0);
	  if ($success){
		  @chromsize = @{$stdout_buf};
		  print STDERR "Printing chromsizes to $file!\n";
		  open(my $fh,">",$file)||die "Could not open $file\n";
		  foreach (@chromsize){
			  next if ($_ =~ /^#/);
			  $chrtag = 1 if ($_ =~ /^chr/);
			  (my $entry = $_) =~ s/^chr//g;
			  print $fh $entry."\n";
		  }
		  close($fh);
	  }
	  else{
		  carp "WARN [$this_function] Using UCSCs fetchChromSizes failed, trying alternative mysql fetch!\n";
		  my $test_fetchChromSizes = can_run('mysql') or
			  croak "ERROR [$this_function] mysql utility not found";
		  $cmd = "mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \"select chrom, size from $genome.chromInfo\"";  ### Alternative to UCSC fetchChromSizes, has mysql dependency
		  my( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
			  run(command => $cmd, verbose => 0);
		  if ($success){
			  @chromsize = @{$stdout_buf};
			  print STDERR "Printing chromsizes to $file!\n";
			  open(my $fh,">",$file)||die "Could not open $file\n";
			  foreach (@chromsize){
				  next if ($_ =~ /^#/);
				  $chrtag = 1 if ($_ =~ /^chr/);
				  (my $entry = $_) =~ s/^chr//g;
				  print $entry $_."\n";
			  }
			  close($fh);
		  }
		  else{
			  carp "ERROR [$this_function] External command call unsuccessful\n";
			  carp "ERROR [$this_function] this is what the command printed:\n";
			  print join "", @$full_buf;
			  confess "Fetching of chromosome sizes failed, please either download fetchChromSizes from the UCSC script collection, or install mysql!\n";
		  }
	  }
  }
  else{
	  print STDERR "Reading chromsizes from $file!\n";
	  open(my $fh,"<:gzip(autopop)",$file)||die "Could not open $file\n";
	  push @chromsize, <$fh>;
	  close($fh);
  }
  foreach (@chromsize){
    chomp($_);
    foreach (split(/\n/,$_)){
		my ($chr,$size)=split (/\t/,$_);
		$chrtag = 1 if ($chr =~ /^chr/);
		$chr =~ s/^chr//g;
		$sizes{$chr}=$size;
    }
  }
  return(\%sizes, $chrtag);
}

sub fetch_chrom_sizes_old{
  my $species = shift;
  my $file    = shift;
  my %sizes;
  my @chromsize;
  my $this_function = (caller(0))[3];

  unless(-f $file){
	  my $test_fetchChromSizes = can_run('fetchChromSizes') or
		  croak "ERROR [$this_function] fetchChromSizes utility not found";

	  my $cmd = "fetchChromSizes $species";
	  my( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
		  run(command => $cmd, verbose => 0);
	  if ($success){
		  @chromsize = @{$stdout_buf};
	  }
	  else{
		  carp "WARN [$this_function] Using UCSCs fetchChromSizes failed, trying alternative mysql fetch!\n";
		  my $test_fetchChromSizes = can_run('mysql') or
			  croak "ERROR [$this_function] mysql utility not found";
		  $cmd = "mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \"select chrom, size from $species.chromInfo\"";  ### Alternative to UCSC fetchChromSizes, has mysql dependency
		  my( $success, $error_message, $full_buf, $stdout_buf, $stderr_buf ) =
			  run(command => $cmd, verbose => 0);
		  if ($success){
			  @chromsize = @{$stdout_buf};
		  }
		  else{
			  carp "ERROR [$this_function] External command call unsuccessful\n";
			  carp "ERROR [$this_function] this is what the command printed:\n";
			  print join "", @$full_buf;
			  confess "Fetching of chromosome sizes failed, please either download fetchChromSizes from the UCSC script collection, or install mysql!\n";
		  }
	  }
  }
  else{
	  open(my $fh,"<",$file)||die "Could not open $file\n";
	  push @chromsize, <$fh>;
	  close($fh);
  }
  foreach (@chromsize){
    chomp($_);
    foreach (split(/\n/,$_)){
      my ($chr,$size)=split (/\t/,$_);
      $sizes{$chr}=$size;
    }
  }
  return(\%sizes);
}

sub dump_vars_global{
	#### Dump all vars in main

	print "Dumping global vars\n";
	{
		no strict 'refs';

		foreach my $entry ( keys %main:: )
		{
			print "$entry\n";
		}
	}
}

sub dump_vars_local{
	### Dump all vars in

	print "Dumping lexical vars\n";
	my $pad = peek_my(0);
	print Dumper ($pad);
}

sub merge_peak_regions{
	### IN
	my $pks	 = shift;
	my $dist = shift;
	### OUT
	my $locations = {};
	my $final	  = {};

	foreach my $pk ( keys %{$pks} ){

		my @tmp     = split(/\_/,$pk);
		push @tmp, split(/\_/,$pks->{$pk});
		my $chrom   = $tmp[0];
		my $strand  = $tmp[1];
		my $start   = $tmp[2];
		my $end     = $tmp[3];
		#	my $which = "plus" if ($strand == +);
		#	my $which = "minus" if ($strand eq "-");
		push @{$locations->{"$chrom\_$strand"}} , "$start,$end";
	}

	for my $key ( keys %{$locations} ) {
		#	print "$chrom: ";
		my @tmp = split (/\_/,$key);
		my $chrom  = $tmp[0];
		my $strand = $tmp[1];
		### Get peaks for region
		my $regions = $locations->{"$key"};
		my @sorted = sort { $a <=> $b } @{$regions};
		my @overlap;
		my ($ts ,$te);
		for (0..$#sorted){
			my ($start,$end) = split (/,/,$sorted[$_]);
			if ($_ == 0){ #BEGIN
				$ts = $start;
				$te = $end;
				if ($#sorted == 0){
					$final->{"$key\_$ts\_$te"}="$ts,$te";
				}
			}
			else{ #NEXT
				if (($te + $dist) >= $start ){# If Overlap, BED_SPK FORMAT, not yet half open, so >= and not >
					my ($checklo, $checkhi);
					foreach my $entry (@overlap){
						$checklo = 1 if ($entry == $_-1);
						$checkhi = 1 if ($entry == $_);
					}
					push @overlap, $_-1 unless ($checklo == 1);
					push @overlap, $_ unless ($checkhi == 1);
					$te = $end if ($end > $te);
					$ts = $start if ($start < $ts);
				}
				else{# No overlap, safe if there was one before
					if (defined $overlap[0] && $overlap[0] ne ""){
						#print STDERR "No Overlap, but one before, so safe: @overlap ";
						#			my @unique_overlap = @{unique_array(@overlap)};
						#			@overlap = @unique_overlap;
						my $ol;
						while (@overlap){
							my $i = shift @overlap;
							$ol .= "$sorted[$i];";
						}
						$ol = substr ($ol,0,-1);

						$final->{"$key\_$ts\_$te"}="$ol" if (defined $ol && $ol ne "");
						$ts = $start;
						$te = $end;
					}
					else{
						$final->{"$key\_$ts\_$te"}="$ts,$te"; #if ($_ == 1);
						$ts = $start;
						$te = $end;
					}
				}

				if ($_ == $#sorted){
					if (defined $overlap[0] && $overlap[0] ne ""){
						#			my @unique_overlap = @{unique_array(@overlap)};
						#			@overlap = @unique_overlap;
						my $ol;
						while (@overlap){
							my $i = shift @overlap;
							$ol .= "$sorted[$i];";
						}
						$ol = substr ($ol,0,-1);

						$final->{"$key\_$ts\_$te"}="$ol" if (defined $ol && $ol ne "");
					}
					else{
						$final->{"$key\_$ts\_$te"}="$start,$end";
					}
				}
			}
		}
	}
	return ($final);
}

sub merge_peak_areas{

	### IN
	my $peaks = shift;
	my $merged = shift;
	### OUT
	my $final = {};

	for my $key ( keys %{$merged} ) {
		my @tmp = split (/\_/,$key);
		my $chrom  = $tmp[0];
		my $strand = $tmp[1];
		my $st = $tmp[2];
		my $en = $tmp[3];
		next if (defined $final->{"$chrom\_$strand\_$st\_$en"}); ## Entry exists, go on
		my @oldpeaks;
		push @oldpeaks , split(/;/,$merged->{"$key"});

		my $poscov;
		my $height = 0;
		my $sum    = 0;
		my $mean   = 'NA';
		foreach my $op (@oldpeaks){
			my ($start,$end) = split(/,/,$op);
			my $oldkey  = "$chrom\_$strand\_$start\_$end";
			my $pkold   = $peaks->{"$oldkey"};
			my @tmp     = split(/\_/,$oldkey);
			push @tmp,  split(/\_/,$pkold);
			my $pst     = $tmp[2];
			my $pe      = $tmp[3];
			my @pro = @{parse_peakprofile([split(/\|/,$tmp[4])], $pst)};
			my $nu      = $pst-1;
			$mean       += $tmp[8];

			foreach my $pos (@pro){
				my @bla = split (/\:/,$pos);
				$bla[0] = ($bla[0]-$nu) if ($bla[0] >= $pst);
				if ($bla[0] == 1){
					$nu++;
					$poscov->{"$nu"} = $bla[1] unless (defined $poscov->{"$nu"});
					if ($bla[1] > $height){
						$height = $bla[1];
						$sum    = $nu;
					}
				}
				else{
					for (my $i = 1;$i <= $bla[0];$i++){
						$nu++;
						$poscov->{"$nu"} = $bla[1] unless (defined $poscov->{"$nu"});
						if ($bla[1] > $height){
							$height = $bla[1];
							$sum    = $nu;
						}
					}
				}
			}
		}
		$mean = nearest(1,$mean/($#oldpeaks+1));

		my ($nprofile, $area) = generate_peakprofile($st, $en, $poscov);
		my $pro = substr($nprofile,0,-1);
		my $pval=Math::CDF::ppois($sum,$mean);

		$final->{"$chrom\_$strand\_$st\_$en"}="$pro\_$height\_$sum\_$area\_$pval" unless (defined $final->{"$chrom\_$strand\_$st\_$en"});
	}
	return ($final);
}

sub parse_peakprofile{
	###IN
	my $tpro  = shift;
	my $start = shift;
	###RET
	my @pro;
	###DO
	my $tstart = $start;
	foreach my $t (@{$tpro}){
		my ($add,$score) = split('\:',$t);
		my $tend = $tstart + $add - 1;
		for ($tstart..$tend){
			push @pro, $_.':'.$score;
		}
		$tstart = $tend + 1;
	}
	return (\@pro);
}

sub generate_peakprofile{
	my $st	   = shift;
	my $en	   = shift;
	my $poscov = shift;
	my $idx    = 1;

	my $nprofile;
	my $area;

	for (my $a = $st;$a < $en;$a++){
		my $h;
		$h = $poscov->{"$a"} if defined $poscov->{"$a"};
		$h = 0 unless ($h);
		$area += $h;
		if (defined $poscov->{$a+1} && $h == $poscov->{$a+1}){
			$idx++;
		}
		else{
			$nprofile .= $idx.':'.$h.'|';
			$idx = 1;
		}
	}
	return ($nprofile, $area);
}

sub split_peaks{
	### IN

	my $peaks  = shift;
	my $ratio  = shift;
	my $limit  = shift;
	my $width  = shift;
	my $cutoff = shift;

	### OUT
	my $newpos;
	my $finalpks;

	foreach my $pks (keys %{$peaks}){
		# $final{"$chrom\_$strand\_$st\_$en"}="$pro\_$height\_$sum\_$area\_$pval" unless (defined $final{"$chrom\_$strand\_$st\_$en"});

		my @tmp     = split(/\_/,$pks);
		push @tmp,  split(/\_/,$peaks->{$pks});
		my $chrom   = $tmp[0];
		my $strand  = $tmp[1];
		my $start   = $tmp[2];
		my $end     = $tmp[3];
		my @profile = @{parse_peakprofile([split(/\|/,$tmp[4])],$start)};
		my $height  = $tmp[5];
		my $summit  = $tmp[6];
		my $area    = $tmp[7];

		next if ($height < $limit);

		my @heights;
		my @nuks;
		my $newpeaks;
		my $poscov;
		my $mean;my $meanlen=0;

#		tie my %poscov, 'Tie::Hash::Indexed'; ### NO reason for that anymore with new split algo

		foreach my $pos (@profile){ ### We need that to find the position of the summit
			my @bla = split (/\:/,$pos);
			$poscov -> {"$bla[0]"} = $bla[1];
			$mean+=$bla[1];
			$meanlen++;
		}

		$mean=nearest(1,$mean/$meanlen); # Get mean for poisson

		### Divide peaks in smaller regions and check

		my ( $nstart, $nend, $nheight ) = (0,0,0);
		my @subpeaks;
		my $round = "firstround";
		my $st = $start;
		my $en = $end;

		if ($round eq "firstround"){
			my $startfrom = $summit;
			for (my $i=1;($startfrom - $i) >= $start;$i++){
				my $pos = $startfrom - $i;
				my $nheight = $poscov -> {"$pos"} || undef;
				if ($nheight < ($height * $ratio) or $nheight == 0){
					$nstart = $pos + 1;
				}
				last if ($nstart != 0);
			}

			for (my $i = 1; ($startfrom + $i) <= $end; $i++){
				my $pos = $startfrom + $i;
				my $nheight = $poscov -> {"$pos"} || undef;
				if ($nheight < ($height * $ratio) or $nheight == 0){
					$nend = $pos - 1;
				}
				last if ($nend != 0);
			}

			$nstart = $start if ($nstart  == 0);
			if ($nend  == 0){
				$nend = $end;
			}
			else{
				$nend = $nend + 1;# make bedformat conform
			}

#			print STDERR join("\t",$nstart,$nend,($nend-$nstart))."\n";

			next if (($nend-$nstart) < $width);

#			print STDERR ($nend-$nstart)."\n";

			my ($nprofile, $splitarea) = generate_peakprofile($nstart, $nend + 1, $poscov);
			my $pro = substr($nprofile,0,-1);

			push @{$newpos} , "$nstart,$nend,$height,$chrom,$strand,$splitarea,$pro,$summit,$mean";

			my ( @loar, @upar);
			tie my %locov, 'Tie::Hash::Indexed';
			tie my %hicov, 'Tie::Hash::Indexed';

			for (my $i=$start;$i<=$nstart-1;$i++){
				push @loar, $i;
			}

			for (my $i=$nend+1;$i<=$end;$i++){
				push @upar, $i;
			}

			%locov = map { $_ => $poscov -> {"$_"} } @loar unless ( $nstart == 0 );
			%hicov = map { $_ => $poscov -> {"$_"} } @upar unless ( $nend == 0 );

			my $loref = \%locov;
			my $hiref = \%hicov;

			push @subpeaks, $loref unless ( $nstart == 0 );
			push @subpeaks, $hiref unless ( $nend == 0 );

			$round = "secondround";
		}

		if ($round eq "secondround"){

			while (@subpeaks){
				my $hash;
				$hash = shift @subpeaks;
				#		my $one_element_warning;
				#		$hash = shift @subpeaks unless ($#subpeaks == 0);
				#		$hash = $subpeaks[0] if ($#subpeaks == 0);
				#		$one_element_warning = 1 if ($#subpeaks == 0);
				tie my %subpeak, 'Tie::Hash::Indexed';
				%subpeak = %{$hash};
				my $max = max(values %subpeak);
				my $check = ($height*$cutoff);
				next if ( $max < $check || $max < $limit);
				my @sf = grep { $subpeak{"$_"} == $max } keys %subpeak;
				my $nmax=$max*$ratio;
				my $startfrom = $sf[0];

				my ($nstart,$nend)=(0,0);

				$en = ( (keys %subpeak)[-1]);
				$st = ( (keys %subpeak)[0]);

				for (my $i=1;($startfrom-$i) >= $st;$i++){
					my $pos = $startfrom-$i;
					my $nheight = $poscov  -> {"$pos"};
					if ($nheight < $nmax  or $nheight == 0){
						$nstart = $startfrom-$i+1;
					}
					last if ($nstart != 0);
				}

				for (my $i = 1; ($startfrom+$i) <= $en; $i++){
					my $pos = $startfrom+$i;
					my $nheight = $poscov -> {"$pos"};
					if ($nheight < $nmax or $nheight == 0){
						$nend = $startfrom+$i-1;
					}
					last if ($nend != 0);
				}

				$nstart = $st if ($nstart == 0);
				if ($nend  == 0){
					$nend = $end;
				}
				else{
					$nend = $nend + 1;# make bedformat conform
				}

				next if (($nend-$nstart) < $width);

				my ($nprofile, $splitarea) = generate_peakprofile($nstart, $nend+1, $poscov);
				my $pro = substr($nprofile,0,-1);
				push @{$newpos} , "$nstart,$nend,$max,$chrom,$strand,$splitarea,$pro,$startfrom,$mean" unless (($nend-$nstart)<$width);
				my @loar = ($st..$nstart-1);
				my @upar = ($nend+1..$en);

				tie my %locov, 'Tie::Hash::Indexed';
				tie my %hicov, 'Tie::Hash::Indexed';

				%locov = map { $_ => $poscov -> {"$_"} } @loar;
				%hicov = map { $_ => $poscov -> {"$_"} } @upar;

				my $loref = \%locov;
				my $hiref = \%hicov;

				push @subpeaks, $loref unless ( $nstart == 0 );
				push @subpeaks, $hiref unless ( $nend == 0 );
				#		push @subpeaks, $loref unless ( $nstart == 0 || ($nend-$nstart) < $width );
				#		push @subpeaks, $hiref unless ( $nend == 0 || ($nend-$nstart) < $width );
				#		$hash = shift @subpeaks if ($one_element_warning == 1);
			}
		}
	}

	while (my $newpeak = shift @{$newpos}){
		my @pos     = split(/,/,$newpeak);
		my $start   = $pos[0];
		my $end     = $pos[1];
		my $height  = $pos[2];
		my $chrom   = $pos[3];
		my $strand  = $pos[4];
		my $area    = $pos[5];
		my $pro     = $pos[6];
		my $summit  = $pos[7];
		my $mean    = $pos[8];

		$finalpks->{"$chrom\_$strand\_$start\_$end"}="$pro\_$height\_$summit\_$area\_$mean" unless (defined $finalpks->{"$chrom\_$strand\_$start\_$end"});
	}

	return ($finalpks);
}

sub unique_array{

	my $arrayref = shift;
	my @array = @{$arrayref};

	my %unique = ();
	foreach my $item (@array)
	{
		$unique{$item} ++;
	}
	my @arrayuid = sort {$a cmp $b} keys %unique;

	return(\@arrayuid);
}

sub mkdircheck {
	my @dirstocreate=();
	my $this_function = (caller(0))[3];
	while(@_){
		push @dirstocreate, shift(@_);
	}
	foreach (@dirstocreate){
		my @total = split(/[\/\\]/,$_);
		my $dir;
		while(@total){
			$dir = dir(shift(@total)) unless (defined $dir);
			$dir = dir($dir,shift(@total));
		}
		my $check = $dir->stringify();
		return if (-d $check);
		make_path($dir,{ verbose => 0 }) or
		    croak "ERROR [$this_function] Cannot create directory $check";
	}
}

sub rmdircheck {
	my @dirstorm=();
	my $this_function = (caller(0))[3];
	while(@_){
		push @dirstorm, shift(@_);
	}
	foreach (@dirstorm){
		my @total = split(/[\/\\]/,$_);
		my $dir;
		while(@total){
			$dir = dir(shift(@total)) unless (defined $dir);
			$dir = dir($dir,shift(@total));
		}
		my $check = $dir->stringify();
		return if (!-d $check);
		remove_tree($dir,{ verbose => 0 }) or
		    croak "ERROR [$this_function] Cannot remove directory $check";
	}
}

sub find_motif {

	my $seq     = uc( $_[0] );
	my $motifref     = $_[1];
	my @motifstoscan = @{$motifref};
	my $length	     = length($seq);
	# return variables
	my %anno;

	foreach my $motif (@motifstoscan){
		my @found;
		my @exact;

		while($seq =~ /(?=($motif))/g){
			my $where = pos($seq);
			push @found, $where;
			push @exact, $1;
			#	    print STDERR "$motif\t$where\t$1\n";
		}

		for (0 .. $#found){
			my $seq_pos	= $found[$_];
			my $ex	= $exact[$_];
			my $mlength = length($ex);
			my $relstrt = ($seq_pos+1)/$length;
			push @{$anno{$motif}{starts}}, $seq_pos + 1;
			push @{$anno{$motif}{ends}},   $seq_pos + $mlength;
			push @{$anno{$motif}{exact}},  $ex;
			push @{$anno{$motif}{relpos}}, $relstrt;
		}
	}
	return ( \%anno );
}

sub connect_db{
	my $species   = shift||undef;
	my $dbversion = shift||undef;
	my $dbtouse   = shift||undef;
	my $mysql     = shift||undef;
	my $spec;

	if ($species eq "Human" or $species eq "human"){$spec			   = "Homo_sapiens";}
	elsif ($species eq "Mouse" or $species eq "mouse"){$spec		   = "Mus_musculus";}
	elsif ($species eq "Zebrafish" or $species eq "zebrafish"){$spec   = "Danio_rerio";}
	elsif ($species eq "Drosophila" or $species eq "drosophila"){$spec = "Drosophila_melanogaster";}
	elsif ($species eq "Worm" or $species eq "worm"){$spec			   = "Caenorhabditis_elegans";}
	else{$spec														   = $species;}

	my ($host,$user,$pass);
	#####Connect with ENSEMBL mysql#####
	($host,$user,$pass) = split(",",$mysql) if (defined $mysql);
	$host = "ensembldb.ensembl.org" unless (defined $host);
	$user = "anonymous" unless (defined $user);
	$pass = undef unless (defined $pass);

	my $registry = 'Bio::EnsEMBL::Registry';
	$registry->load_registry_from_db(
		-host	    => $host,
		-user	    => $user,
		-pass	    => $pass,
		-species    => lc($spec),
		-db_version => $dbversion,
		-dbname	    => $dbtouse,
		-verbose    => 1
	    );

	$registry->set_reconnect_when_lost();
	return $registry;
}

sub connect_compara_db{
	my $species   = shift;
	my $dbversion = shift;
	my $dbtouse   = shift||'';
	my $mysql     = shift||'';
	my $spec;

	if ($species eq "Human" or $species eq "human"){$spec			   = "Homo_sapiens";}
	elsif ($species eq "Mouse" or $species eq "mouse"){$spec		   = "Mus_musculus";}
	elsif ($species eq "Zebrafish" or $species eq "zebrafish"){$spec   = "Danio_rerio";}
	elsif ($species eq "Drosophila" or $species eq "drosophila"){$spec = "Drosophila_melanogaster";}
	elsif ($species eq "Worm" or $species eq "worm"){$spec			   = "Caenorhabditis_elegans";}
	else{$spec														   = $species;}

	my ($host,$user,$pass);
	#####Connect with ENSEMBL mysql#####
	($host,$user,$pass) = split(",",$mysql) if (defined $mysql);
	$host = "ensembldb.ensembl.org" unless (defined $host);
	$user = "anonymous" unless (defined $user);
	$pass = "" unless (defined $pass);

	$dbtouse = 'ensembl_compara_'.$dbversion unless (defined $dbtouse);

	my $compregistry = 'Bio::EnsEMBL::Registry';
	$compregistry->load_registry_from_db(
		-host	    => $host,
		-user	    => $user,
		-pass	    => $pass,
		-species    => 'compara',
		-db_version => $dbversion,
		-dbname	    => $dbtouse,
		-verbose    => 1
	    );

	$compregistry->set_reconnect_when_lost();
	return $compregistry;
}

sub connect_funcgen_db{
	my $species   = shift;
	my $dbversion = shift;
	my $dbtouse   = shift||'';
	my $group     = shift||'';
	my $mysql     = shift||'';
	my $spec;

	if ($species eq "Human" or $species eq "human"){$spec			   = "Homo_sapiens";}
	elsif ($species eq "Mouse" or $species eq "mouse"){$spec		   = "Mus_musculus";}
	elsif ($species eq "Zebrafish" or $species eq "zebrafish"){$spec   = "Danio_rerio";}
	elsif ($species eq "Drosophila" or $species eq "drosophila"){$spec = "Drosophila_melanogaster";}
	elsif ($species eq "Worm" or $species eq "worm"){$spec			   = "Caenorhabditis_elegans";}
	else{$spec														   = $species;}

	my ($host,$user,$pass);
	#####Connect with ENSEMBL mysql#####
	($host,$user,$pass) = split(",",$mysql) if (defined $mysql);
	$host = "ensembldb.ensembl.org" unless (defined $host);
	$user = "anonymous" unless (defined $user);
	$pass = "" unless (defined $pass);

	$dbtouse = lc($spec)."_funcgen_".$dbversion unless (defined $dbtouse);

	my $funcgenregistry = 'Bio::EnsEMBL::Registry';
	$funcgenregistry->load_registry_from_db(
		-host	    => $host,
		-user	    => $user,
		-pass	    => $pass,
		-species    => lc($spec),
		-db_version => $dbversion,
		-dbname	    => $dbtouse,
		-verbose    => 1
	    );

	$funcgenregistry->set_reconnect_when_lost();
	return $funcgenregistry;
}

sub connect_variation_db{
	my $species   = shift;
	my $dbversion = shift;
	my $dbtouse   = shift||'';
	my $group     = shift||'';
	my $mysql     = shift||'';
	my $spec;

	if ($species eq "Human" or $species eq "human"){$spec			   = "Homo_sapiens";}
	elsif ($species eq "Mouse" or $species eq "mouse"){$spec		   = "Mus_musculus";}
	elsif ($species eq "Zebrafish" or $species eq "zebrafish"){$spec   = "Danio_rerio";}
	elsif ($species eq "Drosophila" or $species eq "drosophila"){$spec = "Drosophila_melanogaster";}
	elsif ($species eq "Worm" or $species eq "worm"){$spec			   = "Caenorhabditis_elegans";}
	else{$spec														   = $species;}

	my ($host,$user,$pass);
	#####Connect with ENSEMBL mysql#####
	($host,$user,$pass) = split(",",$mysql) if (defined $mysql);
	$host = "ensembldb.ensembl.org" unless (defined $host);
	$user = "anonymous" unless (defined $user);
	$pass = "" unless (defined $pass);

	$dbtouse = lc($spec)."_variation_".$dbversion unless (defined $dbtouse);

	my $variationregistry = 'Bio::EnsEMBL::Registry';
	$variationregistry->load_registry_from_db(
		-host	    => $host,
		-user	    => $user,
		-pass	    => $pass,
		-species    => lc($spec),
		-db_version => $dbversion,
		-dbname	    => $dbtouse,
		-verbose    => 1
	    );

	$variationregistry->set_reconnect_when_lost();
	return $variationregistry;
}

sub connect_ontology_db{
	my $species   = shift||undef;
	my $dbversion = shift||undef;
	my $dbtouse   = shift||undef;
	my $mysql     = shift||undef;
	my $spec;

	if ($species eq "Human" or $species eq "human"){$spec			   = "Homo_sapiens";}
	elsif ($species eq "Mouse" or $species eq "mouse"){$spec		   = "Mus_musculus";}
	elsif ($species eq "Zebrafish" or $species eq "zebrafish"){$spec   = "Danio_rerio";}
	elsif ($species eq "Drosophila" or $species eq "drosophila"){$spec = "Drosophila_melanogaster";}
	elsif ($species eq "Worm" or $species eq "worm"){$spec			   = "Caenorhabditis_elegans";}
	else{$spec														   = $species;}

	my ($host,$user,$pass);
	#####Connect with ENSEMBL mysql#####
	($host,$user,$pass) = split(",",$mysql) if (defined $mysql);
	$host = "ensembldb.ensembl.org" unless (defined $host);
	$user = "anonymous" unless (defined $user);
	$pass = "" unless (defined $pass);

	$dbtouse = lc($spec)."_core_".$dbversion unless (defined $dbtouse);

	my $ontologyregistry = 'Bio::EnsEMBL::Registry';
	$ontologyregistry->load_registry_from_multiple_dbs(
		{
			    -host	    => $host,
				-user	    => $user,
				-pass	    => $pass,
				-species    => lc($spec),
				-db_version => $dbversion,
				-dbname	    => $dbtouse,
				-verbose    => 1
		},
		{
     			-host	    => $host,
				-user	    => $user,
				-pass	    => $pass,
				-dbname	    => "ensembl_ontology_".$dbversion,
				-verbose    => 1
		}
	    );

	$ontologyregistry->set_reconnect_when_lost();
	return $ontologyregistry;
}

sub parse_bed_linking{

	my $bed = shift;

	#temp vars
	my $ids = ();

	open(BED,"<:gzip(autopop)","$bed") || die ("Could not open $bed!\n");
	while(<BED>){
		chomp (my $raw = $_);
		push my @line , split (/\t/,$raw);
		(my $name	= $line[3])=~ s/_//g ||'u';
		push @$ids, $name;
	}
	close(BED);
	return(unique_array($ids));
}

sub parse_bedgraph{
	my $beds    = shift;
	my $anno    = shift;
	my $chr     = shift;

	#temp vars
	my $ret     = {};
	my $chl     = {};

	foreach my $bedfile (@{$beds}){

		open(my $BED,"<:gzip(autopop)","$bedfile") || die ("Could not open $bedfile!\n");

		while(<$BED>){
			chomp (my $raw = $_);
			push my @line , split (/\t/,$raw);
			(my $chromosome = $line[0]) =~ s/\_/=/g;
			$chromosome =~ s/\s/:/g;
			my $start = $line[1];
			my $end;
			$end = $line[2] if ($line[3]);
			$end = $start+1 unless (defined $end);#in case of per nucleotide bedgraph
			my $score;
			$score = $line[3] if ($line[3]);
			$score = $line[2] unless (defined $score);
			my $strand;
			$strand = $line[4] if (defined $line[4]);
			my $name   = 'u';
			if ($line[5]){
				$name = $line[3];
				$score  = $line[4];
				$strand  = $line[5];
			}
			$strand = 'u' unless (defined $strand);
			my $rest; # should this be more than just a bed 6 file
			if ($line[6]){
				($rest = join("\t",@line[6..$#line]))=~ s/\_/|/g;
			}
			else{
				$rest = "undef";
			}
			my $summit = "$score";
			$ret->{"$chromosome\_$start\_$end\_$strand"} = "$name\_$score\_$summit\_$rest" unless (defined $ret->{"$chromosome\_$start\_$end\_$strand"});

			$chl->{$chromosome}->{exists} = 1 unless (defined $chl->{$chromosome});
			#we need this info only if annotation is from mysql, maybe change that accordingly
			#			if (defined $mysql){
			#push @{$chl->{$chromosome}->{start}}, $start;
			#push @{$chl->{$chromosome}->{end}}  , $end;
		}
		close ($BED);
	}
	return($ret, $chl);
}

sub parse_Prepeaks{
	my $pre = shift;

	#tmpvars
	my @tmpmax;

	#return vars
	my $parsed = {};
	my $maxheight;

	foreach my $file (@{$pre}){
		my @forlimit;
		open(my $pk,"<:gzip(autopop)","$file") || die ("Could not open $file!\n");
		while(<$pk>){
			my @line     = split(/\t/,$_);
			(my $chrom = $line[0]) =~ s/\_/=/g;
			#	next if ($chrom eq "chrM"); ### No Interest in Mitochondria
			my $start   = $line[1];
			my $end     = $line[2];
#			my @profile = split(/\|/,$line[3]);
			my $pro     = $line[3];
			my $height  = $line[4];
			my $strand  = $line[5];
			my $rest; # should this be more than just a bed 6 file
			if ($line[6]){
				($rest = join("\t",@line[6..$#line]))=~ s/\_/|/g;
			}
			else{
				$rest = "undef";
			}
			$strand     = 'u' if ($strand ne '+' && $strand ne '-');
			push @forlimit, $height;
			$parsed->{"$chrom\_$strand\_$start\_$end"}="$pro\_$height\_$rest" unless (defined $parsed->{"$chrom\_$strand\_$start\_$end"});
		}
		my $maxh = max(@forlimit);
		push @tmpmax, $maxh;
	}
	$maxheight = max(@tmpmax);
	return($parsed, $maxheight);
}

sub parse_bed{

	my $beds    = shift;
	my $anno    = shift;
	my $chr     = shift;
	my $ignore  = shift;

	#temp vars
	my $ret     = {};
	my $chl     = {};

	foreach my $bedfile (@{$beds}){

		open(BED,"<:gzip(autopop)","$bedfile") || die ("Could not open $bedfile!\n");

		while(<BED>){
			chomp (my $raw = $_);
			$raw =~ s/[ ]+/\|/g;
			push my @line , split (/\t/,$raw);

			(my $chromosome = $line[0])=~ s/^chr//g;
			$chromosome =~ s/_/=/g;
			my $start	= $line[1];
			my $end		= $line[2];
			$end +=1 if ($end == $start); #Fix for broken bed format from peak callers
			(my $name	= $line[3])=~ s/_/=/g ||'u';
			my $score	= $line[4] || '0';
			my $strand  = $line[5] || 'u';
			$strand		= "u" if (($strand ne "-" && $strand ne "+") || ($ignore eq 'ON'));
			my $rest; # should this be more than just a bed 6 file

			if ($line[6]){
				($rest = join("\t",@line[6..$#line]))=~ s/\_/=/g;
			}
			else{
				$rest = "undef";
			}
			my $summit = $start;

			if ( $bedfile=~/\.pk/ || $bedfile=~/.spk/ ){
				$score  = $line[4];
				$summit = $line[6];
				if ($line[7]){
					($rest = join("\t",@line[7..$#line]))=~ s/\_/=/g;
#					$rest = join("\t",@line[7..$#line]);
				}
			}
			my $annotation = {};
			if(defined $anno && $anno eq 'ON'){

				next unless (defined $chr->{$chromosome});

				$annotation->{ID}     = $name;
				$annotation->{type}   = $rest || $score;
				$annotation->{start}  = $start;
				$annotation->{end}    = $end;
				$annotation->{strand} = $strand;

				$ret->{$chromosome}->{$strand} =  Set::IntervalTree->new unless (defined $ret->{$chromosome}->{$strand});
				$ret->{$chromosome}->{$strand} -> insert($annotation,$start-1,$end+1); #bed is 1 based at end only
			}
			else{
				$ret->{"$chromosome\_$start\_$end\_$strand"} = "$name\_$score\_$summit\_$rest";
			}
			$chl->{$chromosome}->{exists} = 1 unless (defined $chl->{$chromosome});
			#we need this info only if annotation if from mysql, maybe change that accordingly
			#			if (defined $mysql){
			push @{$chl->{$chromosome}->{start}}, $start;
			push @{$chl->{$chromosome}->{end}}  , $end;
			#			}
		}
	}
	close (BED);
	return($ret, $chl);
}

sub parse_annotation{

	my $file    = shift;
	my $check   = shift;
	my $chr     = shift;
	my $species = shift;
	my $host    = shift;
	my $dbv     = shift;
	my $usename = shift || undef;

	if ($file =~/.gff/ || $file =~ /.gtf/){
		open(my $gff,"<:gzip(autopop)",$file);
		chomp(my $line = <$gff>);
#		next unless ($line =~/^#/);
		if ($line =~ /gff-version\s+3/){
			my $out = parse_ensembl_gff3($file, $check, $chr, $usename);
			return($out);
		}
		else{
			my $out = parse_ensembl_gtf($file, $check, $chr);
			return($out);
		}
	}
	elsif($file =~/.bed/){
		my $annobed = [$file];
		my ($out,$tmp) = parse_bed($annobed, 'ON', $chr);
		return($out);
	}
	elsif($file eq 'mysql'){
		my ($out) = parse_mysql($chr, $check, $species, $host, $dbv );
		return($out);
	}
}

sub parse_expression{
	my $file   = shift;
	my $cutoff = shift;
	my $prog   = shift;
	my $unit   = shift;

	#return vars
	my $expression = {};
	open(my $exf,"<:gzip(autopop)",$file);

	if ($prog =~ /cufflinks/i){
		$expression = parse_cufflinks($exf)
	}

	if ($prog =~ /rnacounter/i){
		$expression = parse_rnacounter($exf)
	}

	if ($prog =~ /generic/i){
		$expression = parse_rnacounter($exf)
	}

	close($exf);

	if ($unit eq 'TPM'){
		return(filterTPM($expression),$cutoff);
	}
	elsif ($unit eq 'RPKM' || $unit eq 'FPKM'){
		return(filterTPM(FPKMtoTPM($expression),$cutoff));
	}
}

sub FPKMtoTPM{
	my $ex = shift;
	my $tpms = ();

	my $totalfpkm = sum(values %$ex);
	foreach my $trans (keys %$ex){
		$tpms->{$trans} = exp(log($ex->{$trans}) - log($totalfpkm) + log(1e6));
	}

	return($tpms);
}

sub filterTPM{
	my $ex = shift;
	my $cutoff = shift;

	foreach my $trans (keys %$ex){
		delete $ex->{$trans} if ($ex->{$trans} < $cutoff);
	}

	return($ex);
}

sub parse_cufflinks{
	my $file = shift;

	#return vars
	my $expression = {};

	print STDERR "Parsing Cufflinks output\n";

	while(<$file>){
		chomp($_);
		my @tmp = split("\t",$_);
		next if $tmp[2] ne 'transcript';
		my @ids = split("; ",$tmp[8]);
		my $idh = {};
		while(my $i = shift @ids){
#			next unless $i =~ /transcript_id/;
			my ($k,$v) = split(" ",$i);
			$v =~ s/\"//g;
			$idh->{$k}=$v;
		}
		$expression->{$idh->{'transcript_id'}}=$idh->{'FPKM'} unless ($idh->{'FPKM'} == 0);
	}
	return($expression);
}

sub parse_rnacounter{
	my $file = shift;

	#return vars
	my $expression = {};

	print STDERR "Parsing RNAcounter output\n";

	while(<$file>){
		chomp($_);
		my @tmp = split(' ',$_);
		$expression->{$tmp[0]}=$tmp[1] unless ($tmp[1] == 0);
	}
	return($expression);
}

sub parse_ensembl_gff3{

	my $file       = shift;
	my $check      = shift;
	my $chl        = shift;
	my $usename    = shift;
	my $regions    = {};

	open(my $gff,"<:gzip(autopop)",$file);

	while(<$gff>){
		chomp($_);
		next if ($_ =~ /^#/);

		my ($seqid, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split(/\t/,$_);
		$seqid =~ s/^chr//g;
		$seqid =~ s/_/=/g;
		$strand	= "u" unless ($strand eq "-" || $strand eq "+");

		###See if we need this region for annotation
		next unless (defined $chl->{$seqid}->{exists});

		### See if we have a match
		next unless (defined $check->{lc($type)});

		### Parse annotation information
		#Parent=transcript:ENST00000488147;Name=ENSE00003621279;constitutive=1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00003621279;rank=8;version=1
		my $annotation = {};
		my @a = split(';',$attributes);
		foreach(@a){
			my($k,$val) = split("=",$_);
			$annotation->{$k}=$val;
		}
		if ($annotation->{ID} && !$usename){
			if ($annotation->{ID} =~ /:/){
				$annotation->{ID} = (split(":",$annotation->{ID}))[1];
			}
			else{
				$annotation->{ID} = $annotation->{ID};
			}
		}
		else{
			if ($annotation->{Name} =~ /:/){
				$annotation->{ID} = (split(":",$annotation->{Name}))[1]||(split(":",$annotation->{Name}))[0];
			}
			else{
				$annotation->{ID} = $annotation->{Name};
			}
		}

		$annotation->{start}  = $start-1;
		$annotation->{end}    = $end;
		$annotation->{type}   = $type;#$annotation->{gene_biotype};
		$annotation->{strand} = $strand;
		$annotation->{anno}   = $annotation->{transcript_biotype};

		$regions->{$seqid}->{"$strand"} =  Set::IntervalTree->new unless (defined $regions->{$seqid}->{$strand});
		$regions->{$seqid}->{"$strand"} -> insert($annotation,$start-2,$end+1); #gff is 1 based at start and end
	}
	close($gff);

	return($regions);
}

sub extract_anno_from_gff{

	my $file       = shift;
	my $check      = shift;
	my $anno       = {};

	open(my $gff,"<:gzip(autopop)",$file);

	while(<$gff>){
		chomp($_);
		next if ($_ =~ /^#/);
		#		print $_."\n";
		my ($seqid, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split(/\t/,$_);
		($seqid) =~ s/chr//g;
		$strand	= "u" unless ($strand eq "-" || $strand eq "+");

		### Parse annotation information
		#Parent=transcript:ENST00000488147;Name=ENSE00003621279;constitutive=1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00003621279;rank=8;version=1
		my $annotation = {};
		my @a = split(';',$attributes);
		foreach(@a){
			my($k,$val) = split("=",$_);
			$annotation->{$k}=$val;
		}
		if ($annotation->{ID}){
			$annotation->{ID} = (split(":",$annotation->{ID}))[1];
		}
		else{
			$annotation->{ID} = (split(":",$annotation->{Name}))[1]||(split(":",$annotation->{Name}))[0];
		}

		$annotation->{start}  = $start-1;
		$annotation->{end}    = $end;
		$annotation->{type}   = $type;#$annotation->{gene_biotype};
		$annotation->{strand} = $strand;
		$annotation->{anno}   = $annotation->{transcript_biotype};

		### See if we have a match
		next unless (defined $check->{lc($type)} or defined $check->{lc($annotation->{transcript_biotype})});

		push @{$anno->{$type}}, $annotation->{ID};
		push @{$anno->{$annotation->{transcript_biotype}}} , $annotation->{ID};
	}
	close($gff);

	return($anno);
}

sub extract_anno_from_gtf{

	my $file       = shift;
	my $check      = shift;
	my $anno       = {};

	open(my $gff,"<:gzip(autopop)",$file);

	while(<$gff>){
		chomp($_);
		next if ($_ =~ /^#/);
		#		print $_."\n";
		my ($seqid, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split(/\t/,$_);
		($seqid) =~ s/chr//g;
		$strand	= "u" unless ($strand eq "-" || $strand eq "+");

		### Parse annotation information
		#Parent=transcript:ENST00000488147;Name=ENSE00003621279;constitutive=1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00003621279;rank=8;version=1
		my $annotation = {};
		my @a = split(/\;/,$attributes);
		foreach(@a){
			$_ =~ s/\s//g;
			my($k,$val) = split(/\"/,$_);
			$annotation->{$k}=$val;
		}
		if (!defined $annotation->{ID}){
			if ($type eq "transcript"){
				$annotation->{ID} = $annotation->{transcript_id};
			}
			elsif ($type eq "gene"){
				$annotation->{ID} = $annotation->{gene_id};
			}
			elsif ($type eq "exon"){
				$annotation->{ID} = $annotation->{exon_id};
			}
			else{
				$annotation->{ID} = $annotation->{gene_id};
			}
		}

		$annotation->{start}  = $start-1;
		$annotation->{end}    = $end;
		$annotation->{type}   = $type;#$annotation->{gene_biotype};
		$annotation->{strand} = $strand;
		$annotation->{anno}   = $annotation->{transcript_biotype};

#		print STDERR $annotation->{transcript_biotype};

		### See if we have a match
		next unless (defined $check->{lc($type)} or defined $check->{lc($annotation->{transcript_biotype})});

		push @{$anno->{$type}}, $annotation->{gene_id};
		push @{$anno->{$annotation->{transcript_biotype}}} , $annotation->{gene_id};

	}
	close($gff);

	return($anno);
}

sub parse_ensembl_gtf{

	my $file       = shift;
	my $check      = shift;
	my $chl        = shift;
	my $regions    = {};
	open(my $gtf,"<:gzip(autopop)",$file);

	while(<$gtf>){
		chomp($_);
		next if ($_ =~ /^#/);
		my ($seqid, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split(/\t/,$_);
		($seqid) =~ s/chr//g;
		$strand	= "u" unless ($strand eq "-" || $strand eq "+");

		###See if we need this region for annotation
		next unless (defined $chl->{$seqid}->{exists});

		### See if we have a match
		next unless (defined $check->{lc($type)});

		### Parse annotation information
		my $annotation = {};
		my @a = split(/\;/,$attributes);
		foreach(@a){
			$_ =~ s/\s//g;
			my($k,$val) = split(/\"/,$_);
			$annotation->{$k}=$val;
		}
		if (!defined $annotation->{ID}){
			if ($type eq "transcript"){
				$annotation->{ID} = $annotation->{transcript_id};
			}
			elsif ($type eq "gene"){
				$annotation->{ID} = $annotation->{gene_id};
			}
			elsif ($type eq "exon"){
				$annotation->{ID} = $annotation->{exon_id};
			}
			else{
				$annotation->{ID} = $annotation->{gene_id};
			}
		}

		$annotation->{start}  = $start-1;
		$annotation->{end}    = $end;
		$annotation->{type}   = $type;#$annotation->{gene_biotype};
		$annotation->{strand} = $strand;
		$annotation->{anno}   = $annotation->{transcript_biotype};

		$regions->{$seqid}->{"$strand"} =  Set::IntervalTree->new unless (defined $regions->{$seqid}->{$strand});
		$regions->{$seqid}->{"$strand"} -> insert($annotation,$start-2,$end+1); #gff is 1 based at start and end
	}
	close($gtf);

	return($regions);
}

sub fetch_gene_coords{

	my $species   = shift;
	my $mysql     = shift;
	my $dversion  = shift;
	my $coords    = shift;
	my $region    = "core";

	#return vars
	my $sizes={};

	##################
	# Load Registry
	##################
	my $spec;
	if ($species eq "Human" or $species eq "human"){$spec			   = "Homo_sapiens";}
	elsif ($species eq "Mouse" or $species eq "mouse"){$spec		   = "Mus_musculus";}
	elsif ($species eq "Zebrafish" or $species eq "zebrafish"){$spec   = "Danio_rerio";}
	elsif ($species eq "Drosophila" or $species eq "drosophila"){$spec = "Drosophila_melanogaster";}
	elsif ($species eq "Worm" or $species eq "worm"){$spec			   = "Caenorhabditis_elegans";}
	else{$spec														   = $species;}

	$dversion    = "79" unless ($dversion);
	my $dbtouse  = lc($spec)."_$region\_$dversion";
	my $registry = connect_db($spec, $dversion, $dbtouse, $mysql);

	my $slice_adaptor = $registry -> get_adaptor ("$species","$region",'Slice');
	my $gene_adaptor  = $registry -> get_adaptor ($spec,'core','gene');

	# get all non-redundant slices from the highest possible coordinate
	# systems
	my $slices = $slice_adaptor->fetch_all('chromosome');

	my $chromstoparse = {};
	foreach my $slice (@{$slices}){
		my $chrom  = $slice -> seq_region_name();
		$chromstoparse->{$chrom}=1;
	}

	$coords .= ".gz" unless ($coords =~/.gz$/);
	open(my $FH,">:gzip","$coords") or die "Could not open $coords to write gene coordinates to file!\n";
	my $idx = 1;

	foreach my $chr (sort {$a <=> $b} keys %{$chromstoparse}){
		print STDERR "Parsing chromosome $chr from ENSEMBL\n";
		my $slice = $slice_adaptor->fetch_by_region('chromosome',$chr);

		foreach my $gene ( @{$slice -> get_all_Genes} ) {

			my $genstr   = $gene->display_id();
			my $name     = $gene->external_name();
			my $start    = $gene->seq_region_start();
			my $end	     = $gene->seq_region_end();
			my $geneloci = $gene_adaptor -> fetch_by_stable_id("$genstr");
			my $chrom    = $geneloci -> seq_region_name();
			my $which    = $geneloci -> seq_region_strand();
			my $strand   = "+";
			$strand	     = "-" if ($which < 0);

			$sizes->{$genstr} = join(":", $chrom, $start, $end, $strand);

			my $line = join("\t",$idx,$genstr,join(":",$chrom,join("-",$start,$end),$strand));
			print $FH $line."\n";
			$idx++;
		}
	}
	close($FH);
	return($sizes);
}

sub parse_mysql{

	# get features for gene
	my $chl       = shift;
	my $check     = shift;
	my $species   = shift;
	my $mysql     = shift;
	my $dversion  = shift;
	my $region    = "core";

	# return values
	my $ret={};

	##################
	# Load Registry
	##################
	my $spec;
	if ($species eq "Human" or $species eq "human"){$spec			   = "Homo_sapiens";}
	elsif ($species eq "Mouse" or $species eq "mouse"){$spec		   = "Mus_musculus";}
	elsif ($species eq "Zebrafish" or $species eq "zebrafish"){$spec   = "Danio_rerio";}
	elsif ($species eq "Drosophila" or $species eq "drosophila"){$spec = "Drosophila_melanogaster";}
	elsif ($species eq "Worm" or $species eq "worm"){$spec			   = "Caenorhabditis_elegans";}
	else{$spec														   = $species;}

	$dversion    = "89" unless ($dversion);
	my $dbtouse  = lc($spec)."_$region\_$dversion";
	my $registry = connect_db($spec, $dversion, $dbtouse, $mysql);

	my $slice_adaptor = $registry -> get_adaptor ("$species","$region",'Slice');
	my $gene_adaptor  = $registry -> get_adaptor ($spec,'core','gene');

	foreach my $chr (sort {$a <=> $b} keys %{$chl}){
		next unless $chl->{$chr}->{exists};
		print STDERR "Parsing chromosome $chr from ENSEMBL\n";
		# now we fetch the slice we need, chromsome of interest from lowest start to highest end of bed to annotate
		for (0..$#{$chl->{$chr}->{start}}){
			my $slice  = $slice_adaptor->fetch_by_region("chromosome","$chr",${$chl->{$chr}->{start}}[$_],${$chl->{$chr}->{end}}[$_]);

			foreach my $raw_gene ( @{$slice -> get_all_Genes} ) {
				my $gene = $raw_gene->transfer($slice->seq_region_Slice());
				#			next unless ($gene->strand() && $gene->strand() eq $which) && $gene->display_id eq "$id";
				my $genstr   = $gene->display_id();
				my $name     = $gene->external_name();
				my $type     = $gene->biotype();

				### See if we have a match
				#				next unless ($check->{$type});

				my $start    = $gene->seq_region_start();
				my $end	     = $gene->seq_region_end();
				my $geneloci = $gene_adaptor -> fetch_by_stable_id("$genstr");
				my $chrom    = $geneloci -> seq_region_name();
				my $which    = $geneloci -> seq_region_strand();
				my $strand   = "+";
				$strand	     = "-" if ($which < 0);
				$strand	     = "u" unless ($strand eq "-" || $strand eq "+");

				my $annotation		  =	{};
				$annotation->{ID}     = $genstr;
				$annotation->{start}  = $start-1;
				$annotation->{end}    = $end;
				$annotation->{type}   = $type;
				$annotation->{strand} = $strand;

				$ret->{$chrom}->{$strand} =  Set::IntervalTree->new unless (defined $ret->{$chrom}->{$strand});
				$ret->{$chrom}->{$strand} -> insert($annotation,$start-2,$end+1); #ensembl-mysql is 1 based at start and end

				foreach my $raw_trans ( @{$gene -> get_all_Transcripts} ){
					my $trans   = $raw_trans->transfer($slice->seq_region_Slice());
					next unless ( $trans -> display_id() && $trans -> strand eq $which);
					my $type    = $trans-> biotype();
					my $transtr = $trans-> display_id() ;
					my $trst    = $trans->seq_region_start();
					my $tren    = $trans->seq_region_end();

					my $annotation        = {};
					$annotation->{ID}     = $transtr;
					$annotation->{start}  = $trst-1;
					$annotation->{end}    = $tren;
					$annotation->{type}   = $type;
					$annotation->{strand} = $strand;

#					$ret->{$chrom}->{$strand} =  Set::IntervalTree->new unless (defined $ret->{$chrom}->{$strand});
					$ret->{$chrom}->{$strand} -> insert($annotation,$trst-2,$tren+1); #ensembl-mysql is 1 based at start and end

					if (defined $trans->coding_region_start){
						my $cst = $trans->coding_region_start;
						my $cen = $trans->coding_region_end;

						my $annotation        = {};
						$annotation->{ID}     = $transtr.'_CDS';
						$annotation->{start}  = $cst-1;
						$annotation->{end}    = $cen;
						$annotation->{type}   = $type;
						$annotation->{strand} = $strand;

##						$ret->{$chrom}->{$strand} =  Set::IntervalTree->new unless (defined $ret->{$chrom}->{$strand});
						$ret->{$chrom}->{$strand} -> insert($annotation,$cst-2,$cen+1); #ensembl-mysql is 1 based at start and end

						my ( $five_start, $five_end, $three_start, $three_end ) = ( ''x4 );
						if ( $which == 1 ){
							if ($cst){
								if ( $cst != $trst ){
									$five_start = $trst;
									$five_end   = $cst-1;
								}
								if ( $cen != $tren ){
									$three_start = $cen+1;
									$three_end   = $tren;
								}
							}
						}
						else{
							if ($cen){
								if ( $cen != $tren ){
									$five_start = $cen+1;
									$five_end   = $tren;
								}
								if ( $cst != $trst ){
									$three_start = $trst;
									$three_end   = $cst-1;
								}
							}
						}

						if (defined $five_start && $five_start ne ''){
							my $annotation        = {};
							$annotation->{ID}     = $transtr.'_5UTR';
							$annotation->{start}  = $five_start-1;
							$annotation->{end}    = $five_end;
							$annotation->{type}   = $type;
							$annotation->{strand} = $strand;

##							$ret->{$chrom}->{$strand} =  Set::IntervalTree->new unless (defined $ret->{$chrom}->{$strand});
							$ret->{$chrom}->{$strand} -> insert($annotation,$five_start-2,$five_end+1); #ensembl-mysql is 1 based at start and end
						}

						if (defined $three_start && $three_start ne ''){
							my $annotation        = {};
							$annotation->{ID}     = $transtr.'_3UTR';
							$annotation->{start}  = $three_start-1;
							$annotation->{end}    = $three_end;
							$annotation->{type}   = $type;
							$annotation->{strand} = $strand;

#							$ret->{$chrom}->{$strand} =  Set::IntervalTree->new unless (defined $ret->{$chrom}->{$strand});
							$ret->{$chrom}->{$strand} -> insert($annotation,$three_start-2,$three_end+1); #ensembl-mysql is 1 based at start and end
						}
					}

					my $enc = 1;
					foreach my $exon (@{$trans -> get_all_Exons}){
						if ($exon && $exon ne "" && $exon != 0){
							my ($estart, $estop);
							$estart = $exon -> seq_region_start();
							$estop	= $exon -> seq_region_end();
							if ($estart && $estop && $estart ne "" && $estop ne "" && $estart != 0 && $estop != 0 && $estart ne 'NA' && $estop ne 'NA'){
								my $eid;
								my $elength = $exon->length();
								$eid = $exon -> stable_id() if (defined $exon -> stable_id());
								$eid = "$transtr\_Exon\_$enc" unless $eid;
								my $annotation        = {};
								$annotation->{ID}     = $eid;
								$annotation->{start}  = $estart-1;
								$annotation->{end}    = $estop;
								$annotation->{type}   = $type;
								$annotation->{strand} = $strand;
#								$ret->{$chrom}->{$strand} =  Set::IntervalTree->new unless (defined $ret->{$chrom}->{$strand});
								$ret->{$chrom}->{$strand} -> insert($annotation,$estart-2,$estop+1); #ensembl-mysql is 1 based at start and end
							}
						}
						$enc++;
					}

					my $inc = 1;
					foreach my $intron (@{$trans -> get_all_Introns()}){
						if ($intron && $intron ne "" && $intron != 0){
							my ($istart, $istop);
							$istart = $intron -> start();
							$istop	= $intron -> end();
							if ($istart ne "" && $istop ne "" && $istart != 0 &&$istop != 0){
								my $iid;
								$iid = "$transtr\_Intron\_$inc" unless $iid;
								my $annotation           = {};
								$annotation->{ID}     = $iid;
								$annotation->{start}  = $istart-1;
								$annotation->{end}    = $istop;
								$annotation->{type}   = $type;
								$annotation->{strand} = $strand;

#								$ret->{$chrom}->{$strand} =  Set::IntervalTree->new unless (defined $ret->{$chrom}->{$strand});
								$ret->{$chrom}->{$strand} -> insert($annotation,$istart-2,$istop+1); #ensembl-mysql is 1 based at start and end
							}
						}
						$inc++;
					}
				}
			}
		}
	}
	return($ret,$chl);
}

sub get_features_interval {

	my $intersects   = shift;
	my $coords       = shift;
	#	my $strand      = shift;

	# return values
	my (@strand,@id,@start,@end,@type);

	if (ref($intersects) eq 'HASH'){
		my ($str,$i,$sta,$en,$ty) = get_features($intersects);
		push @strand, @$str; push @id, @$i; push @start, @$sta; push @end, @$en; push @type, @$ty;
	}
	else{
		foreach my $annotation (@{$intersects}){
			if (ref($annotation) ne 'HASH'){
				foreach my $an (@$annotation){
					my ($str,$i,$sta,$en,$ty) = get_features($an);
					push @strand, @$str; push @id, @$i; push @start, @$sta; push @end, @$en; push @type, @$ty;
				}
			}
			else{
				my ($str,$i,$sta,$en,$ty) = get_features($annotation);
				push @strand, @$str; push @id, @$i; push @start, @$sta; push @end, @$en; push @type, @$ty;
			}
		}
	}
	if ($coords){
		return(\@start,\@end,\@type,\@id,\@strand);
	}
	else{
		push my @tmptype, map(split(/\|/,$_),@type);
		return(unique_coords_array(\@start,\@end),unique_array(\@tmptype),unique_array(\@id),unique_array(\@strand));
	}
	#	return(\@start,\@end,\@type,\@id,\@strand);
}

sub get_features{
	my $annotation = shift;

	# return values
	my ($strand,$id,$start,$end,$type)=();

	push @$id     , $annotation->{ID} || 'undef';
	push @$strand , $annotation->{strand} || 'undef';
	push @$start  , $annotation->{start} || 'undef';
	push @$end    , $annotation->{end} || 'undef';
	my @anno;
	push @anno, $annotation->{anno} if (defined $annotation->{anno});
	push @anno, $annotation->{type};
	my $anno = join('|', @anno);
#	my $tpe = $annotation->{type};
#	if ($anno !~ /$tpe/){
#		if ($anno eq 'undef'){
#			$anno = $tpe
#		}
#		else{
#			$anno .= '|'.$tpe
#		}
#	}
	push @$type , $anno;

	return ($strand,$id,$start,$end,$type);
}

sub unique_coords_array{

	my $start  = shift;
	my $end	   = shift;

	#return vars
	my ($s, $e) = ((),());

	#tmp vars
	my @coords;

	if (@$start){
		for (0..$#{$start}){
			push @coords, join("-",${$start}[$_],${$end}[$_]);
		}

		my %unique = ();
		foreach my $item (@coords)
		{
			$unique{$item} ++;
		}
		my @arrayuid = sort {$a cmp $b} keys %unique;

		for (0 ..$#arrayuid){
			(${$s}[$_], ${$e}[$_]) = split("-",$arrayuid[$_]);
		}
	}
	return($s,$e);
}

sub parse_Fasta{
	my $file		  = shift;
	my $lookforheader = shift;

	#RETURN
	my $fa;
	my ($header, $seq);
	while(<$file>){
		chomp(my $line = $_);

		if ($line =~ /Seq/){
			$header = $line;
		}
		else{
			$seq = $line;
		}
		if ($seq && $header){
			$fa->{$header} = $seq;
			($header, $seq) = undef;
		}
	}
	return($fa);
}

sub shuffle_seq{
	my $seq = shift;
	my $shufflestart = shift;
	my $shuffleend   = shift;
	my $shufflerounds = shift;

	my @seqa = split(//,$seq);
	$seq='';
	my @shuffledseq;
	$shufflestart = 1 unless ($shufflestart);
	$shuffleend = $#seqa+1 unless ($shuffleend);
	for ($shufflestart .. $shuffleend-1){
		push @shuffledseq , $seqa[$_];
	}

	for (1..$shufflerounds){
		@shuffledseq = shuffle(@shuffledseq);
	}
	for (0..$shufflestart-1){
		$seq .= $seqa[$_];
	}
	$seq .= join('',@shuffledseq);
	for ($shuffleend..$#seqa){
		$seq .= $seqa[$_];
	}

	return($seq);
}

sub add_constraints{
	my $seq = shift;
	my $constraintstart = shift;
	my $constraintlength   = shift || 8; #default is 8
	my $cons = shift || 'X'; # constraint type, default is make unpairable

	my @seqa = split(//,$seq);
	$seq='';
	my @constraintdseq;
	$constraintstart = 1 unless ($constraintstart);
	my $constraintend = $constraintstart+$constraintlength;
	for ($constraintstart .. $constraintend-1){
		push @constraintdseq , $cons;
	}
	for (0..$constraintstart-1){
		$seq .= $seqa[$_] if (($constraintstart-1) != 0);
	}
	$seq .= join('',@constraintdseq);
	for ($constraintend..$#seqa){
		$seq .= $seqa[$_];
	}
	return($seq);
}

sub parse_bed_plfold{
	my $file	  = shift;
	my $block	  = shift;
	my $rel		  = shift;
	my $gencoords = shift;

	### Parse Bed File
	open(BED,"<:gzip(autopop)",$file) || die ("Could not open $file!\n");
	my $bed={};

	while(<BED>){
        chomp (my $raw = $_);
        push my @line , split (/\t/,$raw);
        push @line, "\." if ( !$line[5] );
        (my $chromosome = $line[0])=~ s/chr//g;
        my $start		= $line[1]+1;
        my $end			= $line[2];
        my $name		= uc($line[3]);
        my $score		= $line[4];
		$score			= 0 if ($score eq $name);
        my $strand		= $line[5];
        my $rest		= join("|",@line[6 .. $#line]);
        $rest			= '' if ($rest eq "|");

		###To process bed and bed12 similar, we define two arrays and push entries there
		my @starts;
		my @ends;

		if ($block){
			my @blocksizes	= split(',',$line[10]);
			my @blockstarts = split(',',$line[11]);
			for (0 .. $#blockstarts){
				push @starts, $start+$blockstarts[$_];
				push @ends, $starts[$_]+$blocksizes[$_];
			}
		}
		else{
			push @starts, $start;
			push @ends,   $end;
		}

		#        if ($file	 =~ /\.pk/){
		#                $score	 =  $line[5];
		#                $strand =  $line[4];
		#        }

		#        if (($end-$start+1) < $minlength){
		#                $wobble = nearest(1,($minlength-($end-$start+1))/2);
		#        }

		while(@starts){

			my $wstart = shift @starts;
			my $wend   = shift @ends;

			unless ($rel){
				my ($gchrom,$gstart,$gend,$gstrand);
				($gchrom,$gstart,$gend,$gstrand) = split(/\:/,$gencoords->{$name}) if (defined $gencoords->{$name});
				($gchrom,$gstart,$gend,$gstrand) = split(/\:/,$gencoords->{$score}) if (defined $gencoords->{$score}); ### Often mixed up, so we check both
				warn "Could not find the gene identifier $name or $score from @line in the coordinates file!\n" unless (defined $gchrom && defined $gstrand);
				warn ("Strand of Peak @line, $strand and Gene $name $score: $gstrand are not the same, please check!!!\n") if ($gstrand ne $strand);

				if ($gstrand eq "+"){
					$wstart -= $gstart; #No +1 as we use 0-based indizes of array later
					$wend   -= $gstart;
				}
				elsif ($gstrand eq "-"){
					$wend	= $gend - $wstart;
					$wstart = $gend - $wend;
				}
				warn "Coordinate fuckup Gene:$gstart,$gend Peak:$wstart,$wend at $chromosome\:$start\:$end\:$name\:$score\:$strand\:$rest\n" if ($wstart < 1 || $wend < 1);
			}
			push @{$bed->{$name}}, "$chromosome\:$start\:$end\:$name\:$score\:$strand\:$rest\:$wstart\:$wend" unless ($wstart < 1 || $wend < 1 || defined $gencoords->{$score});
			push @{$bed->{$score}}, "$chromosome\:$start\:$end\:$name\:$score\:$strand\:$rest\:$wstart\:$wend" unless ($wstart < 1 || $wend < 1 || defined $gencoords->{$name});
		}
	}
	return($bed);
}

sub parse_plfold{

	my $genes	= shift;
	my $RT		= shift;
	my $pdir	= shift;
	my $ustart	= shift;
	my $uend	= shift;
	my $highmem = shift;

	###return vars
	my $folded={};

	foreach my $name (keys %{$genes}){
		print STDERR "Processing $name\n";
		open (L, "<:gzip(autopop)","$pdir\/$name\_lunp.gz") || die ("Could not open $pdir\/$name\_lunp.gz!");
		my %raw;
		if (defined $highmem){
			while(<L>){
				next if ($_ =~ /\#/);
				my @tmp = split('\t',$_);
				die ("Uend not in range of plfold prediction, please choose lower uend or rerun plfold with new params!\n") if ($uend > $#tmp);
				my $i = shift @tmp;
				$raw{$i} = join("\t",@tmp[$ustart-1 .. $uend -1]);
			}
			close(L);

			foreach my $peak (@{$genes->{$name}}){
				my ($chromosome,$start,$end,$gene,$score,$strand,$rest,$wstart,$wend) = split(/\:/,$peak);
				foreach my $pos ($wstart .. $wend){
					my $featurestart = $start+$pos;
					$folded->{$name}->{"Fold"}->{$pos} = $raw{$pos};
					$folded->{$name}->{"Peak"}->{$pos} = join("\t",$chromosome,$featurestart-1,$featurestart,$peak,$pos,$strand);
				}
			}
		}
		else{
			### Parsing on the fly, otherwise memory footprint too large, for cluster use alternative mode
			while(<L>){
				next if ($_ =~ /\#/);
				my @tmp = split('\t',$_);
				die ("Uend not in range of plfold prediction, please choose lower uend or rerun plfold with new params!\n") if ($uend > $#tmp);
				my $i = shift @tmp;
				foreach my $peak (@{$genes->{$name}}){
					#		print STDERR "$name $peak\n";
					my ($chromosome,$start,$end,$gene,$score,$strand,$rest,$wstart,$wend) = split(/\:/,$peak);
					next if ($wstart > $i || $wend < $i);
					my $featurestart = $start;
					foreach my $pos ($wstart .. $wend){
						next unless ($i == $pos);
						$featurestart+=$pos;
						$folded->{$name}->{"Fold"}->{$pos} = join("\t",@tmp[$ustart-1 .. $uend -1]);
						$folded->{$name}->{"Peak"}->{$pos} = join("\t",$chromosome,$featurestart-1,$featurestart,$peak,$pos,$strand);
					}
				}
			}
			close(L);
		}
	}
	return($folded);
}

sub unique_bed{
	my $bed = shift;

	##return vars
	my $unique={};
	my $totalreads=0;

	open (my $in, "<:gzip(autopop)" , $bed);

	while (<$in>){
		chomp (my $raw = $_);
		push my @line , split (/\t/,$raw);
		(my $chromosome = $line[0])=~ s/chr//g;
		$chromosome		= "M" if ($chromosome eq 'MT');
		my $peakstart   = $line[1];
		my $peakend     = $line[2];
		my $name        = $line[3] || 'u';
		my $score       = $line[4] || '0';
		my $strand      = $line[5] || 'u';
		$strand         = "u" if ($strand ne '+' && $strand ne '-');
		my $rest; # should this be more than just a bed 6 file

		if ($line[6]){
			($rest = join("\t",@line[6..$#line]))=~ s/\_/=/g;
		}
		else{
			$rest = "undef";
		}

		$totalreads++ unless (defined $unique->{"$chromosome\_$strand\_$peakstart\_$peakend\_$name"});
		$unique->{"$chromosome\_$strand\_$peakstart\_$peakend\_$name"} = "$score\_$rest" unless (defined $unique->{"$chromosome\_$strand\_$peakstart\_$peakend\_$name"});

	}
	close($in);

	return($unique, $totalreads);
}

sub uniqbed_to_coverage{
	my $unique = shift;
	my $anno   = shift;
	my $sizes  = shift;
	my $peak   = shift;
	my $conv   = shift;

	#return vars
	print STDERR "Transforming Bed with arguments ".join(";",$unique,$anno,$sizes,$peak,$conv)."\n";
	#return vars
	tie my %covplus, 'Tie::Hash::Indexed';
	tie my %covminus, 'Tie::Hash::Indexed';
	tie my %annop, 'Tie::Hash::Indexed';
	tie my %annom, 'Tie::Hash::Indexed';

	foreach my $peaks ( keys %{$unique} ){
		my @tmp= split(/_/,$peaks);
		push @tmp, split(/_/,$unique->{$peaks});
		my $chrom	= $tmp[0];
		(my $strand	= $tmp[1])=~s/\(|\)//g;
		my $cstart	= $tmp[2];
		my $cend	= $tmp[3];
		my $name	= $tmp[4];
		my $score	= $tmp[5];
		my $annotation; my @rest;
		if ($anno && $anno eq "anno"){
			($annotation = join("\t",@tmp[6..$#tmp]))=~ s/\_/=/g;
			$annotation  = '' if ($annotation eq 'undef');
			@rest	 = split(/\t/,$annotation);
		}

		if ( defined $peak && ($peak eq "peak")){
			my @pro;
			if ($name =~ /\:/){
				@pro = split(/\|/,$name) ;
			}
			else{
				@pro = split(/\|/,$rest[0]) if ($annotation && $rest[0] =~ /\:/);
			}
			my %poscov = ();
			if ($conv && $conv eq "on"){ ## If we use a split peak file
				my $nuk = $cstart-1;
				foreach my $pos (@pro){
					my @bla = split (/\:/,$pos);
					my $reg = $bla[0];
					for (1..$reg){
						$nuk++;
						next if ($nuk >= ($sizes->{$chrom}-1));
						$poscov{"$nuk"} = $bla[1] unless (defined $poscov{"$nuk"});
					}
				}
				for (my $i=$cstart;$i<=$cend;$i++){
					if ($strand eq "+" or $strand eq "1" or $strand eq "u" or $strand eq "=" or $strand eq "."){
						$covplus{"$chrom"}{"$i"} += $poscov{"$i"} if (defined $poscov{"$i"});
						$annop{"$chrom"}{"$i"} = $annotation if ($annotation);
					}
					elsif($strand eq "-" or $strand eq "-1"){
						$covminus{"$chrom"}{"$i"} += $poscov{"$i"} if (defined $poscov{"$i"});
						$annom{"$chrom"}{"$i"} = $annotation if ($annotation);
					}
				}
			}
			else{
				my $nu = $cstart;
				foreach my $pos (@pro){
					my @bla = split (/\:/,$pos);
					my $check = $bla[0]+$nu-1;
					next if (defined $poscov{"$check"});
					next if ($check >= ($sizes->{$chrom})-1);
					if ($bla[0] == 1){
						$poscov{"$nu"} = $bla[1] unless (defined $poscov{"$nu"});
						$nu++;
					}
					else{
						for (my $i = 1;$i <= $bla[0];$i++){
							$poscov{"$nu"} = $bla[1] unless (defined $poscov{"$nu"});
							$nu++;
						}
					}
				}

				for (my $i=$cstart;$i<=$cend;$i++){
					if ($strand eq "+" or $strand eq "1" or $strand eq "u" or $strand eq "=" or $strand eq "."){
						$covplus{"$chrom"}{"$i"}+=$poscov{"$i"} if (defined $poscov{"$i"});
						$annop{"$chrom"}{"$i"} = $annotation if ($annotation);
					}
					elsif($strand eq "-" or $strand eq "-1"){
						$covminus{"$chrom"}{"$i"}+=$poscov{"$i"} if (defined $poscov{"$i"});
						$annom{"$chrom"}{"$i"} = $annotation if ($annotation);
					}
				}
			}
		}
		else{
			for (my $i=$cstart;$i<=$cend;$i++){
				next if ($a > ($sizes->{$chrom}));
				if ($strand eq "+" or $strand eq "1" or $strand eq "u" or $strand eq "=" or $strand eq "."){
					$covplus{"$chrom"}{"$i"}+=1;
				}
				elsif($strand eq "-" or $strand eq "-1"){
					$covminus{"$chrom"}{"$i"}+=1;
				}
			}
		}
	}
	return(\%covplus, \%covminus, \%annop, \%annom);
}

sub bed_to_coverage{

	my $bed = shift;
	my $anno   = shift;
	my $sizes  = shift;
	my $peak   = shift;
	my $conv   = shift;
	my $totalreads = 0;

	print STDERR "Transforming Bed with arguments ".join(";",$bed,$anno,$sizes,$peak,$conv)."\n";
	#return vars
	tie my %covplus, 'Tie::Hash::Indexed';
	tie my %covminus, 'Tie::Hash::Indexed';
	tie my %annop, 'Tie::Hash::Indexed';
	tie my %annom, 'Tie::Hash::Indexed';

	my $in;
	my @lines;

	if ($bed =~ /.gz/){
		open ($in, "<:gzip" , $bed) or die "$!"." With file $bed\n";
		chomp(@lines = <$in>);
		close($in)
	}
	else{
		open ($in, "<:gzip(autopop)" , $bed) or die "$!";
		chomp(@lines = <$in>);
		close($in)
	}
	print STDERR "Read ".$#lines." lines from bed ".$bed."\n";
	$totalreads = $#lines;

	foreach my $raw (@lines){
		chomp($raw);
		my @line = split (/\t/,$raw);
		(my $chrom = $line[0])=~ s/^chr//g;
		$chrom		= "M" if ($chrom eq 'MT');
		my $cstart   = $line[1];
		my $cend     = $line[2];
		my $name        = $line[3] || 'u';
		my $score       = $line[4] || '0';
		(my $strand	= $line[5])=~s/\(|\)//g || 'u';
		$strand         = "u" if ($strand ne '+' && $strand ne '-');
		my $annotation; my @rest;
		if ($anno && $anno eq "anno"){
			($annotation = join("\t",@line[6..$#line]))=~ s/\_/=/g if ($line[6]);
			$annotation  = '' if ($annotation eq 'undef');
			@rest	 = split(/\t/,$annotation);
		}

		if ( defined $peak && ($peak eq "peak")){
			my @pro;
			if ($name =~ /\:/){
				@pro = split(/\|/,$name) ;
			}
            elsif($name eq 'X'){ #Special case for piranha
                my $peakwidth = $cend-$cstart;
                for (1..$peakwidth){
                    push @pro, "1:".nearest(1,$score/$peakwidth);
                }
            }
			else{
				@pro = split(/\|/,$rest[0]) if ($rest[0] =~ /\:/);
			}

			my %poscov = ();

			if ($conv && $conv eq "on"){ ## If coordinates already genomic
				my $nuk = $cstart-1;
				foreach my $pos (@pro){
					my @bla = split (/\:/,$pos);
					my $reg = $bla[0];
					for (1..$reg){
						$nuk++;
						next if ($nuk >= ($sizes->{$chrom}-1));
						$poscov{"$nuk"} = $bla[1] unless (defined $poscov{"$nuk"});
					}
				}
				for (my $i=$cstart;$i<=$cend;$i++){
					if ($strand eq "+" or $strand eq "1" or $strand eq "u" or $strand eq "=" or $strand eq "."){
						$covplus{"$chrom"}{"$i"} += $poscov{"$i"} if (defined $poscov{"$i"});
						$annop{"$chrom"}{"$i"} = $annotation if ($annotation);
					}
					elsif($strand eq "-" or $strand eq "-1"){
						$covminus{"$chrom"}{"$i"} += $poscov{"$i"} if (defined $poscov{"$i"});
						$annom{"$chrom"}{"$i"} = $annotation if ($annotation);
					}
				}
			}
			else{
				my $nu = $cstart;
				foreach my $pos (@pro){
					my @bla = split (/\:/,$pos);
					my $check = $bla[0]+$nu-1;
					if (defined $poscov{$check}){
						print STDERR "$check already defined";
						next;
					}
					if ($check >= ($sizes->{$chrom})-1){
						print STDERR "$check greater than size of $chrom with $sizes->{$chrom}";
						next;
					}
					if ($bla[0] == 1){
						$poscov{"$nu"} = $bla[1] unless (defined $poscov{"$nu"});
						$nu++;
					}
					else{
						for (my $i = 1;$i <= $bla[0];$i++){
							$poscov{"$nu"} = $bla[1] unless (defined $poscov{"$nu"});
							$nu++;
						}
					}
				}
				for (my $i=$cstart;$i<=$cend;$i++){
					if ($strand eq "+" or $strand eq "1" or $strand eq "u" or $strand eq "=" or $strand eq "."){
						$covplus{"$chrom"}{"$i"}+=$poscov{"$i"} if (defined $poscov{"$i"});
						$annop{"$chrom"}{"$i"} = $annotation if ($annotation);
					}
					elsif($strand eq "-" or $strand eq "-1"){
						$covminus{"$chrom"}{"$i"}+=$poscov{"$i"} if (defined $poscov{"$i"});
						$annom{"$chrom"}{"$i"} = $annotation if ($annotation);
					}
				}
			}
		}
        elsif (defined $peak && ($peak eq "score") ){
            for (my $i=$cstart;$i<$cend;$i++){
                next if ($i >= ($sizes->{$chrom}));
				if ($strand eq "+" or $strand eq "1" or $strand eq "u" or $strand eq "=" or $strand eq "."){
					$covplus{"$chrom"}{"$i"}+=$score if (defined $score);
					$annop{"$chrom"}{"$i"} = $annotation if ($annotation);
				}
				elsif($strand eq "-" or $strand eq "-1"){
					$covminus{"$chrom"}{"$i"}+=$score if (defined $score);
					$annom{"$chrom"}{"$i"} = $annotation if ($annotation);
				}
            }
        }
		else{
			for (my $i=$cstart;$i<=$cend;$i++){
				next if ($i >= ($sizes->{$chrom}));
				if ($strand eq "+" or $strand eq "1" or $strand eq "u" or $strand eq "=" or $strand eq "."){
					$covplus{"$chrom"}{"$i"}+=1;
					$annop{"$chrom"}{"$i"} = $annotation if ($annotation);
				}
				elsif($strand eq "-" or $strand eq "-1"){
					$covminus{"$chrom"}{"$i"}+=1;
					$annom{"$chrom"}{"$i"} = $annotation if ($annotation);
				}
			}
		}
	}
	return(\%covplus, \%covminus, \%annop, \%annom, $totalreads);

}

1;
