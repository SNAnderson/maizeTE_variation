#filter_blast_output_for_edge_reads_V6.pl by sna
use strict; use warnings;

die "usage: perl filter_blast_output_for_edge_reads_V6.pl <query_windows> <query_flanks_mapped_to_target> <suffix for output file names> <nested te file for query genotype> <query_genotype> <target_genotype> \n" unless @ARGV == 6;

# For all variables, b is query, w is target
my $geno_b = $ARGV[4];
my $geno_w = $ARGV[5];

my $edgeln = 400;

open(my $file1, $ARGV[0]) or die $!; # File 1 defines the window for each TE
my ($TE,$b_te_coordinates,$b_window_genes,$b_window_coordinates,$w_window_genes,$w_window_coordinates,$bchr,$bstartend,$bstart,$bend,$wchr,$wstartend,$wstart,$wend,$btechr,$btestartend,$btestart,$bteend,$b_gap);
my ($TEtmp,$TEfam,$TEind,$wchr1);
my %hash;

while (my $row = <$file1>){
  chomp $row;
  unless ($row =~ m/ctg/){ # contigs do not have annotated genes, so TEs are not anchored
    ($TEtmp,$b_te_coordinates,$b_window_genes,$b_window_coordinates,$w_window_genes,$w_window_coordinates) = split/\t/,$row;
    if ($TEtmp =~ m/ZM00001D/){
      $TEfam = substr $TEtmp,0,8;
      $TEind = substr $TEtmp,16;
      $TE = $TEfam.'Zm00001d'.$TEind;
    }
    else{
      $TE = $TEtmp;
    }
    unless($w_window_coordinates =~ m/NA/){
      ($bchr,$bstartend) = split/:/,$b_window_coordinates;
      ($bstart,$bend) = split/-/,$bstartend;
      ($wchr1,$wstartend) = split/:/,$w_window_coordinates;

      if ($wchr1 =~ m/chr/){
	$wchr = substr $wchr1,3;
      }
      else{
	$wchr = $wchr1;
      }
      
      ($wstart,$wend) = split/-/,$wstartend;
      
      ($btechr,$btestartend) = split/:/,$b_te_coordinates;
      ($btestart,$bteend) = split/-/,$btestartend;
      $b_gap = $bteend - $btestart + 1;

      unless ($wchr =~ m/scaffold/){
	$hash{$TE}{left_boundry} = $wstart;
	$hash{$TE}{right_boundry} = $wend;
	$hash{$TE}{window_length} = $wend - $wstart + 1;
	$hash{$TE}{b_window_length} = $bend - $bstart + 1;
	$hash{$TE}{chr} = $wchr;
	$hash{$TE}{b_window_genes} = $b_window_genes;
	$hash{$TE}{w_window_genes} = $w_window_genes;
	$hash{$TE}{bgap} = $b_gap;
	$hash{$TE}{bcoords} = $b_te_coordinates;
	$hash{$TE}{bwindowcoords} = $b_window_coordinates;
	$hash{$TE}{wwindowcoords} = $w_window_coordinates;
      }
    }
  }
}

close $file1;

open(my $file2, $ARGV[1]) or die $!; # File 2 is the BWA output for the 400 bp spanning the start and end coordinates of each TE
my (@line,@findgap,@remove_);
my ($te,$flank,$chr,$start,$end,$middle,$mapped_length,$matches,$mismatches);
my %hits;

while (my $row = <$file2>){
  unless (($row =~ m/qBeg/) or ($row =~ m/unmapped/) or ($row =~ m/ctg/) or ($row =~ m/scaffold/) or ($row =~ m/UNKNOWN/) or ($row =~ m/Mt/) or ($row =~ m/Pt/) or ($row =~ m/NCVQ/) or ($row =~ m/M99/)){
    chomp $row;
    #print "$row\n";
    @line = split/\t/,$row;
    @findgap = split/\./,$line[0];
    $flank = $findgap[1];

    if ($findgap[0] =~ m/\_/){
      @remove_ = split/\_/,$findgap[0];
      $te = $remove_[1];
    }
    else{
      $te = $findgap[0];
    }
    
    if (exists $hash{$te}{chr}){ #this is to filter out TEs on contigs

      $mapped_length = $line[10];
      if ($mapped_length/$edgeln >= 0.9 and $line[18] >= 0.9){ # check quality--At least 90% maps with at least 90% identity 
	if ($line[5] =~ m/chr/){
	  $chr = substr $line[5], 3;
	}

	# Fix Mo17 NCBI formatting
	elsif ($line[5] =~ m/CM009906.1/){
	  $chr = 1;
	}elsif ($line[5] =~ m/CM009907.1/){
	  $chr = 2;
	}elsif ($line[5] =~ m/CM009908.1/){
	  $chr = 3;
	}elsif ($line[5] =~ m/CM009909.1/){
	  $chr = 4;
	}elsif ($line[5] =~ m/CM009910.1/){
	  $chr = 5;
	}elsif ($line[5] =~ m/CM009911.1/){
	  $chr = 6;
	}elsif ($line[5] =~ m/CM009912.1/){
	  $chr = 7;
	}elsif ($line[5] =~ m/CM009913.1/){
	  $chr = 8;
	}elsif ($line[5] =~ m/CM009914.1/){
	  $chr = 9;
	}elsif ($line[5] =~ m/CM009915.1/){
	  $chr = 10;
	}

	# Fix Mo17 other formatting
	elsif ($line[5] =~ m/M0/){
	  $chr = substr $line[5], 2;
	}
	elsif ($line[5] =~ m/M10/){
	  $chr = 10;
	}
	
	else{
	  $chr = $line[5];
	}
	$start = $line[6];
	$end = $line[7];
	$middle = ($edgeln/2) - $line[1] + $line[6]; # determine expected coordinates of TE based on where the middle of the TE maps
	if ($flank ==1){
	  $middle = $middle + 1;
	}
	$matches = $line[11];
	$mismatches = $line[12];
	
	#Record number of total hits
	if (exists $hits{$te}{total_hits}{$flank}){
	  $hits{$te}{total_hits}{$flank}++;
	}
	else{
	  $hits{$te}{total_hits}{$flank} = 1;
	}
	
	#Check for hits within the right window
	if (($chr == $hash{$te}{chr}) and ($middle >= $hash{$te}{left_boundry}) and ($middle <= $hash{$te}{right_boundry} )){
	  if (exists $hits{$te}{$flank}){
	    $hits{$te}{$flank} = "multi";
	    $hits{$te}{window_hits}{$flank}++;
	  }
	  else{
	    $hits{$te}{$flank} = "$chr\t$start\t$middle\t$end\t$mapped_length\t$matches\t$mismatches";
	    $hits{$te}{window_hits}{$flank} = 1;
	  }
	}

	#Save hits even if in the wrong window
	else {
	  if (exists $hits{$te}{wrongwindow}{$flank}){
	    $hits{$te}{wrongwindow}{$flank} = "multi";
	    $hits{$te}{wrongwindow}{window_hits}{$flank}++;
	  }
	  else{
	    $hits{$te}{wrongwindow}{$flank} = "$chr\t$start\t$middle\t$end\t$mapped_length\t$matches\t$mismatches";
	    $hits{$te}{wrongwindow}{window_hits}{$flank} = 1;
	  }
	}
      }
      elsif($flank == 1 and $line[2] < 205 and $mapped_length > 180 and $line[18] >= 0.9){ # Left flank: if only the outside of the read maps, it might still be able to call a TE absence. Up to 5 bp of mapping into the TE is allowed.
	if ($line[5] =~ m/chr/){
	  $chr = substr $line[5], 3;
	}
	
	# Fix Mo17 NCBI formatting
	elsif ($line[5] =~ m/CM009906.1/){
	  $chr = 1;
	}elsif ($line[5] =~ m/CM009907.1/){
	  $chr = 2;
	}elsif ($line[5] =~ m/CM009908.1/){
	  $chr = 3;
	}elsif ($line[5] =~ m/CM009909.1/){
	  $chr = 4;
	}elsif ($line[5] =~ m/CM009910.1/){
	  $chr = 5;
	}elsif ($line[5] =~ m/CM009911.1/){
	  $chr = 6;
	}elsif ($line[5] =~ m/CM009912.1/){
	  $chr = 7;
	}elsif ($line[5] =~ m/CM009913.1/){
	  $chr = 8;
	}elsif ($line[5] =~ m/CM009914.1/){
	  $chr = 9;
	}elsif ($line[5] =~ m/CM009915.1/){
	  $chr = 10;
	}
	
	# Fix Mo17 other formatting
	elsif ($line[5] =~ m/M0/){
	  $chr = substr $line[5], 2;
	}
	elsif ($line[5] =~ m/M10/){
	  $chr = 10;
	}
	
	else{
	  $chr = $line[5];
	}
	$start = $line[6];
	$end = $line[7];
	$middle = ($edgeln/2) - $line[1] + $line[6];
	$matches = $line[11];
	$mismatches = $line[12];
	if (exists $hits{$te}{total_hits}{$flank}){
	  $hits{$te}{total_outside_hits}{$flank}++;
	}
	else{
	  $hits{$te}{total_outside_hits}{$flank} = 1;
	}
	
	#Check for hits within the right window
	if (($chr == $hash{$te}{chr}) and ($middle >= $hash{$te}{left_boundry}) and ($middle <= $hash{$te}{right_boundry} )){
	  if (exists $hits{$te}{outside}{$flank}){
	    $hits{$te}{outside}{$flank} = "multi";
	    $hits{$te}{window_hits}{outside}{$flank}++;
	  }
	  else{
	    $hits{$te}{outside}{$flank} = "$chr\t$start\t$middle\t$end\t$mapped_length\t$matches\t$mismatches";
	    $hits{$te}{window_hits}{outside}{$flank} = 1;
	  }
	}
      }
      elsif($flank == 2 and $line[1] > 195 and $mapped_length > 180 and $line[18] >= 0.9){ # Right flank: if only the outside of the read maps, it might still be able to call a TE absence. Up to 5 bp of mapping into the TE is allowed
	if ($line[5] =~ m/chr/){
	  $chr = substr $line[5], 3;
	}

	# Fix Mo17 NCBI formatting
	elsif ($line[5] =~ m/CM009906.1/){
	  $chr = 1;
	}elsif ($line[5] =~ m/CM009907.1/){
	  $chr = 2;
	}elsif ($line[5] =~ m/CM009908.1/){
	  $chr = 3;
	}elsif ($line[5] =~ m/CM009909.1/){
	  $chr = 4;
	}elsif ($line[5] =~ m/CM009910.1/){
	  $chr = 5;
	}elsif ($line[5] =~ m/CM009911.1/){
	  $chr = 6;
	}elsif ($line[5] =~ m/CM009912.1/){
	  $chr = 7;
	}elsif ($line[5] =~ m/CM009913.1/){
	  $chr = 8;
	}elsif ($line[5] =~ m/CM009914.1/){
	  $chr = 9;
	}elsif ($line[5] =~ m/CM009915.1/){
	  $chr = 10;
	}
	# Fix Mo17 other formatting
	elsif ($line[5] =~ m/M0/){
	  $chr = substr $line[5], 2;
	}
	elsif ($line[5] =~ m/M10/){
	  $chr = 10;
	}
	
	else{
	  $chr = $line[5];
	}
	$start = $line[6];
	$end = $line[7];
	$middle = ($edgeln/2) - $line[1] + $line[6];
	$matches = $line[11];
	$mismatches = $line[12];
	if (exists $hits{$te}{total_hits}{$flank}){
	  $hits{$te}{total_hits}{$flank}++;
	}
	else{
	  $hits{$te}{total_hits}{$flank} = 1;
	}
	
	#Check for hits within the right window
	if (($chr == $hash{$te}{chr}) and ($middle >= $hash{$te}{left_boundry}) and ($middle <= $hash{$te}{right_boundry} )){
	  if (exists $hits{$te}{outside}{$flank}){
	    $hits{$te}{outside}{$flank} = "multi";
	    $hits{$te}{window_hits}{outside}{$flank}++;
	  }
	  else{
	    $hits{$te}{outside}{$flank} = "$chr\t$start\t$middle\t$end\t$mapped_length\t$matches\t$mismatches";
	    $hits{$te}{window_hits}{outside}{$flank} = 1;
	  }
	}
      }
    }
  }
}

close $file2;

#Make hash of max TSD length per superfamily (Defined in Maize V4 genome paper), then later define te.gone when the gap between flanks is < (2 * max TSD)
my %maxTSD;
$maxTSD{RLC} = 6;
$maxTSD{RLG} = 6;
$maxTSD{RLX} = 6;
$maxTSD{DTT} = 2;
$maxTSD{DTA} = 8;
$maxTSD{DTM} = 9;
$maxTSD{DTH} = 3;
$maxTSD{DTC} = 3;
$maxTSD{DTX} = 9;
$maxTSD{RIL} = 15;
$maxTSD{RIT} = 15;
$maxTSD{DHH} = 1; # actually 0, but give a little buffer
$maxTSD{RST} = 15; 

# Create output files 
my ($out_all, $out_bed, $out_errors) = ("$geno_b.to.$geno_w.output1_$ARGV[2].txt","$geno_b.to.$geno_w.output2_$ARGV[2].bed","$geno_b.to.$geno_w.errors_$ARGV[2].txt");
my ($out1,$out2,$oute);
open ($out1,">",$out_all) or die "Couldn't open $out_all: $!";
open ($out2,">",$out_bed) or die "Couldn't open $out_bed: $!";
open ($oute,">",$out_errors) or die "Couldn't open $out_errors: $!";

# print header
print $out1 "$geno_b.TE\t$geno_b.coordinates\t$geno_b.anchors\t$geno_b.window\t$geno_b.window_length\t$geno_b.TE_length\t$geno_w.coordinates\t$geno_w.anchors\t$geno_w.window\t$geno_w.window_length\t$geno_w.gap_length\tclass_full\tclass_group\n";

# Remember which TEs are not defined as conserved or absent in the first pass so that TEs nested within absent TEs can be designated
my %find_site_gone;
my %poly_tes;

my ($chr1,$start1,$middle1,$end1,$mapped_length1,$matches1,$mismatches1,$chr2,$start2,$middle2,$end2,$mapped_length2,$matches2,$mismatches2);
my ($genesb, $genesw, $gapw,$gapstart,$gapend,$bgap,$bcoords,$superfam);
my ($window_length,$lefttot,$righttot,$leftgood,$rightgood);
for my $tes (keys %hash){
  
  #Get info from window hash
  $genesb = $hash{$tes}{b_window_genes};
  $genesw = $hash{$tes}{w_window_genes};
  $window_length = $hash{$tes}{window_length};
  $chr = $hash{$tes}{chr};
  $bgap = $hash{$tes}{bgap};
  $bcoords = $hash{$tes}{bcoords};
  $superfam = substr $tes, 0, 3;
  
  if(exists $hits{$tes}){ # At least one BWA read hit and passed filter

    #Pull info out of hits hash
    if (exists $hits{$tes}{total_hits}{1}){
      $lefttot = $hits{$tes}{total_hits}{1};
    }
    else {
      $lefttot = 0;
    }
    if (exists $hits{$tes}{total_hits}{2}){
      $righttot = $hits{$tes}{total_hits}{2};
    }
    else{
      $righttot = 0;
    }
    if (exists $hits{$tes}{window_hits}{1}){
      $leftgood = $hits{$tes}{window_hits}{1};
    }
    else{
      $leftgood = 0;
    }
    if (exists $hits{$tes}{window_hits}{2}){
      $rightgood = $hits{$tes}{window_hits}{2};
    }
    else{
      $rightgood = 0;
    }
	
    if ((exists $hits{$tes}{1}) and (exists $hits{$tes}{2})){ # Both left and right hits were recorded
      if(($hits{$tes}{1} !~ m/multi/) and ($hits{$tes}{2} !~ m/multi/)){
	($chr1,$start1,$middle1,$end1,$mapped_length1,$matches1,$mismatches1) = split/\t/,$hits{$tes}{1};
	($chr2,$start2,$middle2,$end2,$mapped_length2,$matches2,$mismatches2) = split/\t/,$hits{$tes}{2};
	$gapw = $middle2-$middle1 + 1;
	$gapstart = $middle1;
	$gapend = $middle2;

	# Print coordinates in an order that bedtools will accept
	if($gapend >= $gapstart){
	  print $out1 "$tes\t$bcoords\t$genesb\t$hash{$tes}{bwindowcoords}\t$hash{$tes}{b_window_length}\t$bgap\t$chr1:$gapstart-$gapend\t$genesw\t$hash{$tes}{wwindowcoords}\t$window_length\t$gapw\tconserved\tshared\n";
	  print $out2 "$chr1\t$gapstart\t$gapend\t$tes\t$bgap\t$gapw\t$bcoords\tconserved\n";
	}
	else{
	  print $out1 "$tes\t$bcoords\t$genesb\t$hash{$tes}{bwindowcoords}\t$hash{$tes}{b_window_length}\t$bgap\t$chr1:$gapstart-$gapend\t$genesw\t$hash{$tes}{wwindowcoords}\t$window_length\t$gapw\tconserved.negative.gap\tshared\n";
	  print $out2 "$chr1\t$gapend\t$gapstart\t$tes\t$bgap\t$gapw\t$bcoords\tconserved\n";
	}
	
      }
      elsif(($hits{$tes}{1} !~ m/multi/) and ($hits{$tes}{2} =~ m/multi/)){
       	($chr1,$start1,$middle1,$end1,$mapped_length1,$matches1,$mismatches1) = split/\t/,$hits{$tes}{1};
	$gapstart = $middle1;
	$gapend = $middle1;
	print $out1 "$tes\t$bcoords\t$genesb\t$hash{$tes}{bwindowcoords}\t$hash{$tes}{b_window_length}\t$bgap\t$chr1:$middle1-$middle1\t$genesw\t$hash{$tes}{wwindowcoords}\t$window_length\tNA\tconserved.left.flank.multi.right\tshared.left\n";
	print $out2 "$chr1\t$gapstart\t$gapend\t$tes.left\t$bgap\tUNK\t$bcoords\tconserved.left.flank\n";
      }
      elsif(($hits{$tes}{2} !~ m/multi/) and ($hits{$tes}{1} =~ m/multi/)){
	($chr2,$start2,$middle2,$end2,$mapped_length2,$matches2,$mismatches2) = split/\t/,$hits{$tes}{2};
	$gapstart = $middle2;
	$gapend = $middle2;
	print $out1 "$tes\t$bcoords\t$genesb\t$hash{$tes}{bwindowcoords}\t$hash{$tes}{b_window_length}\t$bgap\t$chr2:$middle2-$middle2\t$genesw\t$hash{$tes}{wwindowcoords}\t$window_length\tNA\tconserved.right.flank.multi.left\tshared.right\n";
	print $out2 "$chr2\t$gapstart\t$gapend\t$tes.right\t$bgap\tUNK\t$bcoords\tconserved.right.flank\n";
      }
      else{ #gap cannot be calculated due to multiple good hits
	$find_site_gone{$tes} =  "$tes\t$bcoords\t$genesb\t$hash{$tes}{bwindowcoords}\t$hash{$tes}{b_window_length}\t$bgap\tNA\t$genesw\t$hash{$tes}{wwindowcoords}\t$window_length\tNA\tmulti.good.hits\trejected\n";
      }
    }
    elsif ((exists $hits{$tes}{1}) or (exists $hits{$tes}{2})){ #gap cannot be calculated because there are not two good hits
      if((exists $hits{$tes}{1} and $hits{$tes}{1} =~ m/multi/) or (exists $hits{$tes}{2} and $hits{$tes}{2} =~ m/multi/)){
	$find_site_gone{$tes} =  "$tes\t$bcoords\t$genesb\t$hash{$tes}{bwindowcoords}\t$hash{$tes}{b_window_length}\t$bgap\tNA\t$genesw\t$hash{$tes}{wwindowcoords}\t$window_length\tNA\tone.flank.multi.good.hits\trejected\n";
      }
      else{
	if (exists $hits{$tes}{1}){
	  ($chr1,$start1,$middle1,$end1,$mapped_length1,$matches1,$mismatches1) = split/\t/,$hits{$tes}{1};
	  $gapstart = $middle1;
	  $gapend = $middle1;
	  print $out1 "$tes\t$bcoords\t$genesb\t$hash{$tes}{bwindowcoords}\t$hash{$tes}{b_window_length}\t$bgap\t$chr1:$middle1-$middle1\t$genesw\t$hash{$tes}{wwindowcoords}\t$window_length\tNA\tconserved.left.flank\tshared.left\n";
	  print $out2 "$chr1\t$gapstart\t$gapend\t$tes.left\t$bgap\tUNK\t$bcoords\tconserved.left.flank\n";
	}
	else{
	  ($chr2,$start2,$middle2,$end2,$mapped_length2,$matches2,$mismatches2) = split/\t/,$hits{$tes}{2};
	  $gapstart = $middle2;
	  $gapend = $middle2;
	  print $out1 "$tes\t$bcoords\t$genesb\t$hash{$tes}{bwindowcoords}\t$hash{$tes}{b_window_length}\t$bgap\t$chr2:$middle2-$middle2\t$genesw\t$hash{$tes}{wwindowcoords}\t$window_length\tNA\tconserved.right.flank\tshared.right\n";
	  print $out2 "$chr2\t$gapstart\t$gapend\t$tes.right\t$bgap\tUNK\t$bcoords\tconserved.right.flank\n";
	}
      }
    }
    elsif ((exists $hits{$tes}{outside}{1}) and (exists $hits{$tes}{outside}{2})){ # check for absent TE only when full-length pseud-reads do not hit
      if(($hits{$tes}{outside}{1} !~ m/multi/) and ($hits{$tes}{outside}{2} !~ m/multi/)){
	($chr1,$start1,$middle1,$end1,$mapped_length1,$matches1,$mismatches1) = split/\t/,$hits{$tes}{outside}{1};
	($chr2,$start2,$middle2,$end2,$mapped_length2,$matches2,$mismatches2) = split/\t/,$hits{$tes}{outside}{2};
	$gapw = $middle2-$middle1;
	if ($middle1 < $middle2){
	  $gapstart = $middle1;
	  $gapend = $middle2;
	}
	else{
	  $gapstart = $middle2;
	  $gapend = $middle1;
	}
	if (abs($gapw) <= (2 * $maxTSD{$superfam})){ # empty site should lack the TSD length, and deleted TEs should probably be imperfectly excised
	  print $out1 "$tes\t$bcoords\t$genesb\t$hash{$tes}{bwindowcoords}\t$hash{$tes}{b_window_length}\t$bgap\t$chr1:$gapstart-$gapend\t$genesw\t$hash{$tes}{wwindowcoords}\t$window_length\t$gapw\tte.gone\tnon.shared\n";
	  print $out2 "$chr1\t$gapstart\t$gapend\t$tes.site\t$bgap\t$gapw\t$bcoords\tpoly.te.site\n";
	  $poly_tes{$tes} = "TRUE";
	}
	else{
	  $find_site_gone{$tes} = "$tes\t$bcoords\t$genesb\t$hash{$tes}{bwindowcoords}\t$hash{$tes}{b_window_length}\t$bgap\t$chr1:$gapstart-$gapend\t$genesw\t$hash{$tes}{wwindowcoords}\t$window_length\t$gapw\ttruncated.gap.too.big\trejected\n";
	}
      }
      else {
	$find_site_gone{$tes} = "$tes\t$bcoords\t$genesb\t$hash{$tes}{bwindowcoords}\t$hash{$tes}{b_window_length}\t$bgap\tNA\t$genesw\t$hash{$tes}{wwindowcoords}\t$window_length\tNA\tmulti.good.truncated.hits\trejected\n";
      }
    }
    else { # There are no good hits, so check for single hits outside of window.
      if ((exists $hits{$tes}{wrongwindow}{1}) and (exists $hits{$tes}{wrongwindow}{2})){ # Both left and right hits were recorded
	if(($hits{$tes}{wrongwindow}{1} !~ m/multi/) and ($hits{$tes}{wrongwindow}{2} !~ m/multi/)){
	  ($chr1,$start1,$middle1,$end1,$mapped_length1,$matches1,$mismatches1) = split/\t/,$hits{$tes}{wrongwindow}{1};
	  ($chr2,$start2,$middle2,$end2,$mapped_length2,$matches2,$mismatches2) = split/\t/,$hits{$tes}{wrongwindow}{2};
	  $gapw = $middle2-$middle1 + 1;
	  $gapstart = $middle1;
	  $gapend = $middle2;
	  
	  if ($chr1 == $chr2 and abs($gapw) < 1000000){
	    if ($chr1 == $hash{$tes}{chr}){
	      if($gapw >= 0){
		$find_site_gone{$tes} = "$tes\t$bcoords\t$genesb\t$hash{$tes}{bwindowcoords}\t$hash{$tes}{b_window_length}\t$bgap\t$chr1:$gapstart-$gapend\t$genesw\t$hash{$tes}{wwindowcoords}\t$window_length\t$gapw\twrong.window.right.chromosome\trejected\n";
	      }
	      else{
		$find_site_gone{$tes} = "$tes\t$bcoords\t$genesb\t$hash{$tes}{bwindowcoords}\t$hash{$tes}{b_window_length}\t$bgap\t$chr1:$gapstart-$gapend\t$genesw\t$hash{$tes}{wwindowcoords}\t$window_length\t$gapw\twrong.window.negative.gap.right.chromosome\trejected\n";
	      }
	    }
	    else{
	      if($gapw >= 0){
		$find_site_gone{$tes} = "$tes\t$bcoords\t$genesb\t$hash{$tes}{bwindowcoords}\t$hash{$tes}{b_window_length}\t$bgap\t$chr1:$gapstart-$gapend\t$genesw\t$hash{$tes}{wwindowcoords}\t$window_length\t$gapw\twrong.window.wrong.chromosome\trejected\n";
	      }
	      else{
		$find_site_gone{$tes} = "$tes\t$bcoords\t$genesb\t$hash{$tes}{bwindowcoords}\t$hash{$tes}{b_window_length}\t$bgap\t$chr1:$gapstart-$gapend\t$genesw\t$hash{$tes}{wwindowcoords}\t$window_length\t$gapw\twrong.window.negative.gap.wrong.chromosome\trejected\n";
	      }
	    }
	  }
	  else{
	    $find_site_gone{$tes} = "$tes\t$bcoords\t$genesb\t$hash{$tes}{bwindowcoords}\t$hash{$tes}{b_window_length}\t$bgap\tNA\t$genesw\t$hash{$tes}{wwindowcoords}\t$window_length\tNA\trejected.wrong.window.hits\trejected\n";
	  }
	}
	else{ #gap cannot be calculated due to multiple good hits
	  $find_site_gone{$tes} = "$tes\t$bcoords\t$genesb\t$hash{$tes}{bwindowcoords}\t$hash{$tes}{b_window_length}\t$bgap\tNA\t$genesw\t$hash{$tes}{wwindowcoords}\t$window_length\tNA\trejected.wrong.window.hits\trejected\n";
	}
      }
      else{
	$find_site_gone{$tes} = "$tes\t$bcoords\t$genesb\t$hash{$tes}{bwindowcoords}\t$hash{$tes}{b_window_length}\t$bgap\tNA\t$genesw\t$hash{$tes}{wwindowcoords}\t$window_length\tNA\tno.good.hits\trejected\n";
      }
    }
  }
  else{ #no hits in the file
    if (exists $hash{$tes}{chr}){
      $find_site_gone{$tes} = "$tes\t$bcoords\t$genesb\t$hash{$tes}{bwindowcoords}\t$hash{$tes}{b_window_length}\t$bgap\tNA\t$genesw\t$hash{$tes}{wwindowcoords}\t$window_length\tNA\tno.hits\trejected\n";
    }
  }
}

open(my $file3, $ARGV[3]) or die $!; # File 3 has every TE that is nested inside of another TE. If outer TE is defined as absent, inner TE must also be gone.
my ($lilte,$lilc,$bigte,$bigc);
my @printline;

while (my $row = <$file3>){
  chomp $row;
  unless ($row=~m/nested_te/){
    ($lilte,$lilc,$bigte,$bigc) = split/\t/,$row;
    if (exists $poly_tes{$bigte}){
      if (exists $find_site_gone{$lilte}){
	@printline = split/\t/,$find_site_gone{$lilte};
	$printline[-2] = "te.gone.site.gone";
	$printline[-1] = "non.shared";
	$find_site_gone{$lilte} = join("\t",@printline);
      }
      else{
	print $oute "check.call\t$lilte in $bigte\n";
      }
    }
  }
}

close $file3;

for my $reject_tes (keys %find_site_gone){
  chomp $find_site_gone{$reject_tes};
  print $out1 "$find_site_gone{$reject_tes}\n";
}
