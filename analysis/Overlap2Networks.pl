#! usr/bin/perl -w
use strict;

#This code works!! It was tested on identical tissues and gave the correct answer.
#This will take two output directories with Gibbs sampler output files and compare networks across tissues. It returns the overlap across networks.
#Input: Two directories and trans eQTL threshold
#Output: 1 file giving overlap of gene-gene correlation networks and another file givng eQTL sharing
#The eQTL and network similarity thresholds are set here

my $eqtl_thresh = 0.2 ;    #Posterior probability gene is UNAFFECTED by SNP
my $network_fraction = 3/4;    #Fraction of network that must overlap to be considered the 'same' network

###Read Input###

my $Output_dir1 = $ARGV[0];
my $Output_dir2 = $ARGV[1];
my $T1 = $Output_dir1; $T1 =~ s/^.*cis_summary_data\/(.*)\/thresh_[0-9]+.*$/$1/;
my $T2 = $Output_dir2; $T2 =~ s/^.*cis_summary_data\/(.*)\/thresh_[0-9]+.*$/$1/;

opendir(INPUT, $Output_dir1) || die "Cannot read directory $Output_dir1: $!";
my @Output_files1 = grep { (!/^\./) && -f "$Output_dir1/$_" && m/gibbs/ } readdir INPUT;
@Output_files1 = sort @Output_files1;
closedir INPUT;

opendir(INPUT, $Output_dir2) || die "Cannot read directory $Output_dir2: $!";
my @Output_files2 = grep { (!/^\./) && -f "$Output_dir2/$_" && m/gibbs/  } readdir INPUT;
@Output_files2 = sort @Output_files2;
closedir INPUT;

###Create Initial Hashes###

my %Hcorr1; my %Heqtl1;
my %Hcorr2; my %Heqtl2;
foreach my $file1 (@Output_files1) {
  open (FH, "< $Output_dir1/$file1");
  my @lines1 = <FH>;
  close FH;
  
  my $geneline = $lines1[3]; chomp $geneline;
  my @geneline = split("\t", $geneline);
  my $gene = $geneline[1];    #This is a key for %H1

  my $eqtlline = $lines1[6]; chomp $eqtlline;
  my @eqtlline = split("\t", $eqtlline);
  my @corr_ind = grep{ $eqtlline[$_] ne $gene } (0..$#eqtlline);

  @{$Hcorr1{$gene}} = @eqtlline[@corr_ind];

  my $probs = $lines1[7]; chomp $probs;
  my @probs = split("\t", $probs);
  my @eqtl_ind = grep{ $probs[$_] <= $eqtl_thresh } grep{ $probs[$_] =~ m/[0-9]+/ } (0..$#probs);
  if (defined($eqtl_ind[0])) {
    @{$Heqtl1{$gene}} = @eqtlline[@eqtl_ind];
  }
}

foreach my $file2 (@Output_files2) {
  open (FH, "< $Output_dir2/$file2");
  my @lines2 = <FH>;
  close FH;

  my $geneline = $lines2[3]; chomp $geneline;
  my @geneline = split("\t", $geneline);
  my $gene = $geneline[1];    #This is a key for %H2

  my $eqtlline = $lines2[6]; chomp $eqtlline;
  my @eqtlline = split("\t", $eqtlline);
  my @corr_ind = grep{ $eqtlline[$_] ne $gene } (0..$#eqtlline);
  @{$Hcorr2{$gene}} = @eqtlline[@corr_ind];

  my $probs = $lines2[7]; chomp $probs;
  my @probs = split("\t", $probs);
  my @eqtl_ind = grep{ $probs[$_] <= $eqtl_thresh } grep{ $probs[$_] =~ m/[0-9]+/ } (0..$#probs);
  if (defined($eqtl_ind[0])) {
    @{$Heqtl2{$gene}} = @eqtlline[@eqtl_ind];
  }
}


##Look at overlap in both sets of hashes##
my $out_corr = "/Users/Chris/Desktop/Uchicago/StephensGroup/TranseQTLs_EnglehardtStephens/GitWork/Trans-eQTL/analysis/$T1"."_$T2.OverlapCorr.txt";
my $out_eqtl = "/Users/Chris/Desktop/Uchicago/StephensGroup/TranseQTLs_EnglehardtStephens/GitWork/Trans-eQTL/analysis/$T1"."_$T2.OverlapeQTL.txt";

##Correlation Overlap##


#1st Tissue#
my $overlap_eqtl = 0;
my $only1_eqtl = 0;
my $only2_eqtl = 0;
my %SharedeQTLs;    #The keys are genes with sufficient gene-gene correlation overlap. The values are number of eQTLs they share (if may have values of 0)

foreach my $key_corr1 (keys(%Hcorr1)) {
  if (defined($Hcorr2{$key_corr1}[0])) {
    my @genes1 = @{$Hcorr1{$key_corr1}};
    my @genes2 = @{$Hcorr2{$key_corr1}};
    my $overlap_corr = N_overlap(\@genes1, \@genes2);
    my $frac_network = $overlap_corr/(min($#genes1, $#genes2) + 1);

    if ($frac_network >= $network_fraction && defined($Heqtl1{$key_corr1}[0]) && defined($Heqtl2{$key_corr1}[0])) {
      $SharedeQTLs{$key_corr1}[0] = N_overlap(\@{$Heqtl1{$key_corr1}}, \@{$Heqtl2{$key_corr1}});
      $SharedeQTLs{$key_corr1}[1] = $#genes1 + 1;
      $SharedeQTLs{$key_corr1}[2] = $#genes2 + 1;
      $SharedeQTLs{$key_corr1}[3] = $#{$Heqtl1{$key_corr1}} + 1;
      $SharedeQTLs{$key_corr1}[4] = $#{$Heqtl2{$key_corr1}} + 1;
    } elsif ($frac_network >= $network_fraction) {
      $SharedeQTLs{$key_corr1}[0] = 0;
      $SharedeQTLs{$key_corr1}[1] = $#genes1 + 1;
      $SharedeQTLs{$key_corr1}[2] = $#genes2 + 1;
      $SharedeQTLs{$key_corr1}[3] = $#{$Heqtl1{$key_corr1}} + 1;
      $SharedeQTLs{$key_corr1}[4] = $#{$Heqtl2{$key_corr1}} + 1;
    }
  }
}

foreach my $key_eqtl1 (keys(%Heqtl1)) {
  my @genes1 = @{$Heqtl1{$key_eqtl1}};
  if (defined($Heqtl2{$key_eqtl1}[0])) {
    my @genes2 = @{$Heqtl2{$key_eqtl1}};
    $overlap_eqtl += N_overlap(\@genes1, \@genes2);
    $only1_eqtl += N_difference(\@genes1, \@genes2);
    $only2_eqtl += N_difference(\@genes2, \@genes1);
  } else {
    $only1_eqtl += $#genes1 + 1;
  }
}

foreach my $key_eqtl2 (keys(%Heqtl2)) {
  my @genes2 = @{$Heqtl2{$key_eqtl2}};
  if (!defined($Heqtl1{$key_eqtl2}[0])) {
    $only2_eqtl += $#genes2 + 1;
  }
}

open (FH, "> $out_eqtl");
print FH "SourceGene\tShared_trans-eQTLs\tNeighbors_T1\tNeighbors_T2\ttrans-eQTLs_T1\ttrans-eQTLs_T2\n";
foreach my $key (keys %SharedeQTLs) {
  print FH "$key\t";
  my $tmp_line = join("\t", @{$SharedeQTLs{$key}});
  print FH "$tmp_line\n";
}
print FH "#trans-eQTLs_Both\ttrans-eQTLs_T1\ttrans-eQTLs_T2\n";
print FH "#$overlap_eqtl\t$only1_eqtl\t$only2_eqtl";
close FH;



sub N_overlap {
  my ($v1, $v2) = @_;
  my @v1 = @{$v1};
  my @v2 = @{$v2};
  my %v1 = map{$_ => 1} @v1;
  my %v2 = map{$_ => 1} @v2;
  my @overlap_sub = grep( $v1{$_}, @v2 );
  if (defined($overlap_sub[0])) {
    return $#overlap_sub + 1;
  } else {
    return 0;
  }
}

sub N_difference {     #Finds things that are in 1, but NOT in 2
  my ($v1, $v2) = @_;
  my @v1 = @{$v1};
  my @v2 = @{$v2};
  my %v1 = map{$_ => 1} @v1;
  my %v2 = map{$_ => 1} @v2;
  my @notoverlap_sub = grep( !$v2{$_}, @v1 );
  if (defined($notoverlap_sub[0])) {
    return $#notoverlap_sub + 1;
  } else {
    return 0;
  }  
}

sub min {
  my @sorted_vec = sort @_;
  return $sorted_vec[0];
}

