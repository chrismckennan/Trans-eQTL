#! usr/bin/perl -w
use strict;

##This script will determine the consistency of Trans-eQTLs across tissues.
##The main variable is a double hash H with key1 = source and key2 = indirectly effected gene
##The value is the number of tissues this appears in

##The input is taken off of the command line:
##1: The posterior probability a gene is UNAFFECTED by SNP

my $pep = shift;
my $dir_in = '/Users/Chris/Desktop/Uchicago/StephensGroup/TranseQTLs_EnglehardtStephens/GitWork/Trans-eQTL/output/cis_summary_data';

opendir(INPUT, $dir_in) || die "Cannot read directory $dir_in: $!";
my @tissues = grep { (!/^\./) && -d "$dir_in/$_" } readdir INPUT;
closedir INPUT;

my %H;    #Hash to store all trans-eQTLs
foreach my $tissue (@tissues) {
	my $dir_tissue = "$dir_in/$tissue/thresh_80";
	
	opendir(INPUT, $dir_tissue) || die "Cannot read directory $dir_tissue: $!";
	my @files = grep { (!/^\./) && -f "$dir_tissue/$_" && /\.gibbs\.txt/ } readdir INPUT;
	closedir INPUT;
	
	foreach my $file (@files) {
		open (FH, "< $dir_tissue/$file") || die "Check $file\n";
		my @tmp_file = <FH>;
		close FH;
		
		my @source_gene_line = split("\t", $tmp_file[3]);
		my $source_gene = $source_gene_line[1];      #source gene
		chomp $source_gene;
		
		my $tmp_line6 = $tmp_file[6]; chomp $tmp_line6;
		my @network_genes = split("\t", $tmp_line6);
		my $tmp_line7 = $tmp_file[7]; chomp $tmp_line7;
		my @pep_genes = split("\t", $tmp_line7);
		my @trans_eQTLs_index = grep { $network_genes[$_] ne $source_gene && $pep_genes[$_] <= $pep } 0..$#pep_genes;
		my @trens_eQTLs = @network_genes[@trans_eQTLs_index];   #All trans-eQTLs
		
		foreach my $trans (@trens_eQTLs) {
			if (!defined( ${$H{$source_gene}{$trans}}[0] )) {
				@{$H{$source_gene}{$trans}} = ($tissue);
			} else {
				@{$H{$source_gene}{$trans}} = (@{$H{$source_gene}{$trans}}, $tissue);
			}
		}
	}
}

my $n_tissues = $#tissues + 1;
my @overlap = (0) x (2**$n_tissues - 1);
foreach my $key_source (keys( %H )) {
	foreach my $key_trans (keys( %{$H{$key_source}} )) {
		my @tmp_vec = @{ $H{$key_source}{$key_trans} };
		if ($#tmp_vec == 0) {
			my @ind = grep { $tissues[$_] eq $tmp_vec[0] } 0..$#tissues;
			$overlap[$ind[0]]++;
		}
		if ($#tmp_vec == 1) {
			if ( ($tmp_vec[0] =~ m/$tissues[0]/ || $tmp_vec[0] =~ m/$tissues[1]/) && ($tmp_vec[1] =~ m/$tissues[0]/ || $tmp_vec[1] =~ m/$tissues[1]/) ) {
				$overlap[3]++;   #tissues 0 and 1
			}
			if ( ($tmp_vec[0] =~ m/$tissues[0]/ || $tmp_vec[0] =~ m/$tissues[2]/) && ($tmp_vec[1] =~ m/$tissues[0]/ || $tmp_vec[1] =~ m/$tissues[2]/) ) {
				$overlap[4]++;   #tissues 0 and 2
			}
			if ( ($tmp_vec[0] =~ m/$tissues[1]/ || $tmp_vec[0] =~ m/$tissues[2]/) && ($tmp_vec[1] =~ m/$tissues[1]/ || $tmp_vec[1] =~ m/$tissues[2]/) ) {
				$overlap[5]++;   #tissues 1 and 2
			}
		}
		if ($#tmp_vec == 2) {
			$overlap[$#overlap]++;
			print STDOUT "$key_source     $key_trans\n";
		}
	}
}

my $output_file = "/Users/Chris/Desktop/Uchicago/StephensGroup/TranseQTLs_EnglehardtStephens/GitWork/Trans-eQTL/analysis/TissueCount_pep$pep.txt";
open(FH, "> $output_file");
for (my $i = 0; $i <= $#overlap; $i++) {
	if ($i <= 2) {
		print FH "$tissues[$i]\t";
	}
	if ($i == 3) {
		print FH "$tissues[0]_$tissues[1]\t";
	}
	if ($i == 4) {
		print FH "$tissues[0]_$tissues[2]\t";
	}
	if ($i == 5) {
		print FH "$tissues[1]_$tissues[2]\t";
	}	
	if ($i == $#overlap) {
		print FH "All3\t";
	}
	print FH "$overlap[$i]\n";
}

close FH;