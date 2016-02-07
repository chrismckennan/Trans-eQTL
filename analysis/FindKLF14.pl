#! usr/bin/perl -w
use strict;

###Look through a list of gene files to find genes with neighbor KLF14: ENSG00000266265
###KLF14 has been shown to be a master regulator in adipose tissue
###There is no network with anchor gene ENSG00000266265

###The input is a directory with output files (i.e. ../output/..)
###The output is a list of file names with ENSG00000266265 as a neighbor, located in the same directory
###The file name is ENSG00000266265.neigh.txt

my $gene1 = 'ENSG00000266265';
my $gene2 = 'ENSG00000266265.2';
my $dir = shift @ARGV;

opendir(INPUT, $dir) || die "Cannot read directory $dir: $!";
my @dir_files = grep { (!/^\./) && -f "$dir/$_" } readdir INPUT;
closedir INPUT;

my %files;
for (my $i = 0; $i <= $#dir_files; $i++) {
	open (FH, "< $dir/$dir_files[$i]") || die "Check $dir_files[$i]\n";
	my @tmp_file = <FH>;
	close FH;
	
	my $line = $tmp_file[6]; unless (defined $line) {next;}
	chomp $line;
	my @line_array = split "\t", $line;
	my $line2 = $tmp_file[7];
	chomp $line2;
	my @line2_array = split "\t", $line2;
	if ($line =~ m/$gene1/ || $line =~ m/$gene2/) {
		my @gene_index = grep {$line_array[$_] =~ /$gene1/ || $line_array[$_] =~ /$gene2/} 0..$#line_array;
		my $genes =  $line_array[$gene_index[0]];
		my $prob = $line2_array[$gene_index[0]];
		@{$files{$dir_files[$i]}} = ($genes, $prob);
	}
}

open (OUT, "> $dir/ENSG00000266265.neigh.txt") || die;
foreach my $key (sort(keys %files)) {
	my $tmp = join("\t", @{$files{$key}});
	print OUT "$key\t$tmp\n";
}
close OUT;
