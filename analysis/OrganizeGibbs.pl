#! usr/bin/perl -w
use strict;

##This program will organize the output files from the Gibbs sampler
##It simply takes the information in the header of the summary data in the 'data' folder and copies it to the corresponding Gibbs sampler file in the 'output' folder

##The first input is the summary data folder with tissue specific gene files
##The second input is the Gibbs sampler output folder with Gibbs sampler files

###Read Input###

my $Data_dir = $ARGV[0];
my $Output_dir = $ARGV[1];
if ($Data_dir =~ m/\/output\// || $Output_dir =~ m/\/data\//) { die "Re-order inputs"; }

opendir(INPUT, $Data_dir) || die "Cannot read directory $Data_dir: $!";
my @Data_files = grep { (!/^\./) && -f "$Data_dir/$_" } readdir INPUT;
@Data_files = sort @Data_files;
closedir INPUT;

opendir(OUTPUT, $Output_dir) || die "Cannot read directory $Output_dir: $!";
my @Output_files = grep { (!/^\./) && -f "$Output_dir/$_" } readdir OUTPUT;
@Output_files = sort @Output_files;
closedir OUTPUT;

###Concatenate header of summary data with Gibbs sampler output###

for (my $i = 0; $i <= $#Data_files; $i++) {
	my $sumdata = $Data_files[$i]; my $test = $sumdata;
	$test =~ s/(^[^\.]+\.[0-9]+).*/$1/;
	my $gibbsout = $Output_files[$i];
	unless ($gibbsout =~ m/$test/) {die "Problem with directory ordering in $test: $!";} #Make sure I am writing ot the correct file
	
	open (FH, "< $Output_dir/$gibbsout") or die "We're in trouble...: $!";
	my @gibbsout = <FH>;
	my $gibbsline = join('', @gibbsout);
	close FH;
	if ($gibbsline =~ m/^Tissue/) {next;}
	$gibbsline =~ s/[ ]+/\t/g;
	
	open (FH,  "< $Data_dir/$sumdata") or die "Cannot open $Data_dir/$sumdata: $!";
	my @sumdata = <FH>;
	@sumdata = @sumdata[0..4];
	my $sumline = join('', @sumdata);
	chomp $sumline;
	$sumline .= "\nPosterior probability each gene is UNAFFECTED by SNP:";
	close FH;
	
	open (FH, "> $Output_dir/$gibbsout");
	print FH "$sumline\n$gibbsline";
	close FH;
}