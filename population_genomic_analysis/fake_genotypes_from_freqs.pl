# Usage:           perl script_FakeGenotypesFromFreqs.pl freqfile.tsv outfile.tsv N
#
#                  where freqfile.tsv = existing frequency file
#                        outfile.tsv  = output file to be created
#                        N            = number of tab-separated columns in output file
#
# Purpose:    Reads a tab-delimited file where the first column is a frequency
#             (ranging from 0 to 1), and for a user-specified N individuals, creates a
#             file where each corresponding line has N*freq 0's followed by N*(1-freq) 
#             1's, separated by tabs between every character. 
#             This was intended to create a "fake" tab-separated genotype file of 0's
#             and 1's using frequencies calculated from N genotypes per SNP.
#             Note that N*freq and N*(1-freq) should be integers (or nearly so), or
#             the number of 0 and 1 characters outputted will not always equal N.
#
# Run time:   3 seconds per 1,000,000 lines on a MacBook Pro (2018), 
#             with input argument N = 152
#
# Written by: Andy Gloss on Aug 7, 2021

# use strict;
use warnings;

# read input arguments
my $freqfile  = $ARGV[0];
my $outfile   = $ARGV[1];
my $num_indiv = $ARGV[2];

# open files
open (FREQFIL, "$freqfile") or die "Cannot open $freqfile\n";
open (OUTFIL,  ">$outfile") or die "Cannot create file: $outfile\n";

my $counter;

# read each line in the input file
while (my $line = <FREQFIL>) {

	# print update to the terminal every 1,000,000 lines
	$counter++;
	if ($counter =~ /000000$/) {
		print "Finished $counter SNPs...\n";
	}

	# separate contents of the line by tabs and store in an array, where the first
	# element ( $values[0] ) will be the value in the first "column" of the input file
	my @values = split('\t', $line);
	
	# if the frequency is "NA", print N tab-separated NA's to the output file
	if ($values[0] eq "NA") {	
		print OUTFIL "./.\t" x ($num_indiv - 1) . "./.\n";
		next;
	}
	
	# otherwise, print to the output file N*freq 0's and N*(1-freq) 1's, separated by tabs
	# (use sprintf to round to nearest integer in case frequencies are not exact)
	my $min_ct = sprintf( $num_indiv * $values[0] );
	my $maj_ct = sprintf( $num_indiv * (1 - $values[0]) );
	my $geno_string = "\t0|0" x $min_ct . "\t1|1" x $maj_ct;
	$geno_string =~ s/^\t//; # <-- removes the leading tab before the first value
	print OUTFIL "$geno_string\n";

	# to be safe: exit with error message if the number of output columns does not equal N
	my $num_columns = $geno_string =~ tr/\t// + 1;
	if ($num_columns ne $num_indiv) {
		print "Error, line $counter: N*freq and N*(1-freq) does not equal N\n"
	}

}

# close files
close (FREQFIL);
close (OUTFIL);

