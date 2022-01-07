#!/usr/bin/perl

=head1
This script will sum up the frequency of reads inside each
chromosome, and inside windows of length specified by the
user. For example, for a size window of 100000 (100 kb) it
will sum up inside each chromosome all reads that map to 
the first 100kb, 200kb, and so on...

It uses the output of script 'parse_bowtie2.pl' as input 
It would be a '*prs' file

USAGE: perl 'sum_chr_known.pl' <directory> <extension> <interval> 

the adjective 'known' in the name of the script
derives from the fact the hg19 contains sequences
called chrM and chrUn; those are not considered
by this script now, since they constitute ~ 1%
=cut 

use lib "./";
use strict;
use ReadDir;
use warnings;

# Interactive mode
chomp (my $dir       = $ARGV[0]);
chomp (my $ext       = $ARGV[1]);
chomp (my $interval  = $ARGV[2]);

# qsub mode
#my $dir      = '/home/jjovel/jj/dataAnalysis/andy/andys_integration';
#my $ext      = 'prs';
#my $interval = 100000; 

my @files = ReadDir::get_dir_files($dir, $ext);
# Hash containing lenght of each chromosome

my (%human_genome) = (
	'chr1'  => 249250621, #chromosome1
	'chr2'  => 243199373, #chromosome2
	'chr3'  => 198022430, #chromosome3
	'chr4'  => 191154276, #chromosome4
	'chr5'  => 180915260, #chromosome5
	'chr6'  => 171115067, #chromosome6
	'chr7'  => 159138663, #chromosome7
	'chr8'  => 146364022, #chromosome8
	'chr9'  => 141213431, #chromosome9
	'chr10' => 135534747, #chromosome10
	'chr11' => 135006516, #chromosome11
	'chr12' => 133851895, #chromosome12
	'chr13' => 115169878, #chromosome13
	'chr14' => 107349540, #chromosome14
	'chr15' => 102531392, #chromosome15
	'chr16' =>  90354753, #chromosome16
	'chr17' =>  81195210, #chromosome17
	'chr18' =>  78077248, #chromosome18
	'chr19' =>  59128983, #chromosome19
	'chr20' =>  63025520, #chromosome20
	'chr21' =>  48129895, #chromosome21
	'chr22' =>  51304566, #chromosome22
	'chrX'  => 155270560, #chromosomeX
	'chrY'  =>  59373566, #chromosomeY
);

# Create a file containing the intervals 
# specified by the user
my $intervals_file = 'intervals.txt';
open (TEMPLATE, ">$dir/$intervals_file");

print "Creating intervals file: \n";
foreach my $key (sort keys %human_genome){
	for (my $window = $interval; $window <= $human_genome{$key}; $window += $interval){
		print TEMPLATE "$key\t", $window, "\n";
	}
}
close TEMPLATE;

foreach my $infile (@files){
	(my $outfile = $infile) =~ s/\.prs/\_$interval.tsv/;
	print "Processing file: $infile .....\n";
	open (IN, "$dir/$infile");
	open (OUT, ">$dir/$outfile");

	print OUT "Chr\tStart\tEnd\tValue\n";
	
	my %chromosomes = ();
	my $global_sum  = 0;	


	while (my $line = <IN>){

		chomp($line);
		my @temp       = split (/\t/, $line);
		(my $chr       = $temp[0]) =~ s/_.+//;
		my $chr_length = $human_genome{$chr}; 
		my $coord      = $temp[1];
		my $frequency  =  $temp[2] * 200;
		my $window   = 0;

		$global_sum += $temp[2];

		# create an interval
		if ($coord < $interval){$window = $interval}
		else {$window = (int ($coord / $interval) * $interval) + $interval}

		# fill hash of hashes of arrays containing frequency
		# per chromosome and per interval inside chromosomes
		unless (exists $chromosomes{$chr}{$window}){$chromosomes{$chr}{$window} = $frequency}
		else  {$chromosomes{$chr}{$window} += $frequency}

	}
	
	foreach my $prim_key (sort keys %chromosomes){
		foreach my $sec_key (sort keys %{$chromosomes{$prim_key}}){
			(my $new_key = $prim_key) =~ s/chr//;
			my $end = $sec_key + ($interval - 1);
			my $sum = sprintf("%.2f", $chromosomes{$prim_key}{$sec_key});
			print OUT "$new_key\t", $sec_key - $interval + 1, "\t", $sec_key, "\t$sum\n";
		}
	}
	
	# uncomment next line if you want to print total number of hits
	print "Total hits: $global_sum\n";

	close IN;
	close OUT;
}

exit;
