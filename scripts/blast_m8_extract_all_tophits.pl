#This takes files in Liz's blast_m8 format with all hits listed and extracts the top hits.  When an entry has multiple equal hits (by bit score), all equal hits are returned
#Use DIAMOND to create cutoff, then use this script to pull all hits from DIAMOND results

#Usage:
#	perl blast_m8_extract_all_tophits.pl [blast_m8_file list, try using *matches.m8] -b bitscore-cutoff -e ecut -s filesuffix (default "all_tophits")
#(5/13/20 @ 2:00pm) perl /Volumes/OttesenLab_Dropbox2/Protocols/ResultsMgmt_Scripts_fromLiz/blast_m8_extract_all_tophits.pl\
#					*matches.m8  --outfmt 6 qseqid sseqid salltitles pident qstart qend sstart send positive bitscore evalue length

#Expected format (tab-delimited):  Query, Hit ID, Hit description, percent identity, query start, query stop, subject start, subject stop, bitscore, e value, alignment length

use strict;
use Getopt::Long;

my $suffix;
my $e_cut;
my $bit_cut;

GetOptions(
	   "s=s" => \$suffix,
	   "e=f" => \$e_cut,
	   "b=f" => \$bit_cut
);

if (not defined $suffix) { $suffix = "alltophits"; }
if (not defined $e_cut) { $e_cut = 10; }

if ( @ARGV < 1 ) {
	print "This takes files in Liz's blast_m8 format with all hits listed and extracts the top hits.\n";
	print "When an entry has multiple equal hits (by bit score), all equal hits are returned\n";
	print "Usage $0 [list of blast_m8L_files] -b bitscore -e ecut -s file suffix (default alltophits)\n";
}

foreach my $input (@ARGV) {
	open (BLAST, $input);
	open (OUT, ">$input\_$suffix");

	my $last;
	my $counted;
	my $pass;
	my $bit;
	my $total;
	my $printed;
	my %results;

	while (my $line = <BLAST>) {
		chomp $line;
#		print $line;
#		exit;
		my @temp = split (/\t/, $line);
		print "$temp[0]\n";
		if ($temp[0] ne $last) {
			my $printed = 0;
			foreach my $print (sort { $results{$b} <=> $results{$a} } keys %results) {
				if ( $results{$print} >= $bit && $bit >= $bit_cut ) { 
					print OUT "$print\n"; 
					if ($printed == 0) { $pass++; $printed = 1; }
				}
			}
			%results = {};	
			$counted++;
			$last = $temp[0];
			$bit = $temp[9];
			print "$last\n$bit\n";
#			exit;
		}
		next unless ($temp[9] >= $bit && $temp[9] >= $bit_cut && $temp[10] <= $e_cut);
		if ($temp[9] > $bit) { $bit = $temp[9]; }
		$results{$line} = $temp[9];
	}

	my $printed = 0;
	foreach my $print (sort { $results{$b} <=> $results{$a} } keys %results) {
		if ( $results{$print} >= $bit && $bit >= $bit_cut ) { 
			print OUT "$print\n"; 
			if ($printed == 0) { $pass++; $printed = 1; }
		}
	}

	print "\nFile $input parsed\n";
	print "$counted queries counted, $pass queries had one or more hits above the cutoff\n";
	print "Results in File $input\_$suffix\n";
	
	close BLAST;
	close OUT;

}



exit;
