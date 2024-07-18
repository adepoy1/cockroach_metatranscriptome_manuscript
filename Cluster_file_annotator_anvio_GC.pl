# Compile Tables.
# -h to enter path to hitcounts file 
# -o to enter path to output file
#  > will overwrite any existing output file
# annotation files should be remaining parameters, can handle any number of files with any number of annotation columns
# Assumes the first row of each annotation file is the header, and the first column is an accession number.

use strict;
use Getopt::Long;

my $hitcounts_file;
my $output_file; 

GetOptions(
	   "h=s" => \$hitcounts_file,
	   "o=s" => \$output_file
);


my %annotations;
my %headers;
my %blanks;
my %annot_cols;

foreach my $file (@ARGV) {
	open (INPUT, "$file");
	my $header1 = <INPUT>;
	chomp $header1;
	$header1 =~ s/\r+//;
	$header1 =~ s/.*?\t//;
	my @fields = split (/\t/, $header1);
	my $cols = @fields;
	my $blank;
	foreach my $temp (@fields) {
		$blank = $blank . "\t";
	}
	$headers{$file} = $header1;
	$blanks{$file} = $blank;
	$annot_cols{$file} = $cols;
	while (my $line = <INPUT>) {
		chomp $line;
		$line =~ s/\r+//;
		my ($acc, $annot) = $line =~ /(.*?)\t(.*)/;
		$annotations{$file}{$acc}{$annot} = 1;
	}

	close INPUT;
}

open (OUT, ">$output_file\_annotated.txt");

open (INPUT, "$hitcounts_file");

#generate and print header
my $header2 = <INPUT>;
chomp $header2;
my @header_split = split (/\t/, $header2);
my $header_end = @header_split;
$header_end = $header_end - 1;
my $h1 = join("\t", @header_split[0..14]);
my $h2 = join("\t", @header_split[15..$header_end]);
print OUT "$h1\t";
foreach my $file (@ARGV) { 
	print OUT "$headers{$file}\t";
}
print OUT "$h2";

my $found = 0;
my $not_found = 0;

while (my $line = <INPUT>) {
	chomp $line;
	$line =~ s/\r//;
	my @result = split (/\t/, $line);
	my @genes = split (/ /, $result[0]);
	
	print OUT "\n";
	print OUT join("\t", @result[0..14]);
	print OUT "\t";
	foreach my $file (@ARGV) {
		my %ac_annotations;
		foreach my $acc (@genes) {
			if (defined $annotations{$file}{$acc}) {
				foreach my $annot (keys %{$annotations{$file}{$acc}}) {
					$ac_annotations{$annot} = 1;
				}
			}	
		}
		my $entry;
		if ((scalar keys %ac_annotations) >1 ) {
			if ($annot_cols{$file} > 1) {
				
				foreach my $i (0..$annot_cols{$file}-1) {
					my %annot;
					foreach my $j (keys %ac_annotations) {
						my @j_cols = split (/\t/, $j);
						$annot{$j} = 1;
					}
					print OUT (keys %{$annot{$i}});
					print OUT "\t";
				}
			} else {
				print OUT (keys %ac_annotations);
				print OUT "\t";
			}
		} elsif ((scalar keys %ac_annotations) == 1 ) {
			print OUT (keys %ac_annotations);
			print OUT "\t";
		} else {
			print OUT "$blanks{$file}";
		}
	}
	my $result_len = @result;
	print OUT join("\t", @result[15..($result_len-1)]);

}

close INPUT;

close OUTPUT;

exit;