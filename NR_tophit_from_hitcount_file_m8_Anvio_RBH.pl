
#This program takes a list of blast output files and extracts the top match for each query.
#A file is produced with notes for how often each match is found.
#When multiple, equal hits are found, this record is used to identify the most frequently found hit.
#A report file returns notes on how assignments were made
#Results are in file [filename]_single_tophit.txt
#options -b bitcut -e eval -d tophit_database (name of tophit file)
#Please change the path for the HTI

use strict;
use Getopt::Long;

#my $hitcounts_path = "/Users/otte/NR_files";
my $hitcounts_path = "./";

my ($bit_cut,$e_cut,$output_name,$tophit_database,$RBH_file);

GetOptions(
	"b=i" => \$bit_cut,
	   "e=f" => \$e_cut,
	   "d=s" => \$tophit_database,
		"a=s" => \$RBH_file
);

if (not defined $e_cut) { $e_cut = 10; }
if (not defined $tophit_database) { print "\n Please specify database for tophits\n"; exit; }



my %RBH_assign;
if (defined $RBH_file) {
	print "Reading RBH file\n";
	open (RBH, "$RBH_file") || die "could not open RBH file $RBH_file\n";
	
	my $header = <RBH>;
	chomp $header;
	my @header = split (/\t/, $header);

	my $col=12;
	my $acc_col;

	while ($col <= 30) {
		if ($header[$col] =~ m/NCBI_PGAP_ACC/) {
			$acc_col = $col;
		} 
		$col++;
	}

	if (not defined $acc_col) {
		print "WARNING: ACC column not found\n";
	}
	
	while (my $line = <RBH>) {
		chomp $line;
		my @result = split (/\t/, $line);
		$result[$acc_col] =~ s/\.\d+//;
		$RBH_assign{$result[$acc_col]} = $result[1];
	}
	my $temp = scalar keys %RBH_assign;
	print "RBH file read, $temp acc found\n";
}


my %hitcounts;
my $total;
open (HITCOUNTS, "$hitcounts_path\/$tophit_database\_Hitcounts_all");
while (my $line = <HITCOUNTS>) {
	$total++;
	chomp $line;
	my @pep = split (/\t/, $line);	
	$hitcounts{$pep[0]} = $pep[1];
}
close HITCOUNTS;
if ($total > 0) {
	print "Hitcounts file parsed.  $total references found\n";
} else {
	print "Hitcounts file empty.  Starting new hitcounts file\n";
}

my %lookup_hits;
my %lookup_matches;
my %db_info;
my @databases;

print "\nFirst pass\n";
print "Queries\tHits  \tAssign\tNew ACC\tTotal ACC\n";

foreach my $input (@ARGV) {

	open (BLAST, "$input");

	my $m8_name = $input;
#	$m8_name =~ s/.*\///g;
	open (RESULTS, ">$input\_single_tophit.txt");
	
	my $total = 0;
	my $assigned = 0;
	my $found = 0;
	my %matches;
	my @best_match;
	my $num_matches;
	my $last;
	my $new_hits = 0;
	push (@databases, $m8_name);

	while (my $blast = <BLAST>) {
		
		chomp $blast;
		my @result = split (/\t/, $blast);
		
		if ($result[0] ne $last) {  # New hit!  Process last query.
			if ( $num_matches == 1 ) {	# One best hit - assign and add to hit file.
				my $hit_name = $best_match[1];
				$assigned++;
				$found++;
				unless (exists $hitcounts{$hit_name}) { $new_hits++; }
				$hitcounts{$hit_name}++;
				my $temp = join ("\t", @best_match);
				print RESULTS "$temp\n";
			} elsif ($num_matches > 1) {	# Multiple best hits - wait to make assignment.
				$found++;
				foreach my $hit ( keys %matches) {
					$lookup_hits{$hit}{'found'}++;
					$lookup_matches{$m8_name}{$last}{$hit} = $matches{$hit};	
				}
			}
			$num_matches = 0;	#Reset counters and assignment
			$last = $result[0];
		}

		next unless ($result[9] > $bit_cut && $result[10] < $e_cut );		
		my $query = $result[0];
		my $bit = $result[9];
		my $sign = $result[10];
		my $hit_name = $result[1];
		if ($hit_name =~ m/.*ref\|(.*?)\|.*/) { $hit_name = $1; }
		elsif ($hit_name =~ m/.*ref\|(.*?)$/) { $hit_name = $1; }
		elsif ($hit_name =~ m/.*\|(.*?)\|$/) { $hit_name = $1; }
		elsif ($hit_name =~ m/.*\|(.+)$/) { $hit_name = $1; }
		$hit_name =~ s/\.\d+//;		
		if (defined $RBH_file && !exists $RBH_assign{$hit_name}) {
			foreach my $compound_match ( $result[2] =~ m/\|(.*?)\|/g) {
				$compound_match =~ s/\.\d+//;
				if (exists $RBH_assign{$compound_match}) {
					$hit_name = $compound_match;
				}
			}
		}
		$result[1] = $hit_name;
			
		$total++;

		if ($best_match[0] ne $query) { # Initialize new best match.
			@best_match = @result;
			%matches = ( $hit_name => [ @result ] );
			$num_matches = 1;
			next;
		} elsif ($bit > $best_match[9] || $sign < $best_match[10]) {	# Better hit than best match - reinitialize
			@best_match = @result;
			%matches = ( $hit_name => [ @result ] );
			$num_matches = 1;
			next;
		} elsif ($bit == $best_match[9] && $sign == $best_match[10]) {	# Equal hit to best match - add to list
			$num_matches++;
			$matches{$hit_name} = [ @result ];
		}

	}

	if ( $num_matches == 1 ) {	# One best hit - assign and add to hit file.
		my $hit_name = $best_match[1];
		$assigned++;
		$found++;
		$hitcounts{$hit_name}++;
		my $temp = join ("\t", @best_match);
		print RESULTS "$temp\n";
	} elsif ($num_matches > 1) {	# Multiple best hits - wait to make assignment.
		$found++;
		foreach my $hit ( keys %matches) {
			$lookup_hits{$hit}{'found'}++;
			$lookup_matches{$m8_name}{$last}{$hit} = $matches{$hit};
		}
	}


	my $temp = scalar keys %hitcounts;
	print "$found\t$total\t$assigned\t$new_hits\t$temp\t$m8_name\n";
	
	$db_info{$m8_name}{'queries'} = $found;
	$db_info{$m8_name}{'assigned'} = $assigned;

	close BLAST;
	close RESULTS;
	
}

print "\nSecond Pass\n";
print "Assign\tTotal\tNew ACC\tTotal ACC\n";

foreach my $db (@databases) {
	
	my $last_numhits = scalar keys %hitcounts;

	unless (exists $lookup_matches{$db}) {
		print "0\t$db_info{$db}{'queries'}\t0\t$last_numhits\t$db\n";
		next;
	}
	
	open (RESULTS, ">>$db\_single_tophit.txt");

	
	my $count = 0;
	foreach my $query ( keys %{$lookup_matches{$db}} ) {
		
		$count++;
		my $num_matches = scalar keys %{$lookup_matches{$db}{$query}};
		my ($best_hit_name, $best_count);
		foreach my $hit_name (sort { $lookup_hits{$b}{'found'} <=> $lookup_hits{$a}{'found'} } keys %{$lookup_matches{$db}{$query}}) {
			if (not defined $best_hit_name) {
				$best_hit_name = $hit_name;
				$best_count = $hitcounts{$hit_name};
			} elsif ( exists $RBH_assign{$hit_name} && !exists $RBH_assign{$best_hit_name} ) {
				$best_hit_name = $hit_name;
				$best_count = $hitcounts{$hit_name};	
			} elsif ($hitcounts{$hit_name} > $best_count) {
				$best_hit_name = $hit_name;
				$best_count = $hitcounts{$hit_name};
			}
		}
		my $temp = join ("\t", @{$lookup_matches{$db}{$query}{$best_hit_name}});
		print RESULTS "$temp\n";
		delete $lookup_matches{$db}{$query};
		$hitcounts{$best_hit_name}++;
	}

	my $temp = (scalar keys %hitcounts) - $last_numhits;
	$last_numhits = scalar keys %hitcounts;
	print "$count\t$db_info{$db}{'queries'}\t$temp\t$last_numhits\t$db\n";
	
	close RESULTS;

}

open (HITCOUNTS, ">$hitcounts_path\/$tophit_database\_Hitcounts_all");
my %in_RBH;
foreach my $hit ( keys %hitcounts) {
	print HITCOUNTS "$hit\t$hitcounts{$hit}\n";
	if (exists $RBH_assign{$hit} ) {
		$in_RBH{'y'}{'acc'}++;
		$in_RBH{'y'}{'hits'}+= $hitcounts{$hit};
	} else {
		$in_RBH{'n'}{'acc'}{$hit} = 1;
		$in_RBH{'n'}{'hits'}+= $hitcounts{$hit};
	}
}
close HITCOUNTS;

if (defined $RBH_file && scalar keys %{$in_RBH{'n'}{'acc'}} >= 1) {
	my $temp = scalar keys %{$in_RBH{'n'}{'acc'}};
	print "\n$in_RBH{'y'}{'acc'} references with $in_RBH{'y'}{'hits'} hits found in RBH_file\n";
	print "$temp references with $in_RBH{'n'}{'hits'} hits not found\n";
	print "Unassigned Refs printed to file Unassigned_refs\n";
	open (OUT, ">Unassigned_refs");
	foreach my $acc ( keys %{$in_RBH{'n'}{'acc'}} ) {
		print OUT "$acc\t$hitcounts{$acc}\n";
	}
}

exit;





