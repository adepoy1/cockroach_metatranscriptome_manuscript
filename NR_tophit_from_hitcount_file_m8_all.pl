
#This program takes a list of blast output files and extracts the top match for each query.
#A file is produced with notes for how often each match is found.
#When multiple, equal hits are found, this record is used to identify the most frequently found hit.
#A report file returns notes on how assignments were made
#Results are in file [filename]_single_tophit.txt
#options -b bit_cutoff -i ID_cutoff -a aln_cutoff -e eval -d tophit_database (name of tophit file)
#Please change the path for the HTI

use strict;
use Getopt::Long;

#my $hitcounts_path = "/Users/otte/Bioinformatics/NR_files";
my $hitcounts_path = "./";

my ($length_cut,$bit_cut,$e_cut,$id_cut,$num,$aln_cut,$output_name,$tophit_database);

GetOptions(
	   "i=i" => \$id_cut,
           "l=i" => \$length_cut,
           "b=i" => \$bit_cut,
	   "e=f" => \$e_cut,
	   "a=f" => \$aln_cut,
	   "n=i" => \$num,
	   "d=s" => \$tophit_database,

);

if (not defined $e_cut) { $e_cut = 10; }
if (not defined $num) { $num = 1; }
if (not defined $tophit_database) { print "\n Please specify database for tophits\n"; exit; }

foreach my $input (@ARGV) {

	open (BLAST, "$input");

	my $m8_name = $input;
	$m8_name =~ s/.*\///g;
	open (RESULTS, ">$input\_single_tophit.txt");
#	open (NOTES, ">$input\_report.txt");
	
	my $total = 0;
	my $p_count = 0;
	my $assigned = 0;
	my $found = 0;
	my %hitcount;
	my %lookup_hits;
	my %lookup_matches;
	my $last;
	my %matches;
	my @best_match;
	my $num_matches = 0;

	while (my $blast = <BLAST>) {
		
		chomp $blast;
		my @result = split (/\t/, $blast);
		
		if ($result[0] ne $last) {  # New hit!  Process last query.
			if ( $num_matches == 1 ) {	# One best hit - assign and add to hit file.
				my $hit_name = $best_match[1];
				$assigned++;
				$found++;
 				$hitcount{$hit_name} = $hitcount{$hit_name} + 1;
				my $temp = join ("\t", @best_match);
				print RESULTS "$temp\n";
	
			} elsif ($num_matches > 1) {	# Multiple best hits - wait to make assignment.
				$found++;
				foreach my $hit ( keys %matches) {
					if (exists $lookup_hits{$hit}) {
						$lookup_hits{$hit}{'found'} = $lookup_hits{$hit}{'found'} + 1;
					} else {
						$lookup_hits{$hit}{'found'} = 1;
						$lookup_hits{$hit}{'hitcount'} = 0;
					}
					$lookup_matches{$last}{$hit} = $matches{$hit};	
				}
			}
			$num_matches = 0;	#Reset counters and assignment
			$last = $result[0];
			$total++;
		}

		next unless ($result[9] > $bit_cut && $result[4] > $length_cut && $result[3] > $id_cut && $result[10] < $e_cut && $result[11] > $aln_cut);
		my $query = $result[0];
		my $bit = $result[9];
		my $sign = $result[10];
		my $hit_name = $result[1];

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
 		$hitcount{$hit_name} = $hitcount{$hit_name} + 1;
		my $temp = join ("\t", @best_match);
		print RESULTS "$temp\n";
	} elsif ($num_matches > 1) {	# Multiple best hits - wait to make assignment.
		$found++;
		foreach my $hit ( keys %matches) {
			if (exists $lookup_hits{$hit}) {
				$lookup_hits{$hit}{'found'} = $lookup_hits{$hit}{'found'} + 1;
			} else {
				$lookup_hits{$hit}{'found'} = 1;
				$lookup_hits{$hit}{'hitcount'} = 0;
			}
			$lookup_matches{$last}{$hit} = $matches{$hit};	
		}
	}

	print "\nFile $m8_name loaded.  Parsed $total total queries\n";
	print "First pass:  $found hits passed cutoffs, $assigned hits assigned\n";

	print "Parsing hitcounts file, creating temporary hitcounts file.\n";

	my $total = 0;
	my $p_count = 0;
	my $updated = 0;

	open (HITCOUNTS, "$hitcounts_path\/$tophit_database\_Hitcounts_all");
	open (TEMP, ">$hitcounts_path\/Hitcounts_temp");

	while (my $line = <HITCOUNTS>) {

		$total++;
		$p_count++;
#		if ($p_count == 50000) {
#			print "\t$total hitcounts parsed, $updated updated\n";
#			$p_count = 0;
#		}
		chomp $line;
		my @pep = split (/\t/, $line);
		unless (exists $lookup_hits{$pep[0]} || exists $hitcount{$pep[0]}) {
			print TEMP "$pep[0]\t$pep[1]\n";
			next;
		}
		if (exists $hitcount{$pep[0]}) {
			$pep[1] = $pep[1] + $hitcount{$pep[0]};
			delete ( $hitcount{$pep[0]} );
			$updated++;
		}
		if (exists $lookup_hits{$pep[0]}) {
			$lookup_hits{$pep[0]}{'hitcount'} = $pep[1];
		}
		print TEMP "$pep[0]\t$pep[1]\n";
	}
	print "\t$total references found in hitcounts file, $updated updated, ";

	close HITCOUNTS;

#	print "Adding new hits: ";
	$total = 0;
	foreach my $hit ( keys %hitcount ) {
		if (exists $lookup_hits{$hit}) {
			$lookup_hits{$hit}{'hitcount'} = $hitcount{$hit};
		}
		print TEMP "$hit\t$hitcount{$hit}\n";
		$total++;
	}
	print "$total new genes added.\n";
	my %hitcount;

	close TEMP;

	my $to_find = (scalar keys %lookup_matches);
	my $count = 0;
	$p_count = 0;

	print "Second pass:  Assigning queries with multiple hits, $to_find to assign.\n";

	foreach my $query (sort { (scalar keys %{$lookup_matches{$a}}) <=> (scalar keys %{$lookup_matches{$b}}) } keys %lookup_matches ) {
#		if ($p_count == 5000) {
#			print "\t$count searches parsed, $assigned total hits assigned\n";
#			$p_count = 0;
#		}
		$count++;
		$p_count++;
		my $num_matches = scalar keys %{$lookup_matches{$query}};
#		print NOTES "$query had $num_matches equal hits.\n";
		my ($best_hit_name, $best_count);
		foreach my $hit_name (sort { $lookup_hits{$b}{'found'} <=> $lookup_hits{$a}{'found'} } keys %{$lookup_matches{$query}}) {
			if (not defined $best_hit_name) {
				$best_hit_name = $hit_name;
				$best_count = $lookup_hits{$hit_name}{'hitcount'};
			} elsif ($lookup_hits{$hit_name}{'hitcount'} > $best_count) {
				$best_hit_name = $hit_name;
				$best_count = $lookup_hits{$hit_name}{'hitcount'};
			}
			my $temp = join ("\t", @{$lookup_matches{$query}{$hit_name}});
#			print NOTES "$lookup_hits{$hit_name}{'hitcount'} hits\t$temp\n";
		}
#		print NOTES "Using $best_hit_name\n\n";
		my $temp = join ("\t", @{$lookup_matches{$query}{$best_hit_name}});
		print RESULTS "$temp\n";
		$assigned++;
		if (exists $hitcount{$best_hit_name}) {
			$hitcount{$best_hit_name} = $hitcount{$best_hit_name} + 1;
		} else {
			$hitcount{$best_hit_name} = 1;
		}
		$lookup_hits{$best_hit_name}{'hitcount'} = $lookup_hits{$best_hit_name}{'hitcount'} + 1;
	}
	print "Second pass complete.  $count searches parsed, $assigned total hits assigned\n";
#	close NOTES;
	close RESUTS;

	print "Updating hitcounts file.\n";

	my $total = 0;
	my $p_count = 0;
	my $updated = 0;

	open (HITCOUNTS, ">$hitcounts_path\/$tophit_database\_Hitcounts_all");
	open (TEMP, "$hitcounts_path\/Hitcounts_temp");

	while (my $line = <TEMP>) {
		$total++;
		$p_count++;
#		if ($p_count == 50000) {
#			print "\t$total hitcounts parsed, $updated updated\n";
#			$p_count = 0;
#		}
		chomp $line;
		my @pep = split (/\t/, $line);
		if (exists $hitcount{$pep[0]}) {
			$pep[1] = $pep[1] + $hitcount{$pep[0]};
			$updated++;
			delete ( $hitcount{$pep[0]} );
		}
		print HITCOUNTS "$pep[0]\t$pep[1]\n";
	}
	print "\t$updated out of $total references updated, ";
	close TEMP;

#	print "\nAdding new hits: ";
	$total = 0;
	foreach my $hit ( keys %hitcount ) {
		print HITCOUNTS "$hit\t$hitcount{$hit}\n";
		$total++;
	}
	print "$total new genes added.\n";
	close HITCOUNTS;

}

exit;





