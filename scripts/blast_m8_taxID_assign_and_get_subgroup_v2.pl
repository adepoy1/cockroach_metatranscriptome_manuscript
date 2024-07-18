#Takes an NCBI taxonomy ID, identifies all subnodes, and extracts m8 files assigned to that subgroup.
#Default is greedy matching (all reads with at least one hit to that group).  -g N option requires all hits belong to that group.
#Uses m8L alltophits files (files in Liz's m8 format, including the top hit for each read and all equal hits).

#my $NCBI_Taxonomy_filepath = "/Users/otte/Bioinformatics/NCBI_Taxonomy/taxdump_05262010";
my $NCBI_Taxonomy_filepath = "/Volumes/G-DRIVE\ Thunderbolt\ 3/taxdump_SAGs";
#my $NCBI_Taxonomy_filepath = "/Users/otte/NCBI_Taxonomy/taxdump";

use strict;
use Getopt::Long;
# use Bio::Tree::Tree;
# use Bio::Tree::Node;

my ($e_cut, $bit_cut, $suffix, $target, $greedy, $outer_limit);

GetOptions(
	   "e=f" => \$e_cut,
	   "b=f" => \$bit_cut,
	   "s=s" => \$suffix,
	   "t=s" => \$target,
	   "g=s" => \$greedy,
	   "o=s" => \$outer_limit
);

if (not defined $e_cut) { $e_cut = 10; }
if (not defined $suffix) { print "Please enter suffix for output file\n"; exit; }
if (not defined $target) { print "Please enter Taxonomy ID of target\n"; exit; }
if (not defined $greedy) { $greedy = "Y"; }
if (not defined $outer_limit) { 
	$outer_limit = $target; 
} else {
	$greedy = "O";
}

#open (INPUT, "$NCBI_Taxonomy_filepath\nodes.dmp") || die "Could not open $NCBI_Taxonomy_filepath\nodes.dmp file\n";
open (INPUT, "/Volumes/G-DRIVE\ Thunderbolt\ 3/taxdump_SAGs/nodes.dmp") || die "Could not open nodes.dmp file\n";
my %nodes;

while (my $line = <INPUT>) {
	chomp $line;
	my @result = split (/\t\|\t/, $line);
	next if ($result[11] == 1);
	$nodes{$result[1]}{$result[0]} = 0;
}

$nodes{'285892'}{'939390'} = 0; #Add Roseobacter SAG AAA300J04 to unclassified Rhodobacterales
$nodes{'285892'}{'939374'} = 0; #Add Roseobacter SAG AAA298K06 to unclassified Rhodobacterales
$nodes{'285892'}{'939301'} = 0; #Add Roseobacter SAG AAA015O19 to unclassified Rhodobacterales

close INPUT;

# if (not defined $nodes{$target}) { "Target $target does not represent a group, using as \n"; }
if ($greedy eq "O" && not defined $nodes{$outer_limit}) { die "Outer limit $outer_limit not found in nodes.dmp file\n"; }

print "Reading NCBI Taxonomy Names File\n";

#open (INPUT, "$NCBI_Taxonomy_filepath\/names.dmp") || die "Could not open names.dmp file\n";
open (INPUT, "/Volumes/G-DRIVE\ Thunderbolt\ 3/taxdump_SAGs/names_wSAGs.dmp") || die "Could not open names.dmp file\n";

my %names;
my $target_name;
my $subgroup_name;

while (my $line = <INPUT>) {
	chomp $line;
	my @result = split (/\t\|\t/, $line);
	$names{$result[1]} = $result[0];
	next unless ($result[3] =~ m/scientific name/);
	if ($result[0] eq $target && not defined $target_name ) {
		$target_name = $result[1];
	}
	if ($result[0] eq $outer_limit && not defined $subgroup_name ) {
		$subgroup_name = $result[1];
	}
}
$names{" SAR324 cluster bacterium SCGC AAA240-J09 "} = "913331";
$names{"Phaeobacter gallaeciensis BS107"} = "391619";

close INPUT;

#open (INPUT, "/Users/otte/Bioinformatics/NCBI_taxonomy/acc.to.taxid.protein.additions") || die "Could not open names.dmp file\n";
#open (INPUT, "/Users/otte/NCBI_taxonomy/acc.to.taxid.protein.additions") || die "Could not open names.dmp file\n";

my %acc_to_taxID;

#while (my $line = <INPUT>) {
#	chomp $line;
#	my @result = split (/\t/, $line);
#	$acc_to_taxID{$result[0]} = $result[1];
#}
#close INPUT;


print "Building taxonomy subtree\n";

my %subtree;
$subtree{$outer_limit} = 1;
my @to_add = ( $outer_limit );
for (my $i = 0; exists $to_add[$i]; $i++) {
	my $current_id = $to_add[$i];
	foreach my $child_id (keys %{$nodes{$current_id}}) {
		next if (exists $subtree{$child_id});
		push (@to_add, $child_id);
		$subtree{$child_id} = 1;
	}
}

my %target_ids;
if ($target =~ m/ /) {
	@to_add = split(/ /, $target);
	foreach my $id (@to_add) {
		$target_ids{$id} = 1;
		print "Target species $names{$id} ($id\) found\n";
	}
} else {
	$target_ids{$target} = 1;
	@to_add = ( $target );
}
for (my $i = 0; exists $to_add[$i]; $i++) {
	my $current_id = $to_add[$i];
	foreach my $child_id (keys %{$nodes{$current_id}}) {
		next if (exists $target_ids{$child_id});
		push (@to_add, $child_id);
		$target_ids{$child_id} = 1;
	}
}


my $target_species = scalar keys %target_ids;
print "Targeting $target_species taxonomy IDs in subgroup $target_name\n";

print "\nParsing Files.";
if ($greedy =~ m/Y|y/) {
	print "Including all reads with at least one hit to targeted group\n";
	$greedy = 1;
} elsif ($greedy eq "O") {
	my $temp = scalar keys %subtree;
	print "Including reads within $subgroup_name subtree.  $temp acceptable nodes.\n";
} else {
	print "Including reads with all identifiable hits matching targeted group\n";
}

my %hit_assign;
my $matched_queries;
my $matched_hits;


foreach my $input (@ARGV) {

	my $outname = $input;
	$outname =~ s/.*\///g;
	
	open (OUT, ">$outname\_$suffix");

	%hit_assign = ();
	$matched_queries = 0;
	$matched_hits = 0;
	
	my $total = 0;
	my $taxID = 0;
	my $unassigned = 0;
	my $last;

	print "\nFile: $outname\n";
	
	open (BLAST, $input);
	while (my $line = <BLAST>) {

		chomp $line;
		$total++;
		my @match = split (/\t/, $line);
		next unless ($match[9] >= $bit_cut && $match[10] <= $e_cut);
		if ($match[0] ne $last) {
			$taxID++;
			&process_last_hit($last);
			$last = $match[0];
		}
		my %id;
		foreach my $species ($match[2] =~ /\[(.*?)\]/g ) {
			if (exists $names{$species}) {
				my $tax = $names{$species};
				$id{$tax} = $species;
			}
		}
		if (scalar keys %id < 1) {
			foreach my $acc ( $match[1] =~ /\|(.*?)\|/g, $match[1] =~ /.*\|(.*?)$/ ) {
				$acc =~ s/\.\d+//;
				if (exists $acc_to_taxID{$acc}) {
					my $tax = $acc_to_taxID{$acc};
					$id{$tax} = $acc;
				}
			}
		}		
		if (scalar keys %id >= 1) {
			foreach my $tax (keys %id) {
				$hit_assign{$match[0]}{$tax}{$line} = 1;
			}
		}	
		if (scalar keys %id < 1) {
			$unassigned++;
			if ($unassigned < 10) {
#				print "$line\n";
			}
		}
	}
	&process_last_hit($last);
	close OUT;
	close BLAST;
	print "$taxID queries found\t$total hits found, $unassigned did not have recognizable species names\n";
	print "$matched_queries queries matched to subgroup, $matched_hits hits printed\n";
	
}

sub process_last_hit {
	my $last = shift @_;
	my %recognized;
	my %nonredundant_hits;
	my $all_hits;
	my $target_hits;
	foreach my $hit (keys %{$hit_assign{$last}}) {
		foreach my $blast_hit ( keys %{$hit_assign{$last}{$hit}}) {
			$all_hits++;
			$nonredundant_hits{$blast_hit} = 1;
			if (exists $subtree{$hit}) {
				$target_hits++;
				if (exists $target_ids{$hit}) {
					$recognized{$blast_hit} = 1;
				}
			}
		}
	}
	if ( scalar  keys %recognized >= 1 && ( $greedy == 1 || $target_hits == $all_hits ) ) {
		foreach my $hit (keys %recognized) {
#		foreach my $hit (keys %nonredundant_hits) {
#			next if (exists $recognized{$hit});
			print OUT "$hit\n";
			$matched_hits++;
		}
		$matched_queries++;
	}
	%hit_assign = ();
}

exit;
