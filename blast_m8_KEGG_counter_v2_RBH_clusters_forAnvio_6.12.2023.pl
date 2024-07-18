#Uses KeggGenes tophit output to count hits to each level of the KEGG hierarchy
# Please change the path for the KEGG hierarchy file (ko00001.keg)
#   > find it here: https://www.genome.jp/kegg-bin/download_htext?
# Usage:
# perl blast_m8_KEGG_counter_v2_RBH_clusters_forAnvio.pl [single_tophit.txt files] -o [output prefix] -a [gene clusters summary]

use strict;
use Getopt::Long; 


my ($output_name, $assignment_file);

GetOptions(
	"o=s" => \$output_name,
	"a=s" => \$assignment_file
);

if (not defined $output_name) { $output_name = 'temp'; }      

if (not defined $assignment_file) { print "Please enter anvio clusters SUMMARY file with -a\n"; exit; }  

my %ko_hier;
my %ko_info;
my $count = 0;

#open (KEGG, "/Users/otte/Bioinformatics/KEGG_files/ko00001.keg") or die "Cannot open the KEGG hierarchy file!\n";    

open (KEGG, "/Volumes/G-DRIVE\ Thunderbolt\ 3/ko00001.keg") or die "Cannot open the KEGG hierarchy file! Change path in script\n";

my ($first, $second, $third, $fourth);
while (my $line = <KEGG>) {
        chomp $line;
        if ($line =~ /^A/) {
            ($first) = $line =~ /^A(.*)/;
        }
        if ($line =~ /^B\s+\w+/) {
            ($second) = $line =~ /^B\s+(\w+.*)/;
        }
        if ($line =~ /^C\s+\w+/) {
            ($third) = $line =~ /^C\s+(\w+.*)/;
        }
        if ($line =~ /^D\s+/) {
                next if ($first =~ m/Human Diseases/);
                next if ($first =~ m/Organismal Systems/);
#               next if ($third =~ m/\[BR/);
#                my ($ko) = $line =~ /^D\s+(K\d+)/;
                my ($ko, $desc) = $line =~ /^D\s+(K\d+)\s+(.*)/;
                my $heir = "$first\t$second\t$third\t$ko $desc";
                if ($count <5) {print "$heir\n";}
                $ko_hier{$ko}{'first'}{$first} = 1;
                $ko_hier{$ko}{'second'}{$second} = 1;
                $ko_hier{$ko}{'third'}{$third} = 1;
                $ko_hier{$ko}{'desc'} = "$ko $desc";
#               $path_hier{$third} = "$first\t$second";
                $count++;
                
        }
}
close KEGG;

foreach my $ko (keys %ko_hier) {
	my $desc = $ko_hier{$ko}{'desc'};
	my @full_hier;
	
	foreach my $first (keys %{$ko_hier{$ko}{'first'}}) {
		$full_hier[0] = $full_hier[0] . "$first ";
	}

	foreach my $second (keys %{$ko_hier{$ko}{'second'}}) {
		$full_hier[1] = $full_hier[1] . "$second ";
	}  
    	
	foreach my $third (keys %{$ko_hier{$ko}{'third'}}) {
		$full_hier[2] = $full_hier[2] . "$third ";
	}
      	$full_hier[3] = "$ko $desc";
	$ko_hier{$ko}{'desc'} = join ("\t", @full_hier);
}


my %RBH_assign;
my %ko_assign;
my %RBH_info;

my $count;

open (RBH, "$assignment_file") || die "could not open anvio clusters SUMMARY file $assignment_file\n";


my $header = <RBH>;
chomp $header;
my @header = (split /\t/, $header);

my $acc_col;
my $desc_col;
my $kofamacc_col;
my $kofam_col;
my $KEGGclass_acc_col;
my $KEGGclass_col;
my $KEGGmod_acc_col;
my $KEGGmod_col;
my $COGpath_acc_col;
my $COGpath_col;
my $COGcat_acc_col;
my $COGcat_col;
my $COGfunct_acc_col;
my $COGfunct_col;
my $CAZyme_col;
my $KEGG_BRITE_ACC;
my $KEGG_BRITE;
my $col = 12;

while ($col <= 29) {
	if ($header[$col] =~ m/NCBI_PGAP_ACC/) {
		$acc_col = $col;
	} 
	if ($header[$col] =~ m/NCBI_PGAP/) {
		$desc_col = $col;
	} 
	if ($header[$col] =~ m/KOfam_ACC/) {
		$kofamacc_col = $col;
	}  
	if ($header[$col] =~ m/KOfam/) {
		$kofam_col = $col;
	} 
	if ($header[$col] =~ m/KEGG_Class_ACC/) {
		$KEGGclass_acc_col = $col;
	} 
	if ($header[$col] =~ m/KEGG_Class/) {
		$KEGGclass_col = $col;
	} 
	if ($header[$col] =~ m/KEGG_Module_ACC/) {
		$KEGGmod_acc_col = $col;
	} 
	if ($header[$col] =~ m/KEGG_Module/) {
		$KEGGmod_col = $col;
	} 
	if ($header[$col] =~ m/COG20_PATHWAY_ACC/) {
		$COGpath_acc_col = $col;
	} 
	if ($header[$col] =~ m/COG20_PATHWAY/) {
		$COGpath_col = $col;
	} 
	if ($header[$col] =~ m/COG20_CATEGORY_ACC/) {
		$COGcat_acc_col = $col;
	} 
	if ($header[$col] =~ m/COG20_CATEGORY/) {
		$COGcat_col = $col;
	} 
	if ($header[$col] =~ m/COG20_FUNCTION_ACC/) {
		$COGfunct_acc_col = $col;
	} 
	if ($header[$col] =~ m/COG20_FUNCTION/) {
		$COGfunct_col = $col;
	} 
	if ($header[$col] =~ m/CAZyme/) {
		$CAZyme_col = $col;
	}
	if ($header[$col] =~ m/KEGG_BRITE_ACC/) {
		$KEGG_BRITE_ACC = $col;
	}
	if ($header[$col] =~ m/KEGG_BRITE/) {
		$KEGG_BRITE = $col;
	}
	$col++;
}

if (not defined $acc_col) {
	print "WARNING: ACC column not found\n";
}
if (not defined $desc_col) {
	print "WARNING: desc column not found\n";
}
if (not defined $kofamacc_col) {
	print "WARNING: KOfam ACC column not found\n";
}
if (not defined $kofam_col) {
	print "WARNING: KOfam column not found\n";
}
if (not defined $KEGGclass_acc_col) {
	print "WARNING: KEGG class ACC column not found\n";
}
if (not defined $KEGGclass_col) {
	print "WARNING: KEGG class column not found\n";
}
if (not defined $KEGGmod_acc_col) {
	print "WARNING: KEGG module ACC column not found\n";
}
if (not defined $KEGGmod_col) {
	print "WARNING: KEGG module column not found\n";
}
if (not defined $COGpath_acc_col) {
	print "WARNING: COG20 pathway ACC column not found\n";
}
if (not defined $COGpath_col) {
	print "WARNING: COG20 pathway column not found\n";
}
if (not defined $COGcat_acc_col) {
	print "WARNING: COG20 category ACC column not found\n";
}
if (not defined $COGcat_col) {
	print "WARNING: COG20 category column not found\n";
}
if (not defined $COGfunct_acc_col) {
	print "WARNING: COG20 function ACC column not found\n";
}
if (not defined $COGfunct_col) {
	print "WARNING: COG20 function column not found\n";
}
if (not defined $CAZyme_col) {
	print "WARNING: CAZyme column not found\n";
}
if (not defined $KEGG_BRITE_ACC) {
	print "WARNING: KEGG BRITE ACC column not found\n";
}
if (not defined $KEGG_BRITE) {
	print "WARNING: KEGG BRITE column not found\n";
}


while (my $line = <RBH>) {
	chomp $line;
	#next unless ($line =~ m/^>/);

	my @result = split (/\t/, $line);
	my $cluster = $result[1];
	my $acc = $result[$acc_col];
	$acc =~ s/\.\d+//;

if (exists $RBH_info{$cluster}{'acc'}) {
		$RBH_info{$cluster}{'acc'} = $RBH_info{$cluster}{'acc'} . " $acc";
		$RBH_info{$cluster}{'species'}{$result[3]} = 1;
		$RBH_info{$cluster}{'desc'}{$result[$desc_col]} = 1;
		$RBH_info{$cluster}{'KOfam_ACC'}{$result[$kofamacc_col]} = 1;
		$RBH_info{$cluster}{'KOfam'}{$result[$kofam_col]} = 1;
		$RBH_info{$cluster}{'KEGG_Class_ACC'}{$result[$KEGGclass_acc_col]} = 1;
		$RBH_info{$cluster}{'KEGG_Class'}{$result[$KEGGclass_col]} = 1;
		$RBH_info{$cluster}{'KEGG_Module_ACC'}{$result[$KEGGmod_acc_col]} = 1;
		$RBH_info{$cluster}{'KEGG_Module'}{$result[$KEGGmod_col]} = 1;
		$RBH_info{$cluster}{'COG20_pathway_ACC'}{$result[$COGpath_acc_col]} = 1;
		$RBH_info{$cluster}{'COG20_pathway'}{$result[$COGpath_col]} = 1;
		$RBH_info{$cluster}{'COG20_category_ACC'}{$result[$COGcat_acc_col]} = 1;
		$RBH_info{$cluster}{'COG20_category'}{$result[$COGcat_col]} = 1;
		$RBH_info{$cluster}{'COG20_function_ACC'}{$result[$COGfunct_acc_col]} = 1;
		$RBH_info{$cluster}{'COG20_function'}{$result[$COGfunct_col]} = 1;
		$RBH_info{$cluster}{'CAZyme'}{$result[$CAZyme_col]} = 1;
		$RBH_info{$cluster}{'KEGG_BRITE_ACC'}{$result[$KEGG_BRITE_ACC]} = 1;
		$RBH_info{$cluster}{'KEGG_BRITE'}{$result[$KEGG_BRITE]} = 1;

	} else {
		$RBH_info{$cluster}{'acc'} = $acc;
		$RBH_info{$cluster}{'species'}{$result[3]} = 1;
		$RBH_info{$cluster}{'desc'}{$result[$desc_col]} = 1;
		$RBH_info{$cluster}{'KOfam_ACC'}{$result[$kofamacc_col]} = 1;
		$RBH_info{$cluster}{'KOfam'}{$result[$kofam_col]} = 1;
		$RBH_info{$cluster}{'KEGG_Class_ACC'}{$result[$KEGGclass_acc_col]} = 1;
		$RBH_info{$cluster}{'KEGG_Class'}{$result[$KEGGclass_col]} = 1;
		$RBH_info{$cluster}{'KEGG_Module_ACC'}{$result[$KEGGmod_acc_col]} = 1;
		$RBH_info{$cluster}{'KEGG_Module'}{$result[$KEGGmod_col]} = 1;
		$RBH_info{$cluster}{'COG20_pathway_ACC'}{$result[$COGpath_acc_col]} = 1;
		$RBH_info{$cluster}{'COG20_pathway'}{$result[$COGpath_col]} = 1;
		$RBH_info{$cluster}{'COG20_category_ACC'}{$result[$COGcat_acc_col]} = 1;
		$RBH_info{$cluster}{'COG20_category'}{$result[$COGcat_col]} = 1;
		$RBH_info{$cluster}{'COG20_function_ACC'}{$result[$COGfunct_acc_col]} = 1;
		$RBH_info{$cluster}{'COG20_function'}{$result[$COGfunct_col]} = 1;
		$RBH_info{$cluster}{'CAZyme'}{$result[$CAZyme_col]} = 1;
		$RBH_info{$cluster}{'KEGG_BRITE_ACC'}{$result[$KEGG_BRITE_ACC]} = 1;
		$RBH_info{$cluster}{'KEGG_BRITE'}{$result[$KEGG_BRITE]} = 1;

	}
	$RBH_assign{$acc} = $cluster;  
	if ($result[$kofamacc_col] =~ m/\w+/){
		$ko_assign{$cluster}{$result[$kofamacc_col]} = 1;
	}

#	$count++;
#	if ($count < 10) { print "\t$acc\t$cluster\n"; }    

#	my ($cluster) = $result[0] =~ m/>(.*)/;
#	$RBH_info{$cluster}{'acc'} = $result[1];
#	$RBH_info{$cluster}{'species'} = $result[3];
#	$RBH_info{$cluster}{'desc'} = $result[$desc_col];
#	$RBH_info{$cluster}{'ko'} = $result[$kofamacc_col];
#	foreach my $ko ( $result[$kofamacc_col] =~ /(K\d{5})/g ) {
#		$ko_assign{$cluster}{$ko} = 1;
#	} 
#	my @genes = split (/ /, $result[1]);
#	foreach my $acc (@genes) {
#		$acc =~ s/.*\|(.*?)\|$/$1/;
#		$acc =~ s/.*\|(.+)$/$1/;
#		$acc =~ s/\.\d+$//;
#		$RBH_assign{$acc} = $cluster;
#	}      
}            

my %counter_all;
my %counter_db;
my $db_num = 0;
my %db_size;
my @databases;
my %hit_info;
my %unassigned;

foreach my $input (@ARGV) {
	my $db = $input;
	$db =~ s/_sub_cDNA_norRNAs.vs.nr_May312010.bv50_m8_all_NoDup_alltophits_b50_Synechococcus_by_taxID_cyanos/cDNA/;
	$db =~ s/\.VS.*//;
	$db =~ s/\.vs.*//;
	$db =~ s/.*\///g;
	
	if ($input =~ m/NoDup/) {
	$db = $db . "_NoDuplicates";
	}

	$db_num++;

#	print "\nFile $db_num\: $db\n";

	my $total = 0;
	my $p_count = 0;

	
	open (BLAST, $input);
	while (my $line = <BLAST>) {
#		last if ($total == 500);
		chomp $line;
		next if ($line =~ m/ZP_01468898.1/);
		my @match = split (/\t/, $line);
		my $query = $match[0];
		my $acc = $match[1];
		my $desc = $match[2];

		$acc =~ s/.*ref\|(.*?)\.\d+\|.*/$1/;
		$acc =~ s/.*gb\|(.*?)\.\d+\|.*/$1/;
		$acc =~ s/.*\|(.*?)\|$/$1/;
		$acc =~ s/.*\|(.+)/$1/;
		$acc =~ s/\.\d+$//;
		if (exists $RBH_assign{$acc}) {
			$acc = $RBH_assign{$acc};
		} else {
			my $found = 0;
			foreach my $temp ( $desc =~ /\|(.*?)\.\d+\|/g ) {

			if (exists $RBH_assign{$temp}) {
				$acc = $RBH_assign{$temp};
				$found++;
			}
		}
		if ( $found == 0 ) {
			$unassigned{$acc} = $match[2];
		}
	}

	next unless ($match[9] >= 50);
	$total++;

	$counter_all{'hit'}{$acc}++;
	$counter_db{$db}{'hit'}{$acc}++;

}      	

close BLAST;      	

print "$total\t$db\n";
push @databases, $db;
$db_size{$db} = $total;

}

my $num_pep = scalar keys %{$counter_all{'hit'}};

print "\nBlast results parsed.  $num_pep KEGG Genes identified\n";

if (scalar keys %unassigned >= 1) {
	my $not_found = scalar keys %unassigned;
	print "$not_found references not located in clusters file:\n";
	foreach my $acc (keys %unassigned) {
		print "\t$acc $unassigned{$acc}\n";
	} 
}

my %ko_counter;
my %ko_info;

foreach my $acc (keys %{$counter_all{'hit'}}) {
	my %hier;
	my $assigned;
	if (exists $ko_assign{$acc} && $ko_assign{$acc} =~ m/\w+/) {
		foreach my $ko ( keys %{$ko_assign{$acc}} ) {
			$ko_counter{$ko}{$acc} = 1;
			$hier{'ko'}{$ko} = 1;
			if (exists $ko_hier{$ko}) {
				foreach my $h1 ( keys %{$ko_hier{$ko}{'first'}} ) {
					$hier{'first'}{$h1} = 1;
					$assigned++;
				}
				foreach my $h2 ( keys %{$ko_hier{$ko}{'second'}} ) {
					$hier{'second'}{$h2} = 1;
				} 
				foreach my $h3 ( keys %{$ko_hier{$ko}{'third'}} ) {
				$hier{'third'}{$h3} = 1;
				}
			} 
		}
	} else { 
		$hier{'ko'}{'unassigned'} = 1;
		$ko_counter{'unassigned'}{$acc} = 1;
	}
	unless ($assigned >= 1) {
		$hier{'first'}{'unassigned'} = 1;
		$hier{'second'}{'unassigned'} = 1;
		$hier{'third'}{'unassigned'} = 1;
	}
	foreach my $level ( keys %hier) {
		foreach my $kegg ( keys %{$hier{$level}} ) {
			$counter_all{$level}{$kegg} += $counter_all{'hit'}{$acc};
			foreach my $db ( keys %counter_db) {
				if (exists $counter_db{$db}{'hit'}{$acc}) {
					$counter_db{$db}{$level}{$kegg} += $counter_db{$db}{'hit'}{$acc};
				}
			}
			$ko_info{$acc}{$level} = $ko_info{$acc}{$level} . "$kegg ";
		}
	}
}

print "\nProcessing Peptide Hitcounts\n";


my %fact;

foreach my $level (keys %counter_all) {
	my %count;
	my $total = 0;
	my $p_count = 0;
	#next unless ($level eq 'ko');
	print "\nProcessing at $level Hierarchy Level";

	open (OUT, ">$output_name\_Kegg_$level\_hitcounts.txt");

	if ($level eq 'hit') {
		print OUT "Cluster\tAcc\tDesc\tSpecies\tKOfam\tKOfam_ACC\tKEGG_Class\tKEGG_Module\tCOG20_pathway\tCOG20_category\tCOG20_function\tCOG20_function_ACC\tCAZyme\tKEGG_BRITE_ACC\tKEGG_BRITE\t";
	} elsif ($level eq 'ko') {
		print OUT "ko\tNum Assigned\tAssigned clusters\tLevel1\tLevel2\tLevel3\tDesc\t";
	} else {
		print OUT "Hierarchy\tPathway";
	} 
	
	foreach my $db ( @databases ) {
		print OUT "\t$db Hits";
	} 
     	
	print OUT "\t";

	foreach my $db ( @databases ) {
		print OUT "\t$db Percent";
	}      	print OUT "\t";
	
	foreach my $db ( @databases ) {
		print OUT "\t.5 $db Percent";
	}
	print OUT "\t";

	foreach my $hit (sort { $counter_all{$level}{$b} <=> $counter_all{$level}{$a} } keys %{$counter_all{$level}}) {


		if ($level eq 'hit') {
			if (exists $unassigned{$hit}) {
				print OUT "\n$hit\t$unassigned{$hit}\t\t\t\t\t\t\t\t\t\t\t\t\t\t";
			} else {
				print OUT "\n$hit\t$RBH_info{$hit}{'acc'}\t";

				my @cat = ('desc', 'species', 'KOfam','KOfam_ACC', 'KEGG_Class', 'KEGG_Module', 'COG20_pathway', 'COG20_category', 'COG20_function', 'COG20_function_ACC', 'CAZyme', 'KEGG_BRITE_ACC', 'KEGG_BRITE');
				foreach my $to_print (@cat) {
					my $print_col;
					foreach my $entry (keys %{$RBH_info{$hit}{$to_print}}) {
						if ($print_col =~ m/\w+/ && $entry =~ m/\w+/) {
							$print_col = $print_col . " $entry";
						} elsif ($entry =~ m/\w+/) {
							$print_col = $entry;
						}
					}
					print OUT "$print_col\t";	
				}
			}
		} elsif ($level eq 'ko') {
			my $count = scalar keys %{$ko_counter{$hit}};
			my $acc;
			foreach my $temp (sort {$a cmp $b} keys %{$ko_counter{$hit}}) {
				$acc = $acc . "$temp ";
			}
			print OUT "\n$hit\t$count\t$acc\t";

			if (exists $ko_hier{$hit}{'desc'} ) {
			
				print OUT "$ko_hier{$hit}{'desc'}\t";

			} else {
				print OUT "\t\t\t\t";
			}

		} else {

			print OUT "\n$hit\t";
			my $temp = $hit;
			$temp =~ s/ .PATH:.*//;
			$temp =~ s/\d+ //;
			print OUT "$temp";
		}

#		if ($total < 3) {
#			print "$hit\n";
#		}

#		if ($p_count == 10000) {
#			print "$total hits processed\n";
#			$p_count = 0;
#		}

		$total++;
		$p_count++;

		foreach my $db (@databases) {
			if (exists $counter_db{$db}{$level}{$hit}) {
				$count{$db} = $counter_db{$db}{$level}{$hit};
			} else { $count{$db} = 0; }
			print OUT "\t$count{$db}";
      		}

		print OUT "\t";


		foreach my $db (@databases) {
			my $temp = ($count{$db} / $db_size{$db});
			print OUT "\t$temp";
		}

		print OUT "\t";

		foreach my $db (@databases) {
			my $temp = ($count{$db} / $db_size{$db});
			if ($temp == 0) {
				$temp = (0.5 / $db_size{$db});
			}
			print OUT "\t$temp";
		}

	}

	close OUT;
	print " $total peptides processed\n";
}            