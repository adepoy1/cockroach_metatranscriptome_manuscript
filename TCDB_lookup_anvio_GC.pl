


#This program should take a TCDB superfamily file, and a hitcount file, and add in the family information for each hit.


#usage perl TCDB_lookup.pl superfamilies_file families_file input





use strict;


use Getopt::Long;





my %family_lookup;





open (FAMILY, "$ARGV[1]");


my $header = <FAMILY>;





while (my $line = <FAMILY>) {


	


	chomp $line;


	$line =~ s/\t\r//;


	my @result = split (/\t/, $line);


	$family_lookup{$result[0]} = $line;





}



my %superfamily_lookup;

open (SUPERFAMILY, "$ARGV[0]");


my $header = <SUPERFAMILY>;





while (my $line = <SUPERFAMILY>) {


	


	chomp $line;


	$line =~ s/\t*\r//;


	my @result = split (/\t/, $line);


	$superfamily_lookup{$result[0]} = $line;





}





close FAMILY;





open (INPUT, "$ARGV[2]");


open (OUTPUT, ">$ARGV[2]\_annotated.txt");





my $header = <INPUT>;


print OUTPUT "GC\tTCID\tSubfamily\tFamily\tFam_abbreviation\tSuperfamily\n";





while (my $line = <INPUT>) {


	


	chomp $line;


	my @result = split (/\t/, $line);
	
	$result[0] =~ /(GC_\d{8})/;

	my $GC = $1;


	my @title = split (/\|/, $result[1]);


	my $entry = $title[3];

	#print "$entry\n";

	$title[3] =~ /(\w+?\.\w+?\.\w+?)\./;

	my $family = $1;
	
	$title[3] =~ /(\w+?\.\w+?\.\w+?\.\w+?)\./;

	my $subfamily = $1;


	#print "$family\t$title[3]\n";


	if (exists $superfamily_lookup{$entry}) {


		print OUTPUT "$GC\t$superfamily_lookup{$entry}\n";

	} elsif (exists $family_lookup{$family}) {

		print OUTPUT "$GC\t$title[3]\t$subfamily\t$family_lookup{$family}\tnone\n";


	} else {


		print OUTPUT "$GC\tNA\tNA\tNA\tNA\n";


		print "Could not find $family or $entry in database\n";


	}





}





close INPUT;














