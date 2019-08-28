#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my %complete_data;

my $mirbase = "hairpin.fa";
my $dir_rfam = "Fasta_RFAM"; #Here is the Directory
my $mirgenedb = "ALL.gff";

my $existMirbase = existsFile($mirbase); 
my $existsMirgenedb = existsFile($mirgenedb);
my $DIR;
#Here open and process RFAM files
opendir($DIR, $dir_rfam) or die "The RFAM path is not available\n";
my @rfam_files = readdir($DIR);
close $DIR;

foreach my $rfam (@rfam_files){
	next if ($rfam =~ /^\.$|^\.\.$/);
	my $existsRfam = existsFile("$dir_rfam/$rfam");
	if ($existsRfam == 1){
		convert_all_data("$dir_rfam/$rfam", "rfam");
	}
}

if ($existMirbase == 1){
	convert_all_data($mirbase, "mirbase");
}


if ($existsMirgenedb == 1){
	convert_all_data($mirgenedb, "mirgenedb");
}

#print Dumper \%complete_data;

generate_table(\%complete_data);
count_data(\%complete_data);

#####Subs
sub existsFile {
	my $file = shift;
	my $reporter = 1;
	unless (-e $file && !-z $file){ #Checks existence and non-emptyness
		my $msj = "The file: $file has not been created!\n";
		print "$msj\n";
		$reporter = 0;
	}
	return $reporter;
}

sub convert_all_data {
	my ($in, $mode) = @_;
	if ($mode eq "mirbase"){
		open my $IN, "< $in" or die "Please provide input file\n";
		while(<$IN>){
		chomp;
		if ($_ =~ /^\>/){
			#>cel-lin-4 MI0000002 Caenorhabditis elegans lin-4 stem-loop
			my @temporal = split /\t|\s+/, $_;
			my $acc = $temporal[1];
			my $specie = "$temporal[2]_$temporal[3]";
			my $nameFamily = $temporal[4];
			my $database = "miRBase";
			push @{$complete_data{$acc}{"Specie"}}, $specie;
			push @{$complete_data{$acc}{"Name"}}, $nameFamily;
			push @{$complete_data{$acc}{"Database"}}, $database;		
			} else {
				next;
			}
		}
	} elsif ($mode eq "rfam"){
		open my $IN, "< $in" or die "Please provide input file\n";
		my $rnacentral = load_external_databases();
		while(<$IN>){
		chomp;
		if ($_ =~ /^\>/){
			#>CM000866.1/19874310-19874235 Callithrix jacchus chromosome 11, whole genome shotgun sequence
			my @split = split /\s+|\t/, $_;
			my $num = scalar @split;
			my $specie;
			if ($num > 3){
				$specie = "$split[1]_$split[2]";
			} else {
				$specie = "NA";
			}
			my $nameFamily = $in;
			$nameFamily =~ s/(\/.*\/|\.\.\/|.*\/)(.*)(\.fa|\.fasta)/$2/g;
			my $database = "RFAM";
			my $id = $split[0];
			$id =~ s/(\>)(.*)/$2/g;
			my $acc;
			if (exists $$rnacentral{$id}){
				$acc = $$rnacentral{$id};	
			} else {
				$acc = $nameFamily;
			}
			push @{$complete_data{$acc}{"Specie"}}, $specie;
			push @{$complete_data{$acc}{"Name"}}, $nameFamily;
			push @{$complete_data{$acc}{"Database"}}, $database;	
			}
		}
	} elsif ($mode eq "mirgenedb"){
		open my $IN, "< $in" or die "Please provide input file\n";
		my $tagsDB = load_species_names_mirgenedb();
		while(<$IN>){
		chomp;
		next if $_ =~ /^\#|^$/;
		my $type = (split /\s+|\t/, $_)[2];
		if ($type eq "pre_miRNA"){
			my $acc;
			my @split = split /\s+|\t/, $_;
			my $database = "mirgenedb";
			my $nameFamily = $split[-1];
			if ($nameFamily =~ m/\;Alias\=/){
				$acc = (split /\;/, $nameFamily)[1];
				$acc =~ s/(Alias\=)(.*)/$2/g;
			} else {
				$acc = $nameFamily;
				$acc =~ s/(ID\=)(.*)/$2/g;
			}				
			$nameFamily = (split /\;/, $nameFamily)[0];
			$nameFamily =~ s/(ID\=)(.*)/$2/g;
			my $tagSpecie = ( split /\-/,$nameFamily)[0]; 	
			my $specie;
		       	if (exists $$tagsDB{$tagSpecie}){
				$specie = $$tagsDB{$tagSpecie};
			} else {
				die "There is something wrong in the database assigment\n";
			}
			push @{$complete_data{$acc}{"Specie"}}, $specie;
			push @{$complete_data{$acc}{"Name"}}, $nameFamily;
			push @{$complete_data{$acc}{"Database"}}, $database;  	
			}		
		}	
	} else {
		die "$mode is not valid\n";
	}
	return;
}

sub load_external_databases {
	open my $DB, "< /scr/k70san/bogota_unal/miRNAsTunicata/miRNAture/Code/RFAM_Processing/Data/rnacentral_mirnas_mirbase_relation_RFAM_cleaned.txt" or die "Please load the file rnacentral_mirnas_mirbase_relation_RFAM_cleaned.txt\n";
	my %relation;
	while (<$DB>){
		#CM000321.4	135424451	135424388	f8ada79a032efa67ceea370ad2286a9e	URS00004E9304	URS00004E9304_9598	MI0002661
		chomp;
		my $id = (split /\s+|\t/, $_)[0];
		my $mirbaseID = (split /\s+|\t/, $_)[-1];
		my $start = (split /\s+|\t/, $_)[1];
		my $end = (split /\s+|\t/, $_)[2];
		$relation{"$id\/$start-$end"} = $mirbaseID;
	}
	return \%relation;
}

sub load_species_names_mirgenedb {
	open my $DB, "< test_mirgenedb_species_cleaned.txt" or die "Please obtain the file test_mirgenedb_species_cleaned.txt\n";
	my %db;
	while (<$DB>){
		chomp;
		$_ =~ s/^(.*)(\()(.*)(\))(.*)\s*$/$3 $5/g;
		my $tag = $5;
		my $name =  $3;
		$tag =~ s/\s+//g;
		$name =~ s/\s+/\_/g;
		$db{$tag} = $name;
	}
	return \%db;
}

sub generate_table {
	my $data = shift; #Hash of hash of arrays
	foreach my $acc (sort keys %{ $data }) {
		my $all_information = $$data{$acc};
		my $database_complete = $$data{$acc}{"Database"};
		my $species_complete = $$data{$acc}{"Specie"};
		my $names_complete = $$data{$acc}{"Name"};
		#Obtain index, based on database key:		
		my ($index_mirbase, $index_rfam, $index_mirgenedb); #Obtain the order of the databases in the data structure.
		$index_mirbase = get_index_array_match($database_complete, "miRBase");
		$index_rfam = get_index_array_match($database_complete, "RFAM");
		$index_mirgenedb = get_index_array_match($database_complete, "mirgenedb");
		#Database
		$database_complete = test_if_defined($database_complete, $index_mirbase);
		$database_complete = test_if_defined($database_complete, $index_rfam);
		$database_complete = test_if_defined($database_complete, $index_mirgenedb);
		#Species
		$species_complete = test_if_defined($species_complete, $index_mirbase);
		$species_complete = test_if_defined($species_complete, $index_rfam);
		$species_complete = test_if_defined($species_complete, $index_mirgenedb);
		#Names
		$names_complete = test_if_defined($names_complete,$index_mirbase);
		$names_complete = test_if_defined($names_complete,$index_rfam);
		$names_complete = test_if_defined($names_complete,$index_mirgenedb);
		##
		print_line_table($all_information, $acc);
	}	
	return;	
}

sub print_line_table {
	my ($complete, $acc) = @_;
	#Create spaces
	# mirbase, RFAM, mirgenedb
	my @databasesOrder = ("NA,") x 3;
	my @namesOrder = ("NA,") x 3;
	my @speciesOrder = ("NA,") x 3;
	#
	my ($databases, $names, $species);
	my ($index_mirbase, $index_rfam, $index_mirgenedb);
	for (my $i=0; $i <= 2; $i++){
		#Databases
		if ($$complete{"Database"}[$i]){
			if ($$complete{"Database"}[$i] eq "miRBase"){
				$index_mirbase = $i;		
				$databasesOrder[0] .= "$$complete{Database}[$i],";
				$namesOrder[0] .= "$$complete{Name}[$i],";
				$speciesOrder[0] .= "$$complete{Specie}[$i],";
			} elsif ("$$complete{Database}[$i]" eq "RFAM"){
				$index_rfam = $i;
				$databasesOrder[1] .= "$$complete{Database}[$i],";
				$namesOrder[1] .= "$$complete{Name}[$i],";
				$speciesOrder[1] .= "$$complete{Specie}[$i],";
			} elsif ("$$complete{Database}[$i]" eq "mirgenedb"){
				$index_mirgenedb = $i;
				$databasesOrder[2] .= "$$complete{Database}[$i],";
				$namesOrder[2] .= "$$complete{Name}[$i],";
				$speciesOrder[2] .= "$$complete{Specie}[$i],";
			}	
		}
	}
	my $final_databases = clean_array(\@databasesOrder);
	my $final_names = clean_array(\@namesOrder);
	my $final_order = clean_array(\@speciesOrder);
	#Name Databasemirbase DatabaseRFAM DatabaseMirgeneDB ... SpeciesMirgeneDB
	$final_databases = delete_redundant_fields($final_databases);
	$final_names = delete_redundant_fields($final_names);
	$final_order = delete_redundant_fields($final_order);
	$databases =  join "\t", @$final_databases;
	$names = join "\t", @$final_names;
	$species = join "\t", @$final_order;
	print "$acc\t$databases\t$names\t$species\n";
	return;
}

sub clean_array {
	my $array = shift;
	foreach my $ln (@$array){
		if ($ln =~ /^NA\,$/){
			$ln =~ s/^(NA)(,)$/$1/g;	
		} else {
			$ln =~ s/^(NA\,)(.*)(\,)$/$2/g;
		}
	}
	return $array;
}

sub delete_redundant_fields {
	my $array = shift;
	foreach my $ln (@$array){
		my @temp = split /\,/, $ln;
		my @unique = do { my %seen; grep { !$seen{$_}++ } @temp };
		$ln = join ",", @unique;
	}
	return $array;
}

sub get_index_array_match {
	my ($array, $pattern) = @_;
	my $index = "NA";
	my $length_database = (scalar @{$array} - 1);
	for (my $i = 0; $i <= $length_database; $i++ ){
		if ( $$array[$i] eq $pattern ){
			$index = $i;
			last;
		} else {
			;
		}
	}
	return $index;
}


sub test_if_defined {
	my ($array, $index) = @_;
	if ($index eq "NA"){
		push @{$array}, "NA";
	} else {
		;
	}
	return $array;
}

sub count_data {
	my $data = shift;
	my $c_rfam = 0;
	my $c_mirbase = 0;
	my $c_mirgenedb = 0;
	my $complete_families = 0;
	return;	
}

exit;
