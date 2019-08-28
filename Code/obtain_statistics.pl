#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my %complete_data;

my $mirbase = "hairpin.fa";
my $rfam = "RF00027.fa";
my $mirgenedb = "ALL.gff";

my $existMirbase = existsFile($mirbase); 
my $existsRfam = existsFile($rfam); 
my $existsMirgenedb = existsFile($mirgenedb);

if ($existMirbase == 1){
	convert_all_data($mirbase, "mirbase");
}

if ($existsRfam == 1){
	convert_all_data($rfam, "rfam");
}

if ($existsMirgenedb == 1){
	convert_all_data($mirgenedb, "mirgenedb");
}

print Dumper \%complete_data;

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
			$nameFamily =~ s/(\/.*\/|\.\.\/|)(.*)(\.fa|\.fasta)/$2/g;
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
		my $database_complete = $$data{$acc}{"Databases"};
		my $species_complete = $$data{$acc}{"Specie"};
		my $names_complete = $$data{$acc}{"Name"};
		#Obtain index, based on database key:		
		my ($index_mirbase, $index_rfam, $index_mirgenedb); #Obtain the order of the databases in the data structure.
		$index_mirbase = grep { $$database_complete[$_] eq 'miRBase' } (0 .. $$database_complete-1);
		$index_rfam = grep { $$database_complete[$_] eq 'RFAM' } (0 .. $$database_complete-1);
		$index_mirgenedb = grep { $$database_complete[$_] eq 'mirgenedb' } (0 .. $$database_complete-1);
		##	
		#Database
		$database_complete[$index_mirbase] = test_if_defined($database_complete[$index_mirbase]);
		$database_complete[$index_rfam] = test_if_defined($database_complete[$index_rfam]);
		$database_complete[$index_mirgenedb] = test_if_defined($database_complete[$index_mirgenedb]);
		#Species
		$species_complete[$index_mirbase] = test_if_defined($species_complete[$index_mirbase]);
		$species_complete[$index_rfam] = test_if_defined($species_complete[$index_rfam]);
		$species_complete[$index_mirgenedb] = test_if_defined($species_complete[$index_mirgenedb]);
		#Names
		$names_complete[$index_mirbase] = test_if_defined($names_complete[$index_mirbase]);
		$names_complete[$index_rfam] = test_if_defined($names_complete[$index_rfam]);
		$names_complete[$index_mirgenedb] = test_if_defined($names_complete[$index_mirgenedb]);
		##
	}	
	return;	
}

sub test_if_defined {
	my $variable = shift;
	my $return;
	if (length $variable > 0 && $variable){
		$return = $variable;
	} else {
		$return = "NA";
	}
	return $return;
}

sub count_data {
	my $data = shift;
	my $c_rfam = 0;
	my $c_mirbase = 0;
	my $c_mirgenedb = 0;
	my $complete_families = 0;
	
}


























exit;
