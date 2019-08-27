#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my %complete_data;

my $mirbase = "hairpin.fa";
my $rfam = "RF00027.fa";
my $mirgenedb = "ALL.gff";

$existMirbase = existsFile($mirbase); 
$existsRfam = existsFile($rfam); 
$existsMirgenedb = existsFile($mirgenedb);

if ($existMirbase == 1){
	convert_all_data($existMirbase, "mirbase");
}

if ($existsRfam == 1){
	convert_all_data($existsRfam, "rfam");
}

if ($existsMirgenedb == 1){
	convert_all_data($existsMirgenedb, "mirgenedb");
}

print Dumper \%complete_data;

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
		while($IN){
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
			}
		}
	} elsif ($mode eq "rfam"){
		open my $IN, "< $in" or die "Please provide input file\n";
		my $rnacentral = load_external_databases();
		while($IN){
		chomp;
		if ($_ =~ /^\>/){
			#>CM000866.1/19874310-19874235 Callithrix jacchus chromosome 11, whole genome shotgun sequence
			my @split = split /\s+|\t/, $_;
			my $specie = "$split[1]_$split[2]";
			my $nameFamily = $in;
			$nameFamily =~ s/(\/.*\/|\.\.\/)(.*)(\.fa|\.fasta)/$2/g;
			my $database = "RFAM";
			my $id = $split[0];
			$id =~ s/(\>)(.*)/$2/g;
			my $acc;
				if (exists $$rnacentral{$id}){
					$acc = $$rnacentral{$id};	
				} else {
					$acc = $nameFamily;
				}
			}
			push @{$complete_data{$acc}{"Specie"}}, $specie;
			push @{$complete_data{$acc}{"Name"}}, $nameFamily;
			push @{$complete_data{$acc}{"Database"}}, $database;	
		}
	} elsif ($mode eq "mirgenedb"){
		open my $IN, "< $in" or die "Please provide input file\n";
		my $tagsDB = load_species_names_mirgenedb();
		while($IN){
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
			push @{$complete_data{$acc}{"Specie"}}, $specie
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
		$relation{"$id\/$start-$end"} = $relation;
	}
	return \%relation;
}

sub load_species_names_mirgenedb {
	open my $DB, "< test_mirgenedb_species_cleaned.txt" or die "Please obtain the file test_mirgenedb_species_cleaned.txt\n";
	my %db;
	while (<$DB>){
		chomp;
		$_ =~ s/^(.*)(\()(.*)(\))(.*)\s+$/$3 $5/g;
		$db{$5} = $3;
	}
	return \%db;
}

exit;
