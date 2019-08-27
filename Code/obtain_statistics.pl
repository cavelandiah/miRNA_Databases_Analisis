#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;


my $mirbase = "";
my $rfam = "";
my $mirgenedb = "";

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
	} elsif ($mode eq "rfam"){
	} elsif ($mode eq "mirgenedb"){
	} else {
		die "$mode is not valid\n";
	}
	return;
}	
