#!/bin/bash

function get_all() {
	wget $1
	return "Data fetched!"
}

mirbase_db="ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz" #Fasta file
rfam="ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/" #Fasta files with all annotated sequences
mirgenedb="http://mirgenedb.org/gff/ALL?sort=pos&all=1" #All annotated sequences, GFF file

get_all "$mirbase_db"
mkdir Fasta_RFAM; cd Fasta_RFAM
#get_all "$rfam/*.fa"
cd ../
get_all "$mirgenedb"

mv ALL\?sort\=pos\&all\=1 ALL.gff; gunzip hairpin.fa.gz;
