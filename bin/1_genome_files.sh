#!/bin/bash

### Safety and Debugging Settings
# set -e: Exit immediately if a command exits with a non-zero status.
# set -u: Treat unset variables as an error when substituting.
# set -o: Prints current options (Note: Usually this is used with 'pipefail', e.g., 'set -o pipefail')
set -e
set -u
set -o

### declare input arguments
while getopts ":s:g:u:l:e:o:" flag
do
	case $flag in
        	s) 
			species="$OPTARG"
			echo "Option -s with argument: $species"
			;;
        	g) 
			genome="$OPTARG"
			echo "Option -g with argument: $genome"
			;;
        	u) 
			genome_ucsc="$OPTARG"
			echo "Option -u with argument: $genome_ucsc"
			;;
        	l) 
			species_latin="$OPTARG"
			echo "Option -l with argument: $species_latin"
			;;
        	e) 
			ensembl_version="$OPTARG"
			echo "Option -e with argument: $ensembl_version"
			;;
        	o) 
			outdir="$OPTARG"
			echo "Option -o with argument: $outdir"
			;;
	esac
done

### Variable Formatting
# Uses Bash parameter expansion to capitalize the first letter of the species name.
# Required because Ensembl file naming conventions often capitalize the Genus (e.g., 'Homo_sapiens').
species_latin_2=${species_latin^}

### Download and Process Genome Sequence (FASTA)
# 1. Uses rsync to download the .dna.toplevel.fa.gz file from the EBI Ensembl FTP server.
# 2. Renames the downloaded file to a shorter, standard name.
# 3. Decompresses the file using gunzip.
rsync -avzP rsync://ftp.ebi.ac.uk/ensemblorg/pub/release-${ensembl_version}/fasta/${species_latin}/dna/${species_latin_2}.${genome}.dna.toplevel.fa.gz .
mv ./${species_latin_2}.${genome}.dna.toplevel.fa.gz ./${genome}.fa.gz
gunzip ./${genome}.fa.gz

### Download and Process Gene Annotations (GTF)
# 1. Uses rsync to download the gene transfer format (GTF) file from Ensembl.
# 2. Renames the file to include the genome and Ensembl version.
# 3. Decompresses the file.
rsync -avzP rsync://ftp.ebi.ac.uk/ensemblorg/pub/release-${ensembl_version}/gtf/${species_latin}/${species_latin_2}.${genome}.${ensembl_version}.gtf.gz .
mv ./${species_latin_2}.${genome}.${ensembl_version}.gtf.gz ./${genome}.${ensembl_version}.gtf.gz
gunzip ./${genome}.${ensembl_version}.gtf.gz

# download file containing positions of repetitive elements in genome
#wget -L http://hgdownload.soe.ucsc.edu/goldenPath/${genome_ucsc}/database/rmsk.txt.gz .
#gunzip ./rmsk.txt.gz
#mv ./rmsk.txt ./${genome}_rmsk.txt

echo "------------------------------------------------"
echo "Success! Genome FASTA and GTF files have been downloaded."
