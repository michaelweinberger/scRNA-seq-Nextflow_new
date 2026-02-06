#!/bin/bash

#########################################################
### This script:
# generates an input .csv file for cellranger multi
#########################################################

### Safety Settings
# set -e: Exit immediately if a command exits with a non-zero status.
# set -u: Treat unset variables as an error when substituting.
# set -o pipefail: Return value of a pipeline is the status of the last command to exit with a non-zero status.
set -e
set -u
set -o pipefail

### declare input arguments
while getopts ":i:s:r:" flag
do
	case $flag in
        	i) 
			sample_id="$OPTARG"
			echo "Option -i with argument: $sample_id"
			;;
               	s) 
			sample_sheet="$OPTARG"
			echo "Option -s with argument: $sample_sheet"
			;;
        	r) 
			cellranger_index="$OPTARG"
			echo "Option -r with argument: $cellranger_index"
			;;
	esac
done

### Initialize Config File & [gene-expression] Section
# 1. Creates (overwrites) 'cellranger_multi_config.csv'.
# 2. Defines the [gene-expression] header.
# 3. Sets the reference genome path provided by the -r argument.
# 4. Enables BAM file creation.
printf "[gene-expression]\n" > cellranger_multi_config.csv
printf "reference,${cellranger_index}\n" >> cellranger_multi_config.csv
printf "create-bam,TRUE\n" >> cellranger_multi_config.csv

### ID Processing & Directory Lookup
# 1. Stores the Gene Expression ID.
# 2. Creates the CMO (Cell Multiplexing Oligo) ID by performing a string substitution 
#    on $sample_id (replacing the text "GEX" with "CMO").
# 3. Searches the sample sheet for the $sample_id and uses awk to extract the 2nd column
#    (which is assumed to be the FASTQ directory path).
sample_id_gex=$sample_id
sample_id_cmo="${sample_id//GEX/CMO}"
fastq_dir=$(grep $sample_id $sample_sheet | awk 'NR==1{print $2}')   # only print first entry in field

### [libraries] Section Generation
# This block writes the map of FASTQ files to the CSV.
# It defines two rows:
# 1. The Gene Expression data (using the GEX ID).
# 2. The Multiplexing Capture data (using the calculated CMO ID).
printf "[libraries]\n" >> cellranger_multi_config.csv
printf "fastq_id,fastqs,feature_types\n" >> cellranger_multi_config.csv
printf "${sample_id_gex},${fastq_dir},Gene Expression\n" >> cellranger_multi_config.csv
printf "${sample_id_cmo},${fastq_dir},Multiplexing Capture\n" >> cellranger_multi_config.csv

### [samples] Section Generation
# This block defines the actual samples for demultiplexing.
# 1. Writes the headers.
# 2. Greps the sample sheet again for the specific sample ID.
# 3. Uses awk to print columns 3 (Sample Name) and 4 (CMO ID/Tag) separated by a comma.
#    This links the physical CMO tag to the biological sample name.
printf "[samples]\n" >> cellranger_multi_config.csv
printf "sample_id,cmo_ids\n" >> cellranger_multi_config.csv	

# print de-multiplexed sample IDs and CMO IDs
grep $sample_id $sample_sheet | awk '{print $3","$4}' >> cellranger_multi_config.csv

echo "------------------------------------------------"
echo "Success! cellranger_multi_config.csv has been generated."