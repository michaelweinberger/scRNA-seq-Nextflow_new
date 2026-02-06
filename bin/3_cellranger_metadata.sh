#!/bin/bash

#########################################################
### This script:
# creates a metadata file to use with Cell Ranger aligned scRNA-seq data
#########################################################

# === Safety & Error Handling ===
# set -e: Exit immediately if a command exits with a non-zero status.
# set -u: Treat unset variables as an error when substituting.
# set -o pipefail: The return value of a pipeline is the status of the last command to exit with a non-zero status.set -e
set -u
set -o pipefail

### declare input arguments
cellranger_aggr_input_csv=""

while getopts ":i:" flag
do
    case $flag in
        i) cellranger_aggr_input_csv="$OPTARG" ;;
        *) echo "Usage: $0 -i input.csv"; exit 1 ;;
    esac
done

if [[ -z "$cellranger_aggr_input_csv" ]]; then
    echo "Error: Input CSV not provided via -i"
    exit 1
fi

# === Output Header Generation ===
# 1. Reads the first line (header) of the input CSV.
# 2. Skips the first two standard Cell Ranger columns (sample_id, molecule_h5).
# 3. Extracts columns 3 to the end (custom metadata like 'treatment', 'time', etc.).
# 4. Writes 'barcode', 'sample_id', and the extracted headers (tab-separated) to the new output file.
metadata_headers=$(head -n 1 "$cellranger_aggr_input_csv" | awk -F, '{for(i=3;i<=NF;i++) printf "\t%s", $i}')
echo -e "barcode\tsample_id${metadata_headers}" > cellranger_aggr_cell_metadata.tsv

# === Sample Processing Loop ===
# Initialize a counter to track the sample index (required for barcode suffixing).
sample_count=1

# Read the CSV line by line, skipping the header (tail -n +2).
# IFS=',' splits the CSV line into an array named 'cols'.
tail -n +2 "$cellranger_aggr_input_csv" | while IFS=',' read -r -a cols
do
    # Assign variables based on CSV columns
    sample_id="${cols[0]}"
    h5_path="${cols[1]}"
    
    # === Metadata Extraction ===
    # Dynamically loops through columns 3+ of the current row to grab custom metadata values.
    # Adds a tab character before each value to align with the TSV format.
    extra_metadata=""
    for ((i=2; i<${#cols[@]}; i++)); do
        extra_metadata="${extra_metadata}\t${cols[$i]}"
    done

    # === Path Resolution Logic ===
    # This block attempts to locate the 'barcodes.tsv.gz' file.
    # It assumes 'h5_path' points to a molecule_info.h5 file.
    # Strategy 1: Look in 'sample_filtered_feature_bc_matrix' (common in 'multi' pipelines).
    barcodes_file="${h5_path%/*}/sample_filtered_feature_bc_matrix/barcodes.tsv.gz"

    if [[ ! -f "$barcodes_file" ]]; then
        # Strategy 2: Fallback for standard 'cellranger count' paths.
        barcodes_file="${h5_path%/*}/filtered_feature_bc_matrix/barcodes.tsv.gz"
    fi

    # === Barcode Modification & Writing ===
    if [[ -f "$barcodes_file" ]]; then
        echo "Processing Sample ${sample_count}: ${sample_id}"
        
        # 1. Unzip the barcode file.
        # 2. Use AWK to process each barcode:
        #    - gsub: Cell Ranger 'aggr' changes barcode suffixes to match the sample index.
        #      Example: 'AACTG...-1' becomes 'AACTG...-2' if it is the second sample processed.
        #      This command replaces the trailing "-1" with "-" followed by the current loop index.
        #    - print: Output the modified barcode, the sample ID, and any extra metadata.
        # 3. Append (>>) the result to the final metadata TSV.
        gunzip -c "$barcodes_file" | awk -v sid="$sample_id" -v idx="$sample_count" -v meta="$extra_metadata" \
        'BEGIN {FS=OFS="\t"} {
            gsub(/-1$/, "-"idx, $1); 
            print $1, sid meta
        }' >> cellranger_aggr_cell_metadata.tsv
        
    else
        echo "Error: Barcode file not found for ${sample_id} at ${barcodes_file}"
    fi

    sample_count=$((sample_count + 1))
done

echo "Success! cellranger_aggr_cell_metadata.tsv has been generated."