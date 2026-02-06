#!/bin/bash

#########################################################
### This script:
# Generates an input .csv file for cellranger aggr
# based on a sample sheet and a count info file.
#########################################################


# Strict Mode Settings:
# set -e: Exit immediately if a command exits with a non-zero status.
# set -u: Treat unset variables as an error when substituting.
# set -o pipefail: Return value of a pipeline is the status of the last command to exit with a non-zero status.
set -e
set -u
set -o pipefail

### declare input arguments
sample_sheet=""
cellranger_count_info=""

while getopts ":s:i:" flag
do
    case $flag in
        s) sample_sheet="$OPTARG" ;;
        i) cellranger_count_info="$OPTARG" ;;
        *) echo "Usage: $0 -s sample_sheet -i count_info"; exit 1 ;;
    esac
done

if [[ ! -f "$sample_sheet" ]]; then echo "Error: Sample sheet not found."; exit 1; fi
if [[ ! -f "$cellranger_count_info" ]]; then echo "Error: Count info file not found."; exit 1; fi

# ==============================================================================
# ### HEADER GENERATION
# Purpose: Creates the 'cellranger_aggr_input.csv' file and writes the first line.
# Logic:   1. Identifies dynamic metadata headers starting from Column 4 of the sample sheet.
#          2. Adds standard required headers: sample_id, molecule_h5, parent_id.
# Outputs: A new file 'cellranger_aggr_input.csv' containing the header row.
# ==============================================================================
# We take sample_id (Col 3) as a reserved column, then metadata starting from Col 4
metadata_headers=$(head -n 1 "$sample_sheet" | awk -F'\t' '{for(i=4;i<=NF;i++) printf ",%s", $i}')
echo "sample_id,molecule_h5,parent_id${metadata_headers}" > cellranger_aggr_input.csv

# ==============================================================================
# ### PARENT ID (RUN) ITERATION
# Purpose: Identifies unique experimental runs (Parent IDs) to process.
# Logic:   1. Skips the header row (tail -n +2).
#          2. Extracts unique IDs from Column 1.
#          3. Locates the output directory for each run from the count_info file.
# ==============================================================================
tail -n +2 "$sample_sheet" | awk -F'\t' '{print $1}' | sort -u | while read -r run_id
do
    # Clean run ID (replace comma with underscore)
    run_id_clean="${run_id//,/_}"

    # Get output directory for this run
    cellranger_count_outdir=$(awk -F'\t' -v id="$run_id_clean" '$1 == id {print $2}' "$cellranger_count_info")

    if [[ -z "$cellranger_count_outdir" ]]; then
        echo "Warning: No output directory found for $run_id_clean"
        continue
    fi

    # ==============================================================================
    # ### SAMPLE EXTRACTION AND PATH CONSTRUCTION
    # Purpose: Processes specific samples belonging to the current Run ID.
    # Logic:   1. Filters sample sheet for the current Run ID.
    #          2. Extracts Sample ID (Col 3) and Metadata (Col 4+).
    #          3. Constructs file paths for 'sample_filtered_barcodes.csv' and 'sample_molecule_info.h5'.
    # ==============================================================================
    # Using awk to filter ensures we only match the first column
    awk -F'\t' -v rid="$run_id" '$1 == rid' "$sample_sheet" | while IFS=$'\t' read -r -a cols
    do
        parent_id="${cols[0]}"
        sample_id="${cols[2]}" # Column 3

        # Clean parent IDs (replace comma with underscore)
        parent_id_clean="${parent_id//,/_}"
        
        # Collect metadata from Column 4 onwards (index 3 in bash array)
        extra_metadata=""
        for ((i=3; i<${#cols[@]}; i++)); do
            extra_metadata="${extra_metadata},${cols[$i]}"
        done

        # Define paths
        # NOTE: Verify if your file is 'sample_molecule_info.h5' or just 'molecule_info.h5'
        barcodes_path="${cellranger_count_outdir}/outs/per_sample_outs/${sample_id}/count/sample_filtered_barcodes.csv"
        multi_h5_path="${cellranger_count_outdir}/outs/per_sample_outs/${sample_id}/count/sample_molecule_info.h5"

        # ==============================================================================
        # ### VALIDATION AND WRITE TO CSV
        # Purpose: Verifies data integrity and writes the final entry.
        # Logic:   1. Checks if the barcodes file exists.
        #          2. Checks if the barcodes file is not empty (contains > 1 line).
        #          3. Appends the formatted string to 'cellranger_aggr_input.csv'.
        # ==============================================================================
        if [ -f "$barcodes_path" ]; then
            n_barcodes=$(wc -l < "$barcodes_path")

            if [ "$n_barcodes" -gt 1 ]; then
                echo "Found ${n_barcodes} barcodes for ${sample_id}"
                echo "${sample_id},${multi_h5_path},${parent_id_clean}${extra_metadata}" >> cellranger_aggr_input.csv
            else
                echo "Skipping ${sample_id}: No barcodes found."
            fi
        else
            echo "Error: File not found for ${sample_id}: $barcodes_path"
        fi
    done
done

echo "------------------------------------------------"
echo "Success! cellranger_aggr_input.csv has been generated."