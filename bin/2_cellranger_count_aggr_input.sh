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

usage() {
    echo "Usage: $0 -s <sample_sheet> -i <count_info>"
    exit 1
}

while getopts ":s:i:" flag; do
    case $flag in
        s) sample_sheet="$OPTARG" ;;
        i) cellranger_count_info="$OPTARG" ;;
        *) usage ;;
    esac
done

# Check if inputs are provided and exist
[[ -z "${sample_sheet:-}" || -z "${cellranger_count_info:-}" ]] && usage
[[ ! -f "$sample_sheet" || ! -f "$cellranger_count_info" ]] && { echo "Error: Files not found."; exit 1; }

output_file="cellranger_aggr_input.csv"

#########################################################
# Section: 1. Load Count Info
# Description: Reads the count info file into an associative array
#              for fast lookup during the main processing loop.
# Data Structure:
#   count_dirs["SampleID"] = "OutputDirectoryPath"
#########################################################
declare -A count_dirs
while IFS=$'\t' read -r id outdir || [[ -n "$id" ]]; do
    count_dirs["$id"]="$outdir"
done < <(sed '1d' "$cellranger_count_info")

#########################################################
# Section: 2. Generate CSV Header
# Description: Extracts metadata column names (cols 3+) from the
#              sample_sheet header and creates the output CSV header.
# Output: Writes "sample_id,molecule_h5,..." to the output file.
#########################################################
metadata_header=$(awk -F'\t' 'NR==1 {for(i=3;i<=NF;i++) printf ",%s", $i; print ""}' "$sample_sheet")
echo "sample_id,molecule_h5${metadata_header}" > "$output_file"

#########################################################
# Section: 3. Main Processing Loop
# Description: Iterates through the sample sheet row by row.
# Logic:
#   1. Cleans the Sample ID.
#   2. Retrieves the output directory from the `count_dirs` array.
#   3. Validates that the barcode file exists and is not empty.
#   4. Formats metadata (converts tabs to commas).
#   5. Appends the formatted line to the output CSV.
#########################################################
while IFS=$'\t' read -r id col2 metadata_raw || [[ -n "$id" ]]; do
    # Skip if ID is empty or if it's the header (safety check)
    [[ -z "$id" || "$id" == "sample_id" ]] && continue

    # Clean sample ID (replace comma with underscore)
    sample_id_clean="${id//,/_}"

    # Get directory from our pre-loaded array
    cellranger_outdir="${count_dirs[$sample_id_clean]:-}"

    if [[ -z "$cellranger_outdir" ]]; then
        echo "Warning: No count info found for sample ${sample_id_clean}. Skipping."
        continue
    fi

    h5_path="${cellranger_outdir}/outs/molecule_info.h5"
    barcodes_gz="${cellranger_outdir}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"

    # Validate paths and barcode count
    if [[ -f "$barcodes_gz" ]]; then
        # Count lines
        n_barcodes=$(gunzip -c "$barcodes_gz" | wc -l)

        if (( n_barcodes > 0 )); then
            echo "Found ${n_barcodes} mapped barcodes for ${sample_id_clean}"
            
            # Convert remaining metadata tabs to commas
            metadata_csv="${metadata_raw//$'\t'/,}"
            
            # Append to file (only if metadata exists, add leading comma)
            [[ -n "$metadata_csv" ]] && metadata_csv=",${metadata_csv}"
            echo "${sample_id_clean},${h5_path}${metadata_csv}" >> "$output_file"
        else
            echo "Notice: No mapped barcodes found for ${sample_id_clean}"
        fi
    else
        echo "Error: Path not found for ${sample_id_clean}: ${barcodes_gz}"
    fi

done < <(sed '1d' "$sample_sheet")

echo "Success! ${output_file} has been generated."