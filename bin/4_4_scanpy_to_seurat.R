#!/usr/bin/env Rscript



### Run this script to convert a Scanpy AnnData object (saved as an .h5ad file)
# into a Seurat object using the 'AnndataR' package.



### User defined variables
# unpack variables passed from parent shell script
cli <- commandArgs(trailingOnly = TRUE) 
args <- strsplit(cli, "=", fixed = TRUE)

# Parse command line arguments into variables (scanpy_obj, out_dir, out_name, etc.)
for (e in args) {
  argname <- e[1]
  argval <- e[2]
  # regular expression to delete initial \" and trailing \"
  argval <- gsub("(^\\\"|\\\"$)", "", argval)
  assign(argname, argval)
}



### Packages
library(Seurat)
library(anndataR)
library(hdf5r)



### Functions

#' Convert Scanpy AnnData (H5AD) to Seurat Object
#'
#' This function reads a Scanpy .h5ad file using the `anndataR` package and converts it 
#' into a Seurat object. It handles the mapping of data layers (checking for 'raw' data), 
#' subsets the object based on optional cell barcodes, adds external metadata, and 
#' configures default cluster identities.
#'
#' @param h5ad_path Character. The full file path to the source .h5ad file.
#' @param cell_barcodes Character vector (Optional). A list of cell barcodes to retain. 
#'        If NULL, all cells in the H5AD file are kept.
#' @param cell_metadata Data frame (Optional). A data frame containing metadata to 
#'        add to the Seurat object. Row names must match the cell barcodes.
#'
#' @return A Seurat object containing the expression data and metadata from the Scanpy object.
convert_scanpy_to_seurat_anndatar <- function(h5ad_path, cell_barcodes = NULL, cell_metadata = NULL) {
    
    if (!file.exists(h5ad_path)) {
        stop(paste("File not found:", h5ad_path))
    }
    
    cat(paste("--- Reading AnnData object from:", h5ad_path, "using anndataR (as='Seurat') ---\n"))

    # 1. Check if the .raw matrix is present using hdf5r.
    # The presence of the '/raw' group dictates the naming convention.
    has_raw <- tryCatch({
        h5_file <- hdf5r::H5File$new(h5ad_path, mode = "r")
        raw_exists <- "raw" %in% h5_file$ls()$name
        h5_file$close_all()
        raw_exists
    }, error = function(e) {
        cat(paste("Warning: Failed to check H5AD structure with hdf5r. Error:", e$message, "\n"))
        return(TRUE) # Fallback to default assumption
    })

    # 2. Perform the one-step conversion without explicit matrix mapping arguments.
    # Rely on anndataR's defaults, which correctly transfer matrices and meta-data.
    seurat_obj <- anndataR::read_h5ad(
        path = h5ad_path, 
        as = "Seurat"
    )
    
    cat("Seurat object successfully created from H5AD file.\n")

    # 3. Handle cell subsetting
    if (!is.null(cell_barcodes)) {
        cat(paste("Subsetting Seurat object to", length(cell_barcodes), "requested cell barcodes...\n"))
        
        all_cell_barcodes <- Seurat::Cells(seurat_obj)
        valid_barcodes <- intersect(cell_barcodes, all_cell_barcodes)
        
        if (length(valid_barcodes) == 0) {
            stop("No valid cell barcodes found for subsetting.")
        }
        
        seurat_obj <- seurat_obj[, valid_barcodes]
        
        cat(paste("Seurat object subsetted to", ncol(seurat_obj), "cells.\n"))
    }

    # 4. Add barcode metadata
    if (!is.null(cell_metadata)) {
        # Filter metadata to only include cells present in the (potentially subsetted) object
        common_cells <- intersect(rownames(cell_metadata), Cells(seurat_obj))
        cell_metadata <- cell_metadata[common_cells, , drop = FALSE]
        seurat_obj <- AddMetaData(object = seurat_obj, metadata = cell_metadata)
        cat(paste("Barcode metadata added to Seurat object.\n"))
    }

    # 5. Set project name and default identity
    Seurat::Project(seurat_obj) <- basename(h5ad_path)

    if ("leiden" %in% colnames(seurat_obj@meta.data)) {
        Seurat::Idents(seurat_obj) <- "leiden"
        cat("Set default identity to 'leiden'.\n")
        seurat_obj$seurat_clusters <- seurat_obj$leiden
    }

    # 5. create "counts" layer from "X" layer if anndata object does not contain "raw" group
    if (has_raw == FALSE) {
        count_data <- LayerData(seurat_obj, assay = "RNA", layer = "X")
        LayerData(seurat_obj, assay = "RNA", layer = "counts") <- count_data
        cat("Converted 'X' layer to 'counts' layer.\n")
    }

    cat("--- Conversion and processing complete. Seurat object returned ---\n")
    
    return(seurat_obj)
}



### Analysis
if (!exists("scanpy_obj") || !exists("out_dir") || !exists("out_name")) {
    stop("Required variables (scanpy_obj, out_dir, out_name) were not passed as arguments.")
}

# Create output directory if it doesn't exist
if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

# Subset cell barcodes in Scanpy object before conversion to Seurat if cell metadata file has been given
cell_metadata_df <- NULL
cell_barcodes_list <- NULL
if (exists("cell_metadata_csv") && !is.null(cell_metadata_csv) && file.exists(cell_metadata_csv)) {
    # Read row names (cell IDs) from the metadata CSV
    cell_metadata_df <- read.csv(cell_metadata_csv, header = TRUE, row.names = 1)
    cell_barcodes_list <- rownames(cell_metadata_df)
    cat(paste("Identified", length(cell_barcodes_list), "cell barcodes for subsetting.\n"))
}

# Perform the conversion
seurat_obj <- convert_scanpy_to_seurat_anndatar(
    h5ad_path = scanpy_obj, 
    cell_barcodes = cell_barcodes_list,
    cell_metadata = cell_metadata_df
)

# Save the final Seurat object
print(seurat_obj)
print(head(seurat_obj[[]]))
saveRDS(seurat_obj, file.path(out_dir, paste0(out_name, ".rds")))

cat(paste0("Seurat object saved to disk as ", file.path(out_dir, paste0(out_name, ".rds")), ".\n"))