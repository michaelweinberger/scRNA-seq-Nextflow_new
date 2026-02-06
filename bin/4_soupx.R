#!/usr/bin/env Rscript


### Run this script to perform ambient RNA correction via SoupX



### User defined variables
# unpack variables passed from parent shell script
cli <- commandArgs(trailingOnly = TRUE) 
args <- strsplit(cli, "=", fixed = TRUE)

for (e in args) {
  argname <- e[1]
  argval <- e[2]
  # regular expression to delete initial \" and trailing \"
  argval <- gsub("(^\\\"|\\\"$)", "", argval)
  assign(argname, argval)
}



### Packages
library(SoupX)
library(DropletUtils)
library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(writexl)

set.seed(42)



### Functions

#' Batch SoupX Ambient RNA Correction Pipeline
#'
#' This function iterates through samples within a Seurat object, matches them to their 
#' original raw 10X data, and performs SoupX ambient RNA correction on a per-sample basis.
#' It generates corrected count matrices, contamination statistics, and quality control plots.
#'
#' @param seurat_obj A Seurat object containing the filtered cells to be corrected. 
#'                   Must contain clustering information (specifically 'seurat_clusters' in metadata).
#' @param meta_col A string representing the metadata column name in \code{seurat_obj} that holds 
#'                 sample identifiers (e.g., "orig.ident", "sample_id"). Used to split the object.
#' @param raw_paths_list A named vector or list where names match the sample IDs in \code{meta_col} 
#'                       and values are the file paths to the raw 10X feature barcode matrices 
#'                       (directories containing matrix.mtx, barcodes.tsv, features.tsv).
#' @param out_dir String. The path to the directory where results will be saved.
#' @param out_name String. A prefix tag used for naming all output files.
#'
#' @return A combined sparse matrix (dgCMatrix) containing the SoupX-corrected counts for all samples.
#'
#' @details 
#' The function performs the following steps for each sample:
#' 1. Subsets the Seurat object.
#' 2. Loads the corresponding Raw 10X counts.
#' 3. Aligns barcode suffixes between the Seurat object (often modified by merge) and raw data.
#' 4. Runs SoupX: calculates the "soup" profile, adds cluster info, estimates contamination (rho).
#' 5. Adjusts counts to remove ambient RNA.
#' 6. exports a combined Matrix, an Excel summary of genes removed, and QC plots (Rho estimation and Removal histograms).
#'
batch_soupx_pipeline <- function(seurat_obj, 
                                 meta_col, 
                                 raw_paths_list, 
                                 out_dir, 
                                 out_name
                                 ) {
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  samples <- unique(seurat_obj[[meta_col, drop = TRUE]])
  matrix_list <- list()
  prev_list <- list() 
  
  contamination_stats <- data.frame(
    Sample = character(),
    Contamination_Fraction = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (sample_id in samples) {
    message("\n>>> Processing Sample: ", sample_id)

    # Check if sample_id exists in the raw paths list
    if (!(sample_id %in% names(raw_paths_list))) {
      stop(paste0("Error: Sample ID '", sample_id, 
                  "' found in Seurat object but NOT found in the CSV file. ",
                  "Available IDs in CSV are: ", paste(names(raw_paths_list), collapse=", ")))
    }

    # 1. Subset and Load Data
    cell_names <- rownames(seurat_obj@meta.data)[seurat_obj[[meta_col, drop = TRUE]] == sample_id]
    sub_obj <- subset(seurat_obj, cells = cell_names)
    message("\n>>> Reading data from: ", raw_paths_list[[sample_id]])

    raw_counts <- Seurat::Read10X(raw_paths_list[[sample_id]])
    
    # Handle multi-modal data (GEX + Antibody/CRISPR)
    if (is.list(raw_counts)) {
      message("\n>>> Multiple feature types detected. Selecting 'Gene Expression'...")
      raw_counts <- raw_counts[["Gene Expression"]]
    }

    # Align Suffixes
    raw_barcodes <- colnames(raw_counts)
    match_suffix <- regexpr("-[0-9]+$", colnames(sub_obj)[1])
    
    if (match_suffix != -1) {
      sample_suffix <- regmatches(colnames(sub_obj)[1], match_suffix)
      colnames(raw_counts) <- sub("-[0-9]+$", sample_suffix, raw_barcodes)
    }
    message("Barcode names in raw counts are like: ", colnames(raw_counts)[1])
    message("Barcode names in filtered counts are like: ", colnames(sub_obj)[1])
    
    # 2. Pre-SoupX Gene Prevalence
    filtered_counts <- LayerData(sub_obj, assay = "RNA", layer = "counts")
    n_cells_pre <- Matrix::rowSums(filtered_counts > 0)
    
    # Ensure raw_counts has the same genes as filtered_counts
    common_genes <- intersect(rownames(raw_counts), rownames(filtered_counts))
    raw_counts <- raw_counts[common_genes, ]
    filtered_counts <- filtered_counts[common_genes, ]
    message("Number of genes in raw counts: ", length(rownames(raw_counts)))
    message("Number of genes in filtered counts: ", length(rownames(filtered_counts)))
    message("Number of genes after assimilating: ", length(common_genes))
    message("Gene names in raw counts are like: ", rownames(raw_counts)[1])
    message("Gene names in filtered counts are like: ", rownames(filtered_counts)[1])

    # 3. Run SoupX
    message("\n>>> Creating Soup channel.")
    sc <- SoupChannel(raw_counts, filtered_counts)
    if ("seurat_clusters" %in% colnames(sub_obj@meta.data)) {
      clusters <- setNames(sub_obj$seurat_clusters, rownames(sub_obj@meta.data))
      message("\n>>> Adding clustering information.")
      sc <- setClusters(sc, clusters)
    } else {
        stop("Please provide 'seurat_clusters' column in Seurat object metadata.")
    }
    
    sc <- autoEstCont(sc, doPlot = FALSE)
    rho <- sc$fit$rhoEst
    contamination_stats <- rbind(contamination_stats, 
                                 data.frame(Sample = sample_id, Contamination_Fraction = rho))
    
    # 4. Correct Counts and Post-SoupX Prevalence
    message("\n>>> Adjusting gene counts.")
    cleaned_mtx <- adjustCounts(sc, roundToInt = TRUE)

    # Ensure we only compare genes that exist in the cleaned matrix
    common_genes_post <- intersect(names(n_cells_pre), rownames(cleaned_mtx))
    n_cells_pre_aligned <- n_cells_pre[common_genes_post]
    n_cells_post <- Matrix::rowSums(cleaned_mtx[common_genes_post, ] > 0)
    
    # 5. Compute Per-Sample Prevalence & Fraction Removed
    df_prev <- data.frame(
      Gene = common_genes_post,
      Cells_Pre = as.numeric(n_cells_pre_aligned),
      Cells_Post = as.numeric(n_cells_post),
      Sample = sample_id
    )
    
    df_prev$Difference <- df_prev$Cells_Pre - df_prev$Cells_Post
    df_prev$Fraction_Removed <- ifelse(df_prev$Cells_Pre > 0, 
                                       df_prev$Difference / df_prev$Cells_Pre, 
                                       0)
    
    # Sort: Descending order by Fraction_Removed
    df_prev <- df_prev %>% arrange(desc(Fraction_Removed))

    # Save individual sample CSV
    sample_csv_name <- file.path(out_dir, paste0(out_name, "_", sample_id, "_gene_prevalence.csv"))
    write.csv(df_prev, file = sample_csv_name, row.names = FALSE)
    
    prev_list[[sample_id]] <- df_prev
    matrix_list[[sample_id]] <- cleaned_mtx
    
    rm(sub_obj, raw_counts, filtered_counts, sc)
    gc()
  }
  
  # --- FINAL EXPORTS ---
  # A. Save Concatenated 10X Matrix
  message("\n>>> Combining single sample results.")
  combined_matrix <- do.call(cbind, matrix_list)
  matrix_out_path <- file.path(out_dir, paste0(out_name, "_cleaned_10X"))
  write10xCounts(path = matrix_out_path, x = combined_matrix, overwrite = TRUE, version = "3")
  
  # B. Combined Summary EXCEL (One tab per sample)
  excel_path <- file.path(out_dir, paste0(out_name, "_TOTAL_gene_prevalence.xlsx"))
  message("\n>>> Saving multi-tab Excel file to: ", excel_path)
  write_xlsx(prev_list, path = excel_path)
  
  # C. Contamination Stats CSV
  write.csv(contamination_stats, file = file.path(out_dir, paste0(out_name, "_contamination_stats.csv")), row.names = FALSE)
  
  # D. PLOT 1: Contamination Bar Chart (PDF)
  p_rho <- ggplot(contamination_stats, aes(x = Sample, y = Contamination_Fraction, fill = Sample)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    geom_text(aes(label = scales::percent(Contamination_Fraction, accuracy = 0.1)), vjust = -0.5) +
    theme_minimal() +
    labs(title = paste(out_name, ": Contamination Fractions"), y = "rho", x = "Sample") +
    scale_y_continuous(labels = scales::percent, limits = c(0, max(contamination_stats$Contamination_Fraction) * 1.3))
  
  ggsave(file.path(out_dir, paste0(out_name, "_contamination_rho_plot.pdf")), plot = p_rho, width = 7, height = 5)
  
  # E. PLOT 2: Histogram of Fraction Removed (PDF)
  # We filter out genes with 0 expression pre-SoupX to avoid cluttering the plot
  full_prev_df <- dplyr::bind_rows(prev_list)
  p_hist <- ggplot(full_prev_df[full_prev_df$Cells_Pre > 0, ], aes(x = Fraction_Removed, fill = Sample)) +
    geom_histogram(bins = 50, color = "white", alpha = 0.8) +
    facet_wrap(~Sample) +
    theme_bw() +
    labs(title = paste(out_name, ": Gene-wise Fraction of Cells Removed"),
         subtitle = "Fraction of cells where gene expression was identified as ambient RNA",
         x = "Fraction of Cells Removed (0 = kept, 1 = fully removed)",
         y = "Number of Genes") +
    scale_x_continuous(labels = scales::percent) +
    theme(legend.position = "none")

  ggsave(file.path(out_dir, paste0(out_name, "_fraction_removed_histogram.pdf")), plot = p_hist, width = 10, height = 6)
  
  message("\n>>> All files saved to: ", out_dir)
  return(combined_matrix)
}



### Analysis
required_vars <- c("mapping_mode", "cellranger_count_info", "seurat_obj_file", "out_dir", "out_name", "harmony_var")
for (v in required_vars) {
  if (!exists(v)) stop(paste("Missing required argument:", v))
}

# create output directory if it does not exist
if (!dir.exists(out_dir)) {
    dir.create(out_dir)
}

seurat_obj <- readRDS(seurat_obj_file)
cellranger_df <- read.table(cellranger_count_info, sep = "\t", header=TRUE)

if (mapping_mode == "cell ranger count") {
    cellranger_df$cellranger_dir <- paste0(cellranger_df$cellranger_dir, "/outs/raw_feature_bc_matrix")
} else if (mapping_mode == "cell ranger multi") {
    cellranger_df$cellranger_dir <- paste0(cellranger_df$cellranger_dir, "/outs/multi/count/raw_feature_bc_matrix")
} else {
  stop("Argument mapping_mode must be set to 'cell ranger count' or 'cell ranger multi'")
}

raw_paths_list <- cellranger_df$cellranger_dir
names(raw_paths_list) <- cellranger_df$sample_id
raw_paths_list <- raw_paths_list[!duplicated(names(raw_paths_list))]

soupx_matrix <- batch_soupx_pipeline(seurat_obj=seurat_obj, 
                                     meta_col=harmony_var, 
                                     raw_paths_list=raw_paths_list, 
                                     out_dir=out_dir, 
                                     out_name=out_name
                                     )