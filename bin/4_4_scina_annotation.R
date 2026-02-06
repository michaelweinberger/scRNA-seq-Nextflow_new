#!/usr/bin/env Rscript



### Run this script to perform automated cell type annotation with SCINA



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
library(Seurat)
library(readr)
library(SCINA)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(writexl)
library(stringr)



### Functions

#' Perform Automated Cell Type Annotation using SCINA
#'
#' @description
#' This function takes a Seurat object containing raw count data, performs
#' normalization and log-transformation using Seurat's \code{NormalizeData}, 
#' formats the marker genes, runs the SCINA algorithm, and exports prediction 
#' results and diagnostic plots.
#' 
#' @details
#' The function performs the following steps:
#' 1. Validates the existence of the output directory and raw count data.
#' 2. Normalizes the data using \code{Seurat::NormalizeData} (LogNormalize, scale 10k).
#' 3. Converts the markers dataframe into the list format required by SCINA.
#' 4. Runs the SCINA algorithm.
#' 5. Generates outputs:
#'    - A CSV of predicted labels and confidence scores.
#'    - A histogram of confidence scores.
#'    - A correlation matrix and heatmap of cell type profiles.
#'
#' @param seurat_object A Seurat object initialized with raw counts in the specified assay/slot.
#' @param markers_df A dataframe where columns represent Cell Types and rows represent Marker Genes.
#' @param out_dir String. The directory where results and plots will be saved.
#' @param out_name String. A prefix for the output filenames.
#' @param assay_name String. The name of the assay containing the raw counts (default is "RNA").
#' @param slot_name String. The name of the slot/layer containing the raw counts (default is "counts").
#'
#' @return A named list containing:
#'   \item{results_df}{A dataframe containing cell barcodes, SCINA predicted labels, and confidence scores.}
#'   \item{markers_list}{The formatted list of marker genes used for the annotation.}
#'
#' @export
seurat_scina_annotate <- function(seurat_object, markers_df, out_dir, out_name, assay_name = "RNA", slot_name = "counts") {

  # 1. Input Validation and Initial Check
  cat("Starting SCINA annotation...\n")

  # 1.1 Create Output Directory if it doesn't exist
  if (!dir.exists(out_dir)) {
    cat(paste("Creating output directory:", out_dir, "\n"))
    dir.create(out_dir, recursive = TRUE)
  }

  # 1.2 Check for raw counts presence before normalization
  cat("Checking for raw count matrix in the Seurat object...\n")
  tryCatch({
    # We don't extract it fully yet, but ensure it's accessible for NormalizeData
    check_matrix <- LayerData(seurat_object, assay = assay_name, layer = slot_name)
    if (nrow(check_matrix) == 0 || ncol(check_matrix) == 0) {
      stop("Raw count matrix is empty or missing in the specified assay/slot.")
    }
  }, error = function(e) {
    stop(paste("Error accessing raw count data:", e$message))
  })

  # 2. Normalization and Log-Transformation using Seurat
  cat("Normalizing and log-transforming raw counts using Seurat::NormalizeData...\n")

  # Ensure the assay is set as the active assay for Seurat functions to work smoothly
  DefaultAssay(seurat_object) <- assay_name

  # create layer "counts" needed for NormalizeData, if slot_name is not "counts"
  if (slot_name != "counts") {
    count_data <- LayerData(seurat_object, assay = assay_name, layer = slot_name)
    LayerData(seurat_object, assay = assay_name, layer = "counts") <- count_data
  }

  # Run Seurat's normalization. This updates the 'data' slot of the assay.
  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)

  # Extract the resulting normalized matrix from the 'data' slot
  normalized_matrix <- LayerData(seurat_object, assay = assay_name, layer = "data")

  # SCINA requires the matrix to have cells (barcodes) as columns and genes as rows.
  # Seurat matrices are typically in this format.

  # 3. Convert the tibble/dataframe into a named list required by SCINA
  # The list format should be: list(CellTypeA = c("Marker1", "Marker2"), CellTypeB = c("Marker3", "Marker4"))
  markers_list <- list()
  for (col_name in colnames(markers_df)) {
    # Remove NA values and convert to character vector
    markers <- na.omit(markers_df[[col_name]])
    markers_list[[col_name]] <- as.character(markers)
  }

  # The status message is updated to list the cell types found, instead of just counting them.
  cat(paste("Successfully loaded markers for", paste(names(markers_list), collapse=", "), ".\n"))

  # 4. Run SCINA Annotation
  cat("Running SCINA algorithm (this may take a few minutes for large datasets)...\n")
  scina_results <- SCINA(
    normalized_matrix, # Use the normalized matrix extracted after Seurat's NormalizeData
    markers_list,
    # log_file = "SCINA_log.txt", # Optional: logs the process
    allow_unknown = TRUE,
    rm_overlap = FALSE
  )

  cat("SCINA analysis complete.\n")
  
  # 5. Save SCINA Results and Plot Heatmap
  # 5.1 Extract and save SCINA results dataframe
  # Use apply() to iterate over columns (MARGIN = 2)
  #    - max() returns the maximum value in the current column.
  #    - na.rm = TRUE ensures that if there are any NA values, max() still returns the highest non-NA number.
  max_values <- apply(
    X = as.data.frame(scina_results$probabilities), 
    MARGIN = 2, 
    FUN = max, 
    na.rm = TRUE
  )

  results_df <- data.frame(
      cell_barcode = colnames(scina_results$probabilities),
      SCINA_predicted_celltype = scina_results$cell_labels,
      SCINA_confidence = max_values,
      row.names = NULL
  )
  results_file <- file.path(out_dir, paste0(out_name, "_scina_predictions.csv"))
  cat(paste("Saving SCINA prediction results to:", results_file, "\n"))
  write.csv(results_df, results_file, row.names = FALSE)
  
  # plot histogram of SCINA confidence values
  gg <- ggplot(results_df, aes(x = SCINA_confidence)) +
    geom_histogram(binwidth = 0.05, # Adjust binwidth as needed, e.g., 0.05 or 0.1
                   fill = "steelblue") +
    labs(title = "Distribution of SCINA Confidence Scores",
         x = "SCINA Confidence [0, 1]",
         y = "Frequency (Number of Cells)") +
    theme_minimal() + # A clean background theme
    # Set x-axis limits and breaks to cover the range [0, 1]
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1))

  ggsave(gg, filename=file.path(out_dir, paste0(out_name, "_scina_confidence_histo.pdf")), 
         device="pdf", dpi=300, useDingbats=FALSE)

  # 5.2 Generate and save SCINA marker heatmap using plotheat.SCINA()
  #heatmap_file <- file.path(out_dir, paste0(out_name, "_scina_marker_heatmap.png"))
  #cat(paste("Saving SCINA marker heatmap to:", heatmap_file, "\n"))

  # Plotting function must be wrapped in a graphics device call (e.g., pdf or png)
  #png(heatmap_file, width = 10, height = 7)
  # The plotheat.SCINA function plots a heatmap of marker gene expression across
  # predicted cell types.
  #plotheat.SCINA(normalized_matrix, scina_results, markers_list)
  #dev.off()
  
  # 5.3 Generate and save cell type correlation heatmap (based on SCINA probabilities)
  cat("Calculating and saving cell type correlation heatmap...\n")
  
  # Calculate Pearson correlation between cell type profiles based on their probability scores
  correlation_matrix <- cor(t(scina_results$probabilities))
  
  # Save the correlation matrix itself
  cor_file <- file.path(out_dir, paste0(out_name, "_scina_celltype_correlation.csv"))
  cat(paste("Saving cell type correlation matrix to:", cor_file, "\n"))
  write.csv(as.data.frame(correlation_matrix), cor_file, row.names = TRUE)
  
  # Plot the correlation heatmap
  cor_plot_file <- file.path(out_dir, paste0(out_name, "_scina_celltype_correlation_plot.pdf"))
  cat(paste("Saving correlation heatmap to:", cor_plot_file, "\n"))
  
  # Plotting function must be wrapped in a graphics device call
  pdf(cor_plot_file, width = 10, height = 10)
  pheatmap(
    correlation_matrix,
    display_numbers = TRUE,
    number_format = "%.2f",
    color = colorRampPalette(c("blue", "white", "red"))(50),
    main = paste(out_name, "\nCell Type Correlation (SCINA Scores)"),
    # Remove borders for a cleaner look
    border_color = NA
  )
  dev.off()

  cat("SCINA run complete.\n")
  
  # 6. Return the results
  return(list(results_df = results_df, markers_list = markers_list))
}



### Analysis
required_vars <- c("seurat_obj_file", "cell_type_markers_csv", "out_dir", "out_name")
for (v in required_vars) {
  if (!exists(v)) stop(paste("Missing required argument:", v))
}

# 1. Run the annotation function
seurat_obj <- readRDS(seurat_obj_file)
cell_type_markers_df <- read.csv(cell_type_markers_csv, header = TRUE) 

layers <- Layers(seurat_obj[["RNA"]])
if ("X" %in% layers && !"counts" %in% layers) {
    cat("Annotating using layer X.\n")
    scina_results_list <- seurat_scina_annotate(
        seurat_object = seurat_obj,
        markers_df = cell_type_markers_df,
        out_dir = out_dir,
        out_name = out_name,
        slot_name = "X"
        )
} else {
    cat("Annotating using layer counts.\n")
    scina_results_list <- seurat_scina_annotate(
        seurat_object = seurat_obj,
        markers_df = cell_type_markers_df,
        out_dir = out_dir,
        out_name = out_name
        )
}