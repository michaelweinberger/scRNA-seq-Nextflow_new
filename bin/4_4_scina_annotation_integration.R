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

#' Generate a Custom Annotated UMAP Plot
#'
#' This function creates a UMAP plot where cluster labels are placed dynamically
#' outside the cluster centroids (to avoid obscuring data points) and styled with
#' specific colors and a white outline for readability.
#'
#' @param seurat_obj A Seurat object containing the data and UMAP reduction.
#' @param group_key A string. The column name in metadata to group cells by (e.g., "ID").
#' @param color_key A string. The key in `seurat_obj@misc` containing the color mapping list.
#' @param out_dir A string. The directory where the plot will be saved.
#' @param out_name A string. The base name for the output file.
#' @param resolution A numeric value. The clustering resolution used, for file naming purposes.
#' @return NULL. The function saves a PDF file to the disk as a side effect.
plot_umap_with_annotations <- function(seurat_obj, group_key, color_key, out_dir, out_name, resolution) {
  cat(paste0("\n--- Generating UMAP plot with external annotations for ", group_key, " ---\n"))
  
  if (!("umap" %in% names(seurat_obj@reductions))) {
    warning("UMAP reduction not found in the Seurat object. Skipping annotated UMAP plot.")
    return(NULL)
  }
  
  # 1. Prepare Data
  umap_coords <- Embeddings(seurat_obj, reduction = "umap")
  plot_df <- as.data.frame(umap_coords)
  colnames(plot_df) <- c("UMAP1", "UMAP2")
  
  # Check if group key exists and add group metadata
  if (!group_key %in% colnames(seurat_obj@meta.data)) {
      warning(paste0("Group key '", group_key, "' not found in Seurat object metadata. Skipping plot."))
      return(NULL)
  }
  plot_df[[group_key]] <- seurat_obj@meta.data[[group_key]]
  
  # Get colors from Seurat object's miscellaneous slot
  color_map <- NULL
  if (color_key %in% names(seurat_obj@misc)) {
    colors_vec <- seurat_obj@misc[[color_key]]$colors
    categories <- seurat_obj@misc[[color_key]]$categories
    # Ensure categories and colors are aligned
    if (length(categories) == length(colors_vec)) {
        color_map <- setNames(colors_vec, categories)
    } else {
        warning(paste0("Mismatch in colors and categories for '", color_key, "'. Using default colors."))
    }
  } else {
    warning(paste0("Color key '", color_key, "' not found in seurat_obj@misc. Using default colors."))
  }
  
  # 2. Calculate Cluster Centroids and Spread
  group_stats <- plot_df %>%
    dplyr::group_by(!!sym(group_key)) %>%
    dplyr::summarise(
      cx = mean(UMAP1, na.rm = TRUE),
      cy = mean(UMAP2, na.rm = TRUE),
      sx = sd(UMAP1, na.rm = TRUE),
      sy = sd(UMAP2, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    # Replace NA std dev (e.g., single cell cluster) with a small, non-zero value
    dplyr::mutate(
      sx = replace(sx, is.na(sx), 0.5), 
      sy = replace(sy, is.na(sy), 0.5)
    )
  
  # Calculate the global center
  global_center_x <- mean(plot_df$UMAP1, na.rm = TRUE)
  global_center_y <- mean(plot_df$UMAP2, na.rm = TRUE)
  
  shift_multiplier <- 2.0 # Adjust this value to push labels further out/in
  
  # 3. Calculate External Annotation Position (Vector Calculation)
  annot_df <- group_stats %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      # 1. Direction vector from global center to cluster centroid
      dx = cx - global_center_x,
      dy = cy - global_center_y,
      
      # 2. Determine shift distance based on cluster spread
      # Use the maximum standard deviation to determine the shift scale, ensuring it clears the cluster's widest axis.
      shift_distance = shift_multiplier * max(sx, sy),
      
      # 3. Normalize direction vector
      norm = sqrt(dx^2 + dy^2),
      
      # Handle case where norm is 0 (cluster exactly at global center)
      dx_n = ifelse(norm > 0, dx / norm, 1.0), 
      dy_n = ifelse(norm > 0, dy / norm, 1.0),
      
      # 4. Calculate final annotation position
      x_annot = cx + dx_n * shift_distance,
      y_annot = cy + dy_n * shift_distance
    ) %>%
    dplyr::ungroup()
  
  # 4. Generate Plot
  # Base DimPlot (no labels, no legend)
  gg <- DimPlot(
    seurat_obj, 
    reduction = "umap", 
    group.by = group_key, 
    cols = color_map,
    label = FALSE, 
    pt.size = 0.4
  ) +
    NoLegend() +
    ggplot2::geom_text(
      data = annot_df, 
      aes(x = x_annot, y = y_annot, label = !!sym(group_key)),
      color = if (!is.null(color_map)) {
          color_map[as.character(annot_df[[group_key]])]
      } else {
          "black"
      },
      size = 3.5, 
      fontface = "bold",
      inherit.aes = FALSE
    ) +
    # Re-draw the labels with a white outline
    ggplot2::geom_text(
      data = annot_df, 
      aes(x = x_annot, y = y_annot, label = !!sym(group_key)),
      color = "white",
      size = 4, # Slightly larger for the outline
      fontface = "bold",
      inherit.aes = FALSE,
      alpha = 0.8
    ) +
    # Redraw the actual labels on top
    ggplot2::geom_text(
      data = annot_df, 
      aes(x = x_annot, y = y_annot, label = !!sym(group_key)),
      color = if (!is.null(color_map)) {
          color_map[as.character(annot_df[[group_key]])]
      } else {
          "black"
      },
      size = 3.5, 
      fontface = "bold",
      inherit.aes = FALSE
    ) +
    ggplot2::ggtitle(paste0(out_name, " UMAP Annotated (", group_key, ")")) +
    Seurat::NoAxes()

  # 5. Save the plot
  # Sanitize name for seurat_clusters case
  ident_name <- gsub("seurat_clusters", paste0("res_0_", resolution * 100), group_key)
  plot_file <- file.path(out_dir, paste0(out_name, "_umap_", ident_name, "_annotated.pdf"))
  
  ggplot2::ggsave(
    gg, 
    filename = plot_file, 
    device = "pdf", 
    dpi = 300, 
    useDingbats = FALSE,
    width = 8, 
    height = 8
  )
  cat(paste("Annotated UMAP plot saved to:", plot_file, "\n"))
  
  return(NULL)
}



#' Perform Automated Cell Type Annotation using SCINA
#'
#' This function takes a Seurat object and a dataframe of marker genes, normalizes the data,
#' runs the SCINA algorithm to predict cell types, and generates initial confidence plots.
#'
#' @param seurat_object A Seurat object containing raw counts.
#' @param markers_df A dataframe where columns are cell types and rows are marker genes.
#' @param out_dir A string. Directory to save outputs.
#' @param out_name A string. Base name for output files.
#' @param assay_name A string. The Seurat assay to use (default "RNA").
#' @param slot_name A string. The slot containing raw counts (default "counts").
#' @return A list containing:
#'   \item{results_df}{A dataframe of SCINA predictions and confidence scores.}
#'   \item{markers_list}{A list format of the input markers used for SCINA.}
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
         x = "SCINA Confidence [0, 1])",
         y = "Frequency (Number of Cells)") +
    theme_minimal() + # A clean background theme
    # Set x-axis limits and breaks to cover the range [0, 1]
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1))

  ggsave(gg, filename=file.path(out_dir, paste0(out_name, "_scina_confidence_histo.pdf")), 
         device="pdf", dpi=300, useDingbats=FALSE)

  # 5.2 Generate and save cell type correlation heatmap (based on SCINA probabilities)
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
  pdf(cor_plot_file, width = 8, height = 8)
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
  
  # 7. Return the results
  return(list(results_df = results_df, markers_list = markers_list))
}



#' Generate UMAP Feature Plots for Marker Genes
#'
#' Iterates through the cell types provided in `markers_list`, finds the corresponding
#' marker genes that exist in the Seurat object, and generates a composite FeaturePlot
#' (gene expression on UMAP) for each cell type.
#'
#' @param seurat_object A Seurat object with a valid UMAP reduction.
#' @param markers_list A named list where names are cell types and values are character vectors of gene names.
#' @param out_dir A string. Directory to save plots.
#' @param out_name A string. Base name for output files.
#' @return NULL. Saves PNG files to the disk.
plot_umap_feature_plots <- function(seurat_object, markers_list, out_dir, out_name) {
  cat("\n--- Generating UMAP FeaturePlots for cell-type-specific marker genes ---\n")

  # Ensure the UMAP reduction exists before trying to plot
  if (!("umap" %in% names(seurat_object@reductions))) {
    warning("UMAP reduction not found in the Seurat object. Skipping UMAP FeaturePlot generation. Please run 'RunUMAP()' on your Seurat object before calling this function.")
    return(NULL)
  }

  # Get all genes present in the Seurat object
  all_genes_in_seurat <- rownames(seurat_object)

  for (cell_type in names(markers_list)) {

    # Get markers for the current cell type
    markers_for_type <- markers_list[[cell_type]]

    # Filter markers: keep only those present in the Seurat object
    valid_markers <- markers_for_type[markers_for_type %in% all_genes_in_seurat]

    if (length(valid_markers) > 0) {

      # Determine the number of features for optimal layout
      n_features <- length(valid_markers)
      n_cols <- min(4, n_features) # Plot up to 4 columns
      n_rows <- ceiling(n_features / n_cols)

      # Calculate appropriate width and height (e.g., 6 inches per column/row)
      plot_width <- 6 * n_cols
      plot_height <- 6 * n_rows

      # Sanitize cell_type name for file saving
      sanitized_cell_type <- gsub("[^A-Za-z0-9_]", "_", cell_type)
      plot_file <- file.path(out_dir, paste0(out_name, "_UMAP_Markers_", sanitized_cell_type, ".png"))

      cat(paste("Saving UMAP FeaturePlots for", cell_type, "to:", plot_file, "\n"))

      # Generate the FeaturePlot
      umap_plot <- FeaturePlot(
        seurat_object,
        features = valid_markers,
        reduction = "umap",
        combine = TRUE,
        ncol = n_cols
      ) + ggplot2::ggtitle(paste("UMAP FeaturePlot for", cell_type, "Markers"))

      # Save the plot as PNG
      ggplot2::ggsave(
        umap_plot,
        filename = plot_file,
        device = "png",
        width = plot_width,
        height = plot_height,
        units = "in",
        dpi = 300
      )
    } else {
      warning(paste("No valid marker genes found for cell type:", cell_type, ". Skipping UMAP FeaturePlot."))
    }
  }
  cat("Finished generating UMAP FeaturePlots for cell-type markers.\n")
}



#' Integrate SCINA results and plot cell type frequency per cluster
#'
#' This function merges SCINA predictions into the Seurat metadata, calculates the
#' consensus cell type for each unsupervised cluster, standardizes naming conventions,
#' finds differential markers for the new annotations, and generates summary visualizations
#' (Bar charts, Heatmaps, Violin plots).
#'
#' @param seurat_object A Seurat object.
#' @param scina_results_df Dataframe containing SCINA predictions.
#' @param cell_type_meta_df Dataframe with metadata for ordering and colors (cols: cell_type, ID, color).
#' @param out_dir String. Output directory.
#' @param out_name String. Base output filename.
#' @param n_markers Integer. Number of markers to save to Excel.
#' @param heatmap_n Integer. Number of markers to include in the Heatmap.
#' @param violin_n Integer. Number of markers to include in the Violin plot.
#' @param umap_marker_n Integer. Number of consensus markers to return for UMAP plotting.
#' @param leiden_res Numeric. The resolution used for Leiden/Seurat clustering (for naming).
#' @return A list containing:
#'   \item{seurat_object}{A modified Seurat object with the SCINA annotation results added.}
#'   \item{umap_markers_list}{A named list of the top UMAP markers (genes) for each SCINA consensus label.}
#' @export
integrate_scina_results <- function(seurat_object, scina_results_df, cell_type_meta_df, out_dir, out_name, n_markers, heatmap_n, violin_n, umap_marker_n, leiden_res) {
  
  cat("\n--- Integrating SCINA labels and calculating cluster frequencies ---\n")
  
  # --- 1. Integrate Results back into Seurat Object Metadata ---
  predicted_labels <- setNames(scina_results_df$SCINA_predicted_celltype, scina_results_df$cell_barcode)

  if (!all(colnames(seurat_object) %in% names(predicted_labels))) {
    warning("Predicted cell labels do not match the Seurat object cell barcodes.")
  }
  
  predicted_labels <- predicted_labels[colnames(seurat_object)]

  seurat_object <- AddMetaData(
    object = seurat_object,
    metadata = predicted_labels,
    col.name = "SCINA_predicted_celltype"
  )
  
  cat("SCINA labels successfully added to Seurat object metadata.\n")

  # --- 2. Check for Clusters and Prepare Data ---
  if (!"seurat_clusters" %in% colnames(seurat_object@meta.data)) {
    warning("Metadata column 'seurat_clusters' not found. Skipping frequency analysis and consensus labeling.")
    return(list(seurat_object = seurat_object, umap_markers_list = list()))
  }
  
  meta_data <- seurat_object@meta.data
  
  counts_df <- meta_data %>%
    dplyr::group_by(seurat_clusters, SCINA_predicted_celltype) %>%
    dplyr::summarise(Count = dplyr::n(), .groups = 'drop')
  
  freq_df <- counts_df %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::mutate(Frequency = Count / sum(Count)) %>%
    dplyr::ungroup()
  
  # --- 3. DETERMINE AND INTEGRATE CONSENSUS LABEL ---
  consensus_labels <- freq_df %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::slice_max(order_by = Frequency, n = 1, with_ties = FALSE) %>% 
    dplyr::ungroup() %>%
    dplyr::select(seurat_clusters, SCINA_predicted_celltype) %>%
    dplyr::rename(cell_type = SCINA_predicted_celltype)
  
  old_metadata <- seurat_object[[]]
  new_metadata <- old_metadata %>% left_join(consensus_labels, by = join_by(seurat_clusters))
  rownames(new_metadata) <- Cells(seurat_object)

  seurat_object[[]] <- new_metadata

  cat("Consensus SCINA label ('cell_type') added to Seurat object metadata.\n")
  
  # Re-level 'cell_type' based on cell_type_meta_df
  cat(paste("Re-leveling 'cell_type' categories based on 'cell_type' column order in metadata file...\n"))

  # 1. Standardize labels in Seurat object (replace '_' with ' ')
  seurat_object$cell_type <- as.character(seurat_object$cell_type)
  seurat_object$cell_type <- str_replace_all(seurat_object$cell_type, "_", " ")

  # 2. Standardize labels in metadata file
  cell_type_meta_df$cell_type <- as.character(cell_type_meta_df$cell_type)
  cell_type_meta_df$cell_type <- str_replace_all(cell_type_meta_df$cell_type, "_", " ")
  
  # Get the desired order from the metadata file
  new_cell_type_order <- unique(cell_type_meta_df$cell_type)

  # Filter the order to only include cell types present in the Seurat object
  present_categories <- unique(seurat_object$cell_type)
  
  final_order <- new_cell_type_order[new_cell_type_order %in% present_categories]

  # 3. Re-level the 'cell_type' column
  if (length(final_order) > 0) {
    seurat_object$cell_type <- factor(seurat_object$cell_type, levels = final_order, ordered = TRUE)
    cat(paste0("'cell_type' categories re-leveled to: ", paste(head(final_order, 5), collapse=", "), "... (", length(final_order), " categories total)\n"))
  } else {
    warning("Could not re-level 'cell_type'. No matching cell types found.")
  }

  # Extract and Set Color Map from cell_type_meta_df
  cell_type_colors <- NULL
  
  # 1. Create color dictionary based on the final order
  meta_df_standardized <- cell_type_meta_df %>%
      dplyr::select(cell_type, ID, color) %>%
      dplyr::distinct() %>%
      dplyr::filter(cell_type %in% final_order) %>%
      # Ensure the dataframe is in the order of final_order
      dplyr::mutate(cell_type = factor(cell_type, levels = final_order, ordered = TRUE)) %>%
      dplyr::arrange(cell_type)
      
  colors_for_plot <- meta_df_standardized$color
  
  if (length(colors_for_plot) > 0 && !any(is.na(colors_for_plot))) {
      # Set the color palette for full cell type names
      seurat_object@misc$cell_type_colors <- list(categories = final_order, colors = colors_for_plot)
      cell_type_colors <- setNames(colors_for_plot, final_order)
      cat("Set 'cell_type' color map in seurat_object@misc based on provided metadata colors.\n")
  } else {
      warning("Could not set 'cell_type_colors' from metadata. Default colors will be used.")
  }

  # 2. Map shortened IDs to a new column 'ID'
  old_metadata <- seurat_object[[]]
  new_metadata <- old_metadata %>% left_join(meta_df_standardized, by = join_by(cell_type))
  rownames(new_metadata) <- Cells(seurat_object)
  seurat_object[[]] <- new_metadata
  
  # Filter the order of IDs based on the final order of cell types
  id_map <- setNames(meta_df_standardized$ID, meta_df_standardized$cell_type)
  id_order <- id_map[final_order]
  if (length(id_order) > 0) {
      seurat_object$ID <- factor(seurat_object$ID, levels = id_order, ordered = TRUE)
      
      # Also set the color map for the ID column for plotting
      seurat_object@misc$cell_type_ID_colors <- list(categories = id_order, colors = colors_for_plot)
      
      cat("Shortened cell type labels ('ID') and corresponding colors added to Seurat object metadata.\n")
      
      # Call the new custom plotting function
      plot_umap_with_annotations(
          seurat_obj = seurat_object, 
          group_key = 'ID', 
          color_key = 'cell_type_ID_colors', 
          out_dir = out_dir, 
          out_name = out_name,
          resolution = leiden_res
      )
  } else {
      warning("Could not map 'cell_type' to 'ID' from metadata. Skipping 'cell_type_ID' creation.")
  }

  # Plot regular UMAP Feature plot with cell type overlay (using full name)
  # This plot remains for completeness alongside the custom annotated one
  seurat_object <- SetIdent(seurat_object, value = "cell_type")
  gg <- DimPlot(seurat_object, reduction = "umap", group.by="cell_type", 
                cols=colors_for_plot, # Use custom colors if available
                label = FALSE, pt.size = 0.4) + NoLegend()
  ggsave(gg, filename=file.path(out_dir, paste0(out_name, "_umap_cell_type_no_labels.png")), 
         device="png", dpi=300, width=11)
  
  # --- 4. Filter for Top Cell Labels per Cluster (for plotting) ---
  top_freq_df <- freq_df %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::slice_max(order_by = Frequency, n = 5) %>%
    dplyr::ungroup()
    
  write.csv(top_freq_df, file.path(out_dir, paste0(out_name, "_top_freq_scina_celltypes_per_cluster.csv")), row.names = FALSE)
    
  # --- 5. Plot Cluster Frequency Bar Chart ---
  freq_plot_file <- file.path(out_dir, paste0(out_name, "_scina_cluster_frequency_plot.pdf"))
  cat(paste("Saving cluster frequency plot to:", freq_plot_file, "\n"))
  
  pdf(freq_plot_file, width = 15, height = 12)
  cluster_freq_plot <- ggplot2::ggplot(top_freq_df, ggplot2::aes(x = str_replace_all(SCINA_predicted_celltype, "_", " "), y = Frequency, fill = str_replace_all(SCINA_predicted_celltype, "_", " "))) +
    ggplot2::geom_bar(stat = "identity", position = "dodge", color = "black") +
    ggplot2::facet_wrap(~seurat_clusters, scales = "free_x") +
    ggplot2::labs(
      title = paste("Top 5 SCINA Cell Type Frequencies per Seurat Cluster:", out_name),
      x = "SCINA Predicted Cell Type",
      y = "Relative Frequency (Proportion)",
      fill = "Cell Type"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.title = ggplot2::element_text(face = "bold"))
  print(cluster_freq_plot)
  dev.off() 

  cat("Cluster frequency plot saved.\n")
  
  # --- 6. FIND MARKERS, EXPORT, AND GENERATE PLOTS ---
  cat("Finding markers for SCINA consensus labels...\n")
  
  seurat_object <- SetIdent(seurat_object, value = "cell_type")

  markers <- FindAllMarkers(
    seurat_object,
    only.pos = TRUE,
    min.pct = 0.1,
    logfc.threshold = 0.25,
    assay = DefaultAssay(seurat_object) 
  )
  
  # --- 6.1 EXCEL EXPORT (Top N Markers) ---
  top_excel_markers <- markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(order_by = avg_log2FC, n = n_markers)

  cat(paste("Saving top", n_markers, "markers to Excel file...\n"))
  top_markers_df_list <- split(as.data.frame(top_excel_markers), f=top_excel_markers$cluster)
    
  marker_excel_file <- file.path(out_dir, paste0(out_name, "_scina_cell_type_top_markers.xlsx"))
    
  writexl::write_xlsx(
    top_markers_df_list, 
    path = marker_excel_file,
    col_names = TRUE,
    format_headers = TRUE,
    use_zip64 = FALSE
  )
  cat(paste("Top markers saved to:", marker_excel_file, "\n"))
  
  # --- 6.2 HEATMAP PLOT (Using Top heatmap_n Markers) ---
  top_heatmap_markers <- markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(order_by = avg_log2FC, n = heatmap_n)
    
  heatmap_genes <- unique(top_heatmap_markers$gene)

  if (length(heatmap_genes) == 0) {
    warning("No significant positive markers found for SCINA consensus clusters. Skipping heatmap plot.")
  } else {
    heatmap_plot_file <- file.path(out_dir, paste0(out_name, "_scina_cell_type_heatmap_top", heatmap_n, "_markers.png"))
    cat(paste("Saving marker heatmap (Top", heatmap_n, "markers) to:", heatmap_plot_file, "\n"))
    
    gg_heatmap <- DoHeatmap(seurat_object, features = heatmap_genes) + 
      NoLegend() +
      ggplot2::labs(title = paste("Top", heatmap_n, "Markers Expression Across SCINA Consensus Labels"))
    
    ggplot2::ggsave(
      gg_heatmap, 
      filename = heatmap_plot_file, 
      device = "png", 
      dpi = 300
    )
    cat("Marker heatmap saved.\n")
  }

  # --- 6.3 STACKED VIOLIN PLOT (Using Top violin_n Markers) ---
  top_violin_markers <- markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(order_by = avg_log2FC, n = violin_n)
    
  plot_genes <- unique(top_violin_markers$gene) 
  
  if (length(plot_genes) == 0) {
    warning("No significant positive markers found. Skipping stacked violin plot.")
  } else {
    vln_plot_file <- file.path(out_dir, paste0(out_name, "_scina_cell_type_stacked_violin_top", violin_n, "_markers.pdf"))
    cat(paste("Saving stacked violin plot for top", violin_n, "markers to:", vln_plot_file, "\n"))

    pdf(vln_plot_file, width = 15, height = 10)
    
    # Use custom colors if set, otherwise Seurat's default
    plot_colors <- if (!is.null(cell_type_colors)) { colors_for_plot } else { Seurat:::DiscretePalette(length(unique(seurat_object$cell_type))) }

    stacked_vln_plot <- VlnPlot(
      seurat_object,
      features = plot_genes,
      group.by = "cell_type",
      stack = TRUE,
      flip = TRUE,
      cols = plot_colors
    ) +
      ggplot2::labs(title = paste("Top", violin_n, "Markers Expression by SCINA Consensus Label:", out_name)) +
      ggplot2::theme(legend.position = "none") 
    
    print(stacked_vln_plot)
    dev.off()
    cat("Stacked violin plot saved.\n")
  }
  
  # --- 6.4 Consensus Markers ---
  top_umap_markers <- markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(order_by = avg_log2FC, n = umap_marker_n)
  
  umap_markers_list <- split(top_umap_markers$gene, top_umap_markers$cluster)
  
  # --- 7. Return the Seurat Object and the UMAP Markers List ---
  return(list(seurat_object = seurat_object, umap_markers_list = umap_markers_list))
}



#' Generate Generic UMAP Plots for Seurat Identifiers
#'
#' A flexible wrapper to generate UMAP plots for various metadata columns or gene features.
#' Handles automatic coloring and file saving.
#'
#' @param seurat_obj A Seurat object.
#' @param out_dir String. Output directory.
#' @param out_name String. Base output filename.
#' @param idents_list Vector of strings. Metadata column names or genes to plot.
#' @param labels Logical. Whether to print text labels on the clusters.
#' @param resolution Numeric. Resolution used for 'seurat_clusters' if present.
#' @param file_type String. Output device (e.g., "pdf", "png").
#' @return NULL. Saves files to disk.
umap_seurat <- function(seurat_obj, out_dir, out_name, idents_list, labels = TRUE, resolution = 0.4, file_type = "pdf") {

  # colours for DimPlot
  colours <- c(brewer.pal(n=9,"Set1"), brewer.pal(n=8,"Set2"), 
               brewer.pal(n=12,"Set3"), brewer.pal(n=8,"Dark2"))
  colours <- colours[-13]

  for (ident in idents_list) {
    
    # Check for custom colors in Seurat object's misc slot
    color_map <- NULL
    if (paste0(ident, "_colors") %in% names(seurat_obj@misc)) {
        colors_vec <- seurat_obj@misc[[paste0(ident, "_colors")]]$colors
        categories <- seurat_obj@misc[[paste0(ident, "_colors")]]$categories
        if (length(categories) == length(colors_vec)) {
            color_map <- setNames(colors_vec, categories)
        }
    }
    
    # Fallback to default brewer colors if no custom colors are found
    if (is.null(color_map) && (is.character(seurat_obj[[ident]][,]) || is.factor(seurat_obj[[ident]][,]))) {
        n_unique <- nrow(unique(seurat_obj[[ident]]))
        if (n_unique > length(colours)) {
            warning(paste0("Not enough default colors for ", ident, ". Using a wider default palette."))
            # Use a function to get more colors if necessary
            colour_func <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))
            color_map <- colour_func(n_unique)
        } else {
            color_map <- colours[1:n_unique]
        }
    }

    if (labels == TRUE) {    
      if (ident == "seurat_clusters") {
        resolution <- resolution * 100
        gg <- DimPlot(seurat_obj, reduction = "umap", group.by="seurat_clusters", 
                      cols=color_map,
                      label = TRUE, pt.size = 0.4) + NoLegend()
        ggsave(gg, filename=file.path(out_dir, paste0(out_name, "_umap_res_0_", resolution, ".", 
                                  file_type)), device=file_type, dpi=300)
      } else { 
        if (is.character(seurat_obj[[ident]][,]) || is.factor(seurat_obj[[ident]][,])) {
          gg <- DimPlot(seurat_obj, reduction = "umap", group.by=ident, 
                        cols=color_map,
                        label = TRUE, pt.size = 0.4) + NoLegend()
          ggsave(gg, filename=file.path(out_dir, paste0(out_name, "_umap_", ident, ".", 
                                    file_type)), device=file_type, dpi=300)
        } else {
          gg <- FeaturePlot(seurat_obj, reduction = "umap", features=ident,
                            label = FALSE, pt.size = 0.4)
          ggsave(gg, filename=file.path(out_dir, paste0(out_name, "_umap_", ident, ".", 
                                    file_type)), device=file_type, dpi=300)
        }
      }
    } else if (labels == FALSE) {
      if (ident == "seurat_clusters") {
        resolution <- resolution * 100
        gg <- DimPlot(seurat_obj, reduction = "umap", group.by="seurat_clusters",
                      cols=color_map,
                      label = FALSE, pt.size = 0.4)
        ggsave(gg, filename=file.path(out_dir, paste0(out_name, "_umap_res_0_", resolution, 
                                  "_no_labels.", file_type)), device=file_type, dpi=300, width=11)
      } else { 
        if (is.character(seurat_obj[[ident]][,]) || is.factor(seurat_obj[[ident]][,]))  {
          gg <- DimPlot(seurat_obj, reduction = "umap", group.by=ident,
                        cols=color_map,
                        label = FALSE, pt.size = 0.4)
          ggsave(gg, filename=file.path(out_dir, paste0(out_name, "_umap_", ident, "_no_labels.", 
                                    file_type)), device=file_type, dpi=300, width=11)
        } else {
          gg <- FeaturePlot(seurat_obj, reduction = "umap", features=ident,
                            label = FALSE, pt.size = 0.4)
          ggsave(gg, filename=file.path(out_dir, paste0(out_name, "_umap_", ident, "_no_labels.", 
                                    file_type)), device=file_type, dpi=300, width=11)
        }
      }
    }
  }
}



### Analysis
required_vars <- c("seurat_obj_file", "cell_type_markers_csv", "cell_type_metadata_csv", "leiden_res", "out_dir", "out_name")
for (v in required_vars) {
  if (!exists(v)) stop(paste("Missing required argument:", v))
}

leiden_res <- as.numeric(leiden_res)

# 1. Run the annotation function
seurat_obj <- readRDS(seurat_obj_file)
cell_type_markers_df <- read.csv(cell_type_markers_csv, header = TRUE) 
cell_type_meta_df <- read.csv(cell_type_metadata_csv, header = TRUE) 


scina_results_list <- seurat_scina_annotate(
    seurat_object = seurat_obj,
    markers_df = cell_type_markers_df,
    out_dir = out_dir,
    out_name = out_name
    )

# Unpack results
scina_results_df <- scina_results_list$results_df
markers_list <- scina_results_list$markers_list

# 2. Plot feature plots for original cell type markers
plot_umap_feature_plots(
    seurat_object = seurat_obj, 
    markers_list = markers_list, 
    out_dir = out_dir, 
    out_name = out_name
    )

# 3. Integrate Results back into Seurat Object and Plot Cluster Frequencies
scina_integration_list <- integrate_scina_results(
    seurat_object = seurat_obj,
    scina_results_df = scina_results_df,
    cell_type_meta_df = cell_type_meta_df,
    out_dir = out_dir,
    out_name = out_name,
    n_markers = 200,        # For Excel
    heatmap_n = 5,          # For Heatmap plot
    violin_n = 3,           # For Stacked Violin plot
    umap_marker_n = 6,      # For UMAP FeaturePlot plot based on consensus labels
    leiden_res = leiden_res # Pass resolution for file naming in custom plot
    )

# Unpack results
seurat_obj_anno <- scina_integration_list$seurat_object
umap_markers_list <- scina_integration_list$umap_markers_list

# 4. Plot UMAP feature plots for annotated cell type markers
if (length(umap_markers_list) > 0) {
    # Call the dedicated function with the new list of consensus markers
    plot_umap_feature_plots(
      seurat_object = seurat_obj_anno, 
      markers_list = umap_markers_list, 
      out_dir = out_dir, 
      out_name = paste0(out_name, "_scina_cell_type_top_markers")
    )
} else {
    warning("No significant positive markers found for consensus clusters. Skipping UMAP FeaturePlot generation based on consensus labels.")
}

# 6. Plot umap with cell cycle marker expression
if("Cycling" %in% colnames(cell_type_markers_df)) {
    all_genes_in_seurat = unique(rownames(seurat_obj_anno))
    cycling_markers = cell_type_markers_df$Cycling[cell_type_markers_df$Cycling %in% all_genes_in_seurat]

    umap_seurat(
        seurat_obj = seurat_obj_anno, 
        out_dir = out_dir, 
        out_name = out_name,
        idents_list = cycling_markers, 
        labels = FALSE, 
        file_type = "png"
        )
} else {
    warning("'Cycling' key not found in markers_dict. Skipping cell cycle marker UMAP plot.")
}

# 7. Plot UMAP plots with metadata overlays
idents_list = colnames(seurat_obj_anno@meta.data)
idents_list = idents_list[!idents_list %in% c("barcode", "cell_id", "cell_type_ID")] # Exclude SCINA raw and custom ID plot

umap_seurat(
    seurat_obj = seurat_obj_anno, 
    out_dir = out_dir, 
    out_name = out_name, 
    resolution = leiden_res, 
    idents_list = idents_list, 
    labels = FALSE, 
    file_type = "png"
    )

# 8. Generate blue-print csv file for manual cell type annotation
annotation_df <- seurat_obj_anno@meta.data[, c("seurat_clusters", "cell_type")]
annotation_df <- annotation_df[!duplicated(annotation_df), ]
colnames(annotation_df) <- c("cluster", "cell_type")
annotation_df$cluster <- as.numeric(as.character(annotation_df$cluster))
annotation_df <- annotation_df[order(annotation_df$cluster), ]
annotation_df$order <- as.numeric(annotation_df$cell_type)
write.csv(annotation_df, file.path(out_dir, paste0(out_name, "_annotation.csv")), row.names = FALSE)

# 9. Export metadata as csv
write.csv(seurat_obj_anno@meta.data, file.path(out_dir, paste0(out_name, "_metadata.csv")), row.names = TRUE)

# 10. Save final object
print(head(seurat_obj_anno@meta.data))
saveRDS(seurat_obj_anno, file.path(out_dir, paste0(out_name, "_scRNAseq_no_doublets_annotated.rds")))

cat("Annotation and cluster frequency analysis complete. Results stored in metadata column 'SCINA_predicted_celltype' and saved to disk.\n")