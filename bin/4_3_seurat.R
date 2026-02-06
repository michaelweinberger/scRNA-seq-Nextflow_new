#!/usr/bin/env Rscript


### Run this script to generate Seurat object, remove doublets + identify cluster markers



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
library(patchwork)
library(Seurat)
library(Matrix)
library(sctransform)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(dplyr)
library(cowplot)
library(writexl)
library(DoubletFinder)
library(harmony)

set.seed(42)

# Increase the global size limit to ~10 GiB (value is in bytes)
options(future.globals.maxSize = 15000 * 1024^2)



### Functions

#' Extract Gene Names from GTF
#'
#' Reads a GTF file, filters for gene features, and extracts unique "gene_name" attributes.
#' This is used to ensure the Seurat object only contains valid/annotated genes.
#'
#' @param gtf_file Path to the .gtf file.
#' @return A character vector of unique gene names.
extract_genes_from_gtf <- function(gtf_file) {
  message(paste0("Reading GTF file: ", gtf_file))
  
  # Read GTF (header=FALSE, comment.char="#" handles standard GTF)
  gtf <- read.table(gtf_file, sep = "\t", header = FALSE, comment.char = "#", quote = "")
  
  # Filter for rows representing genes (column 3)
  gtf_genes <- gtf[gtf$V3 == "gene", ]
  
  # Extract gene_name attribute using regex
  # Assumes format: ... gene_name "NAME"; ...
  attributes <- gtf_genes$V9
  gene_names <- regmatches(attributes, regexec('gene_name "([^"]+)"', attributes))
  
  # Extract the capture group (index 2)
  gene_names_clean <- sapply(gene_names, function(x) {
    if (length(x) >= 2) return(x[2]) else return(NA)
  })
  
  # Remove NAs and duplicates
  valid_genes <- unique(gene_names_clean[!is.na(gene_names_clean)])
  
  message(paste0("  > Found ", length(valid_genes), " unique gene names in GTF."))
  return(valid_genes)
}



#' Create Seurat Object (Directory Input)
#'
#' Loads 10X data from a directory, handles multi-modal data (selecting Gene Expression),
#' attaches metadata if provided, and filters genes against a GTF whitelist.
#'
#' @param data_dir Directory containing Cell Ranger matrix files.
#' @param out_dir Output directory path.
#' @param project_name Name for the Seurat project.
#' @param min_cells Minimum cells a gene must be expressed in (default: 3).
#' @param min_features Minimum features per cell (default: 200).
#' @param metadata_file Path to external cell metadata (optional).
#' @param gtf_file Path to GTF for gene filtering (optional).
#' @return A Seurat object.
create_seurat <- function(data_dir, 
                          out_dir, 
                          project_name,
                          min_cells=3,
                          min_features=200, 
                          metadata_file=NULL,
                          gtf_file=NULL
                          ) {
    
    message("Creating data set from ", data_dir)
    
    # read in the 10X data
    seurat_data <- Read10X(data.dir=data_dir)
    
    # Handle multi-modal data (GEX + Antibody/CRISPR)
    if (is.list(seurat_data)) {
      message("\n>>> Multiple feature types detected. Selecting 'Gene Expression'...")
      seurat_data <- seurat_data[["Gene Expression"]]
    }
  
    if (!is.null(metadata_file)) {
        # read in metadata created using cellranger bash script
        meta_data <- read.table(metadata_file, header=TRUE)
        rownames(meta_data) <- meta_data$barcode
        meta_data$barcode <- NULL
        #print(head(meta_data))
    
        # create a new Seurat object
        seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                         project = project_name, 
                                         min.cells = min_cells, 
                                         min.features = min_features, 
                                         meta.data=meta_data
                                         )    
    } else {
        # create a new Seurat object
        seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                         project = project_name, 
                                         min.cells = min_cells, 
                                         min.features = min_features
                                         )
    }
    
    # add project name to cell barcodes
    #seurat_obj <- RenameCells(seurat_obj, new.names = paste(names(Idents(seurat_obj)), 
    #                                                        "_", project_name, sep=""))
    
    # Filter genes based on GTF
    if (!is.null(gtf_file)) {
      message("Filtering Seurat object using GTF file...")
      valid_genes <- extract_genes_from_gtf(gtf_file)
      
      # Intersect GTF genes with Seurat object genes
      genes_to_keep <- intersect(rownames(seurat_obj), valid_genes)
      message(paste0("  > Keeping ", length(genes_to_keep), " genes present in both Seurat object and GTF."))
      
      seurat_obj <- subset(seurat_obj, features = genes_to_keep)
    }

    return(seurat_obj)
}



#' Create Merged Seurat Object (H5 Input)
#'
#' Reads multiple H5 files, assigns batch and sample IDs, merges raw counts, 
#' and joins external metadata via the 'sample_id' column using dplyr::left_join.
#'
#' @param h5_file_paths Vector of paths to .h5 files.
#' @param project_name Name of the merged project.
#' @param out_dir Output directory.
#' @param min_cells Minimum cells per gene.
#' @param min_features Minimum features per cell.
#' @param metadata_file External metadata file path (must contain 'sample_id').
#' @param gtf_file GTF file for gene filtering.
#' @return A merged Seurat object.
create_seurat_h5 <- function(h5_file_paths, 
                             project_name,
                             out_dir,
                             min_cells=3,
                             min_features=200,
                             metadata_file=NULL,
                             gtf_file=NULL) {
    
    seurat_list <- list()
    
    # Process each H5 file individually (minimal initial filtering)
    for (i in seq_along(h5_file_paths)) {
        h5_path <- h5_file_paths[i]
        
        # Define IDs
        current_batch_id <- as.character(i)
        base_name <- basename(h5_path)
        current_sample_id <- sub("\\.h5$", "", base_name) 
        
        message(paste0("\n--- Reading Sample ", i, " ---"))
        message(paste("  > Batch ID:", current_batch_id, "| Sample ID:", current_sample_id))
        
        if (!file.exists(h5_path)) {
            warning(paste("File not found, skipping:", h5_path))
            next
        }
        
        # Read data and handle multiple assays
        data_list <- Read10X_h5(filename = h5_path)
        count_matrix <- if (is.list(data_list)) data_list[[1]] else data_list
        
        # Manually replace underscores with dashes to ensure consistency across all batches
        if (!is.null(rownames(count_matrix))) {
            #print(head(rownames(count_matrix)))
            rownames(count_matrix) <- gsub("_", "-", rownames(count_matrix))
        }

        # Create individual Seurat Object with MINIMAL filtering
        seurat_obj <- CreateSeuratObject(
            counts = count_matrix, 
            min.cells = min_cells, 
            min.features = min_features,
            project = current_batch_id
        )
        
        # Add internal metadata IDs
        seurat_obj$batch_id <- current_batch_id
        seurat_obj$sample_id <- current_sample_id
        
        seurat_list[[i]] <- seurat_obj
    }
    
    if (length(seurat_list) == 0) {
        stop("No valid Seurat objects were created. Check input paths.")
    }
    
    # Merge the Seurat Objects
    message("\n--- Merging All Raw Count Seurat Objects ---")
    
    merged_seurat_obj <- merge(
        x = seurat_list[[1]], 
        y = seurat_list[-1], 
        # Add the numerical batch ID to the start of each cell's barcode for uniqueness
        add.cell.ids = sapply(seurat_list, function(obj) unique(obj$batch_id)),
        project = project_name
    )
    
    message("\n--- Joining Layers ---")
    merged_seurat_obj <- JoinLayers(merged_seurat_obj)

    message(paste0("  > Successfully merged raw data into a single object with ", 
                   ncol(merged_seurat_obj), " total cells."))

    # START EXTERNAL METADATA MERGE
    if (!is.null(metadata_file)) {
        message("\n--- Merging External Metadata using dplyr::left_join ---")

        # read in metadata created using cellranger bash script
        meta_data <- read.table(metadata_file, header=TRUE)

        # Check for 'sample_id' in external metadata
        if (!"sample_id" %in% colnames(meta_data)) {
            stop("External 'metadata' dataframe must contain a column named 'sample_id' for joining.")
        }
    
        # Extract the cell-level metadata from the Seurat object
        seurat_meta <- merged_seurat_obj@meta.data
    
        # Explicitly join the Seurat metadata (cell-level) with the external metadata (sample-level)
        # The 'by = "sample_id"' ensures the external data is broadcast to all matching cells.
        # The 'row.names = FALSE' argument ensures the row names (cell barcodes) are retained.
        merged_meta_data <- dplyr::left_join(
            x = seurat_meta, 
            y = meta_data, 
            by = "sample_id"
        )
    
        # Restore the cell barcode row names to the merged metadata dataframe
        rownames(merged_meta_data) <- rownames(seurat_meta)
    
        # Update the Seurat object metadata
        merged_seurat_obj@meta.data <- merged_meta_data

        # Check for missing external metadata
        new_cols <- setdiff(colnames(meta_data), colnames(seurat_meta))
        if (length(new_cols) > 0) {
            if (any(is.na(merged_seurat_obj@meta.data[, new_cols]))) {
                warning("Some cells have missing external metadata. Check for unmatching 'sample_id' values.")
            }
            message(paste("  > External metadata merged. New columns:", paste(new_cols, collapse = ", ")))
        } else {
             message("  > No new columns were added. External metadata may contain redundant columns.")
        }
    }

    # Filter genes based on GTF
    if (!is.null(gtf_file)) {
      message("Filtering Seurat object using GTF file...")
      valid_genes <- extract_genes_from_gtf(gtf_file)
      valid_genes <- gsub("_", "-", valid_genes)
      
      # Intersect GTF genes with Seurat object genes
      genes_to_keep <- intersect(rownames(merged_seurat_obj), valid_genes)
      message(paste0("  > Keeping ", length(genes_to_keep), " genes present in both Seurat object and GTF."))
      
      merged_seurat_obj <- subset(merged_seurat_obj, features = genes_to_keep)
    }

    message("\nMerged Seurat object created.")

    return(merged_seurat_obj)
}



#' Identify Outliers using MAD
#'
#' Calculates the median and Median Absolute Deviation (MAD) of a metric.
#' Returns a boolean vector indicating values falling outside median +/- nmads * MAD.
#'
#' @param df Dataframe containing the metric.
#' @param metric Column name of the metric to test.
#' @param nmads Number of MADs to use as the threshold.
#' @return Boolean vector (TRUE if outlier).
is_outlier <- function(df, metric, nmads) {
  # Extract the metric column
  M <- df[, metric]
  
  # Calculate median and MAD
  m <- median(M, na.rm = TRUE)
  m_dev <- mad(M, constant = 1, na.rm = TRUE)
  if (m_dev == 0) m_dev <- 0.001
  
  lower_threshold <- m - nmads * m_dev
  upper_threshold <- m + nmads * m_dev
  outlier <- (M < lower_threshold) | (M > upper_threshold)
  
  cat(paste0(metric, ": Lower bound of included range is ", lower_threshold, 
             ", upper bound is ", upper_threshold, ".\n"))
  
  return(outlier)
}



#' QC and Preprocessing Wrapper
#'
#' 1. Calculates mitochondrial percentage.
#' 2. Generates QC Violin and Scatter plots.
#' 3. Filters cells based on outliers (using is_outlier) and %MT cutoff.
#' 4. Performs Cell Cycle Scoring.
#' 5. Normalizes data (LogNormalize or SCTransform).
#' 6. Finds Variable Features (can exclude ribosomal or custom genes).
#' 7. Runs PCA and generates Elbow/JackStraw plots.
#'
#' @param seurat_obj The Seurat object.
#' @param s_genes Vector of S-phase genes.
#' @param g2m_genes Vector of G2M-phase genes.
#' @param assay Normalization method ("RNA" or "SCT").
#' @param pct_mt_cutoff Maximum allowed mitochondrial percentage.
#' @param drop_ribo_prot If "Yes", removes ribosomal genes from variable features.
#' @param genes_drop Vector of specific genes to remove from variable features.
#' @return Processed Seurat object.
preprocess_seurat <- function(seurat_obj, 
                              project_name,
                              out_dir,
                              s_genes,
                              g2m_genes,
                              assay = "RNA",
                              pct_mt_cutoff = 5,
                              jack_straw = FALSE,
                              drop_ribo_prot = "No",
                              genes_drop = NULL
                              ) {
    
    message("\n--- Starting QC and Preprocessing on Seurat Object ---")
    
    # Calculate QC metrics
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-|^Mt-|^mt-")
    
    # Visualize QC metrics
    gg <- VlnPlot(
        object=seurat_obj, 
        features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        combine = FALSE
    )
    gg <- patchwork::wrap_plots(gg, ncol = 3)
    
    ggsave(
        gg, 
        filename=file.path(out_dir, paste0(project_name, "_merged_qc_violin.pdf")), 
        device="pdf", 
        dpi=300
    )
    
    # Visualize feature scatter plots
    plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    
    pdf(file = file.path(out_dir, paste0(project_name, "_merged_qc_scatter.pdf")))
        print(plot1 + plot2)
    dev.off()
    
    # Subsetting (Filtering) data set
    message(paste0("  > Filtering data"))
    
    seurat_obj$outlier <- (
        is_outlier(seurat_obj[[]], "nCount_RNA", 5)
        | is_outlier(seurat_obj[[]], "nFeature_RNA", 5)
        )
    
    seurat_obj$mt_outlier <- seurat_obj$percent.mt > pct_mt_cutoff
    
    n_outliers = sum(seurat_obj$outlier)
    n_mt_outliers = sum(seurat_obj$mt_outlier)
    print(paste("Removed", n_outliers, "cells due to outlier total counts or number of expressed genes"))
    print(paste("Removed", n_mt_outliers, "cells due to outlier mitochondrial gene expression ratio"))
    
    seurat_obj <- subset(seurat_obj, subset = outlier == FALSE)
    seurat_obj <- subset(seurat_obj, subset = mt_outlier == FALSE)
    
    message(paste("  > Post-filtering size:", ncol(seurat_obj), "cells."))
    
    # score cell cycle
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat_obj <- CellCycleScoring(seurat_obj, 
                                   assay = "RNA",
                                   s.features = s_genes, 
                                   g2m.features = g2m_genes, 
                                   set.ident = FALSE
                                   )
    
    # Preprocess data
    if (assay == "RNA") {
        seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
        seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj), vars.to.regress = c("nCount_RNA", "percent.mt"))
    } else if (assay == "SCT") {
        seurat_obj <- SCTransform(seurat_obj, vars.to.regress = c("nCount_RNA", "percent.mt"), conserve.memory = TRUE)
    } else {
        print("Please set assay to either 'RNA' or 'SCT'.")
    }
  
    # Exclude cell cycle genes from HVG
    current_hvg <- VariableFeatures(seurat_obj)
    cc_genes <- c(s_genes, g2m_genes)

    if (length(cc_genes) > 0) {
        print(paste("Excluding", length(cc_genes), "cell cycle genes from HVG."))
        new_hvg <- setdiff(current_hvg, cc_genes)
        VariableFeatures(seurat_obj) <- new_hvg
    }

    # Exclude ribosomal protein genes from HVG
    if (drop_ribo_prot == "Yes") {
        current_hvg <- VariableFeatures(seurat_obj)
        ribo_genes <- grep(pattern = "^RP[SL]", x = current_hvg, value = TRUE, ignore.case = TRUE)
  
        if (length(ribo_genes) > 0) {
            print(paste("Identified", length(ribo_genes), "ribosomal protein genes. Excluding them from HVG."))
            new_hvg <- setdiff(current_hvg, ribo_genes)
            VariableFeatures(seurat_obj) <- new_hvg
        }
    }
    
    if (!is.null(genes_drop)) {
        current_hvg <- VariableFeatures(seurat_obj)
        print("Excluding custom genes from HVG.")
        new_hvg <- setdiff(current_hvg, genes_drop)
        VariableFeatures(seurat_obj) <- new_hvg
    }
    
    # Identify the 10 most highly variable genes
    top10 <- head(VariableFeatures(seurat_obj), 10)
    plot1 <- VariableFeaturePlot(seurat_obj)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    pdf(file = file.path(out_dir, paste0(project_name, "_qc_variable_features.pdf")))
        print(plot2)
    dev.off()
  
    # compute PCs
    seurat_obj <- RunPCA(seurat_obj, assay = assay, features = VariableFeatures(object = seurat_obj))
    print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)
    
    # Determine number of PCs to be used in downstream analysis
    if (jack_straw == TRUE) {
        message("  > Computing JackStraw score (can be time-consuming)...")
        seurat_obj <- JackStraw(seurat_obj, num.replicate = 100)
        seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:40)
        
        pdf(file = file.path(out_dir, paste0(project_name, "_jack_straw.pdf")))
            print(JackStrawPlot(seurat_obj, dims = 1:40))
        dev.off()
    }
    
    pdf(file = file.path(out_dir, paste0(project_name, "_elbow_plot.pdf")))
        print(ElbowPlot(seurat_obj, ndims = 40, reduction = "pca"))
    dev.off()
    
    return(seurat_obj)
}



#' Add Gene Symbols from BLAST Results
#'
#' Merges BLAST alignment results into the Seurat feature metadata.
#' Useful for cross-species analysis or adding common names to ID-based features.
#'
#' @param seurat_obj Seurat object.
#' @param blast_df Dataframe containing BLAST output (qseqid, sseqid, etc.).
#' @return Seurat object with updated feature metadata.
seurat_add_gene_symbols <- function(seurat_obj, 
                                    blast_df, 
                                    out_dir, 
                                    project_name,
                                    assay="RNA") {
  
  # compute the absolute number of matching basepairs in each row
  blast_df$matches <- blast_df$length - blast_df$mismatch
    
  # create column indicating whether uppercase query gene symbol matches uppercase search gene symbol
  blast_df$name_match_bool <- toupper(blast_df$qseqid) == toupper(blast_df$sseqid)
    
  # sort dataframe + drop duplicated query gene symbol rows
  blast_df <- blast_df[order(blast_df$qseqid, 
                             -blast_df$matches, 
                             blast_df$evalue, 
                             -blast_df$name_match_bool),]
  blast_df <- blast_df[!duplicated(blast_df$qseqid),]
  
  # sort dataframe + drop duplicated target gene symbol rows
  blast_df <- blast_df[order(blast_df$sseqid, 
                             -blast_df$matches, 
                             blast_df$evalue, 
                             -blast_df$name_match_bool),]
  blast_df <- blast_df[!duplicated(blast_df$sseqid),]
  
  # merge Seurat object feature metadata and blast results
  #seurat_obj@assays$RNA@meta.features$qseqid <- rownames(seurat_obj@assays$RNA@meta.features)
  feature_meta <- as.data.frame(rownames(seurat_obj[[assay]]))
  colnames(feature_meta) <- "qseqid"
  rownames(feature_meta) <- feature_meta$qseqid
  seurat_obj[[assay]][["qseqid"]] <- rownames(feature_meta)
  feature_meta <- feature_meta %>% 
                  left_join(blast_df[,c("qseqid", "sseqid")], by = join_by(qseqid))
  rownames(feature_meta) <- seurat_obj[[assay]]@meta.features["qseqid"][,1]
  seurat_obj[[assay]] <- AddMetaData(seurat_obj[[assay]], metadata = feature_meta)
  
  return(seurat_obj)
  
  #keep <- seurat_obj@assays$RNA@meta.features[!is.na(seurat_obj@assays$RNA@meta.features$sseqid), "qseqid"]
  #seurat_obj <- subset(seurat_obj, features = keep)
  #rownames(seurat_obj@assays$RNA@meta.features) <- seurat_obj@assays$RNA@meta.features$sseqid
  #seurat_obj@assays$RNA@counts@Dimnames[[1]] <- seurat_obj@assays$RNA@meta.features$sseqid
  #seurat_obj@assays$RNA@data@Dimnames[[1]] <- seurat_obj@assays$RNA@meta.features$sseqid
}




#' Cluster Cells and Run UMAP
#'
#' 1. If 'harmony_var' is provided, runs Harmony integration to correct batch effects.
#' 2. Constructs the Nearest Neighbor graph (FindNeighbors).
#' 3. Clusters cells (FindClusters).
#' 4. Runs UMAP for visualization.
#'
#' @param seurat_obj Seurat object.
#' @param n_dims Number of PCs/Harmony dimensions to use.
#' @param harmony_var Metadata column for batch correction (pass NULL to skip Harmony).
#' @param resolution Clustering resolution (higher = more clusters).
#' @return Seurat object with cluster IDs and UMAP reduction.
cluster_seurat <- function(seurat_obj, assay, out_dir, project_name, n_dims=40, harmony_var=NULL, 
                           resolution=0.4) {
  
  print("Clustering cells")
  
  if (!is.null(harmony_var)) {
    print("Performing harmonisation")
  
    pdf(file = file.path(out_dir, paste0(project_name, "_harmony_run.pdf")))
      options(repr.plot.height = 2.5, repr.plot.width = 6)
      seurat_obj <- seurat_obj %>% RunHarmony(harmony_var, assay.use = assay, plot_convergence = TRUE, max_iter=50)
    dev.off()
  
    # compare PCs before and after harmonisation
    png(file = file.path(out_dir, paste0(project_name, "_PC_without_harmony.png")))
      options(repr.plot.height = 5, repr.plot.width = 12)
      p1 <- DimPlot(object = seurat_obj, reduction = "pca", pt.size = .1, group.by = harmony_var, raster = FALSE)
      p2 <- VlnPlot(object = seurat_obj, features = "PC_1", group.by = harmony_var, pt.size = .1, raster = FALSE)
      print(cowplot::plot_grid(p1,p2))
    dev.off()
    
    png(file = file.path(out_dir, paste0(project_name, "_PC_with_harmony.png")))
      options(repr.plot.height = 5, repr.plot.width = 12)
      p1 <- DimPlot(object = seurat_obj, reduction = "harmony", pt.size = .1, group.by = harmony_var, raster = FALSE)
      p2 <- VlnPlot(object = seurat_obj, features = "harmony_1", group.by = harmony_var, pt.size = .1, raster = FALSE)
      print(cowplot::plot_grid(p1,p2))
    dev.off()
  
    # plot elbow plot again
    pdf(file = file.path(out_dir, paste0(project_name, "_elbow_plot_harmony.pdf")))
      print(ElbowPlot(seurat_obj, ndims = 40, reduction = "pca"))
    dev.off()
    
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:n_dims)
    seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
     
    # Run non-linear dimensional reduction (UMAP/tSNE)
    # If you haven't installed UMAP, you can do so via reticulate::py_install(packages = 'umap-learn')
    seurat_obj <- RunUMAP(seurat_obj, assay = assay, reduction = "harmony", dims = 1:n_dims)
    
  } else {
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:n_dims)
    seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
    
    # Run non-linear dimensional reduction (UMAP/tSNE)
    # If you haven't installed UMAP, you can do so via reticulate::py_install(packages = 'umap-learn')
    seurat_obj <- RunUMAP(seurat_obj, assay = assay, dims = 1:n_dims)
  }
  
  return(seurat_obj)
}



#' Filter Samples by Cell Count
#'
#' plots cell counts per sample (based on harmony_var) and removes samples 
#' that have fewer cells than the specified threshold.
#'
#' @param seurat_obj Seurat object.
#' @param harmony_var Metadata column used to group samples (e.g., "sample_id").
#' @param threshold Minimum number of cells required to keep a sample.
#' @return Filtered Seurat object.
filter_samples_by_count <- function(seurat_obj, harmony_var, out_dir, project_name, threshold = 50) {
  
  if (is.null(harmony_var) || !(harmony_var %in% colnames(seurat_obj@meta.data))) {
    message("Warning: harmony_var not found in metadata. Skipping sample filtering.")
    return(seurat_obj)
  }
  
  # 1. Calculate counts per sample
  # Using table() and converting to dataframe for ggplot
  counts_df <- as.data.frame(table(seurat_obj@meta.data[[harmony_var]]))
  colnames(counts_df) <- c("SampleID", "Count")
  counts_df <- counts_df[order(counts_df$Count), ] # Sort for the plot
  
  # 2. Plotting (Horizontal Bar Plot)
  gg <- ggplot(counts_df, aes(x = reorder(SampleID, Count), y = Count)) +
    geom_bar(stat = "identity", fill = "skyblue", color = "black") +
    coord_flip() + # Makes it horizontal
    geom_hline(yintercept = threshold, color = "red", linetype = "dashed", size = 1) +
    annotate("text", x = 1, y = threshold, label = paste("Threshold:", threshold), 
             color = "red", vjust = -1, hjust = -0.1) +
    labs(title = "Number of Cells per Sample", 
         x = "Sample ID", 
         y = "Cell Count") +
    theme_minimal()
  
  ggsave(plot = gg, filename = file.path(out_dir, paste0(project_name, "_cell_counts_horizontal.pdf")), 
         width = 12, height = 8, device = "pdf")
  
  # 3. Filter the Seurat Object
  samples_to_keep <- counts_df$SampleID[counts_df$Count >= threshold]
  samples_to_remove <- counts_df$SampleID[counts_df$Count < threshold]
  
  message("Removed samples: ", paste(samples_to_remove, collapse = ", "))
  
  # Subset Seurat object to keep only cells from valid samples
  # Using the binary %in% operator on the metadata column
  cells_to_keep <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[harmony_var]] %in% samples_to_keep, ])
  seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
  
  # 4. Clean up categorical metadata (droplevels)
  seurat_obj@meta.data[[harmony_var]] <- droplevels(as.factor(seurat_obj@meta.data[[harmony_var]]))
  
  return(seurat_obj)
}



#' Identify Doublets
#'
#' Splits the object by a metadata variable (e.g., sample) and runs DoubletFinder 
#' on each subset. Returns a dataframe of identified doublets for removal.
#'
#' @param seurat_obj Seurat object.
#' @param metadata_var Variable to split by before doublet detection.
#' @param pK The PC neighborhood size.
#' @param dfr Expected doublet formation rate (default 0.03 or 3%).
#' @return Dataframe containing cell IDs identified as doublets.
doublet_finder <- function(seurat_obj, out_dir, project_name, metadata_var, 
                           n_dims=20, resolution=0.4, pK=0.09, dfr=0.03) {
  
  Idents(seurat_obj) <- metadata_var
  
  # initiate dataframe for doublet results
  doublets <- data.frame(cell_id=character(), 
                         metadata_var=character(), 
                         seurat_clusters=character())
  
  loop_vars <- if(is.null(metadata_var)) "All" else unique(seurat_obj[[metadata_var]][[1]])
  
  for (i in loop_vars) {
    print(paste("Identifying doublets in", i))
    
    if (i == "All") {
        seurat_tmp <- seurat_obj
    } else {
        seurat_tmp <- subset(x = seurat_obj, idents = i)
    }

    if ("SCT" %in% Assays(seurat_tmp)) {
        DefaultAssay(object = seurat_tmp) <- "RNA"
        seurat_tmp[["SCT"]] <- NULL
    }

    # Re-process subset here for better accuracy
    seurat_tmp <- NormalizeData(seurat_tmp, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat_tmp <- FindVariableFeatures(seurat_tmp, selection.method = "vst", nfeatures = 3000)
    seurat_tmp <- ScaleData(seurat_tmp, features = rownames(seurat_tmp))

    seurat_tmp <- RunPCA(seurat_tmp, assay = "RNA", features = VariableFeatures(object = seurat_tmp))
    seurat_tmp <- FindNeighbors(seurat_tmp, dims = 1:n_dims)
    seurat_tmp <- FindClusters(seurat_tmp, resolution = resolution)
    
    # pK Identification (no ground-truth)
    #sweep.res.list <- paramSweep(seurat_tmp, PCs = 1:n_dims, sct = FALSE)
    #sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    #bcmvn <- find.pK(sweep.stats)
    
    #pdf(file = paste(out_dir, "/", project_name, "_", i, "_doubletfinder_pk_bar.pdf", sep=""))
    #print(barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2))
    #dev.off()
    
    nExp_poi <- round(dfr*nrow(seurat_tmp[[]]))
    
    # Homotypic Doublet Proportion Estimate 
    #homotypic.prop <- modelHomotypic(seurat_tmp[["seurat_clusters"]])
    #nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    
    # Run DoubletFinder with varying classification stringencies
    seurat_tmp <- doubletFinder(seurat_tmp, PCs = 1:n_dims, pN = 0.25, pK = pK, 
                                nExp = nExp_poi, reuse.pANN = NULL, sct = FALSE)
    
    # name of the DF prediction can change, so extract the correct column name.
    DF.names = colnames(seurat_tmp[[]])[grepl("DF.classification", colnames(seurat_tmp[[]]))]
    DF.name <- tail(DF.names, 1)
    
    pdf(file = file.path(out_dir, paste0(project_name, "_", i, "_doubletfinder_umap.pdf")))
    print(cowplot::plot_grid(ncol = 2, DimPlot(seurat_tmp, group.by = "seurat_clusters") + NoAxes(),
                             DimPlot(seurat_tmp, group.by = DF.name) + NoAxes()))
    dev.off()
    
    # extract IDs of doublets
    cols_to_extract <- c("cell_id", "seurat_clusters")
    if (!is.null(metadata_var)) cols_to_extract <- c(cols_to_extract, metadata_var)
    
    tmp <- as.data.frame(seurat_tmp[[]])
    tmp$cell_id <- rownames(tmp)
    tmp <- tmp[tmp[,DF.name]=="Doublet", cols_to_extract]
    doublets <- rbind(doublets, tmp)
    
    print(paste("Identified", length(seurat_tmp[[DF.name]][seurat_tmp[[DF.name]]=="Doublet"]), 
                "doublets"))
  }
  
  return(doublets)
}



#' Plot UMAP
#'
#' Iterates through a list of metadata columns ('idents_list') and generates 
#' UMAP plots for each. Can save as PDF (vector) or PNG (raster).
#'
#' @param seurat_obj Seurat object.
#' @param idents_list List of metadata columns to overlay on the UMAP.
#' @param labels Boolean; if TRUE, adds text labels to clusters.
#' @param file_type Output format ("pdf" or "png").
#' @return The Seurat object (returns unchanged).
umap_seurat <- function(seurat_obj, out_dir, project_name, idents_list, labels=TRUE, 
                        resolution=0.4, file_type="pdf") {
  
  # colours for DimPlot
  colours <- c(brewer.pal(n=9,"Set1"), brewer.pal(n=8,"Set2"), 
               brewer.pal(n=12,"Set3"), brewer.pal(n=8,"Dark2"))
  colours <- colours[-13]
  
  for (ident in idents_list) {
    
    if (labels == TRUE) {    
      if (ident == "seurat_clusters") {
        resolution <- resolution * 100
        gg <- DimPlot(seurat_obj, reduction = "umap", group.by="seurat_clusters", 
                      cols=colours[1:nrow(unique(seurat_obj[[ident]]))],
                      label = TRUE, pt.size = 0.4) + NoLegend()
        ggsave(gg, filename=file.path(out_dir, paste0(project_name, "_umap_res_0_", resolution, ".", 
                                  file_type)), device=file_type, dpi=300)
      } else { 
        if (is.character(seurat_obj[[ident]][,]) || is.factor(seurat_obj[[ident]][,])) {
          gg <- DimPlot(seurat_obj, reduction = "umap", group.by=ident, 
                        cols=colours[1:nrow(unique(seurat_obj[[ident]]))],
                        label = TRUE, pt.size = 0.4) + NoLegend()
          ggsave(gg, filename=file.path(out_dir, paste0(project_name, "_umap_", ident, ".", 
                                    file_type)), device=file_type, dpi=300)
        } else {
          gg <- FeaturePlot(seurat_obj, reduction = "umap", features=ident,
                            label = FALSE, pt.size = 0.4)
          ggsave(gg, filename=file.path(out_dir, paste0(project_name, "_umap_", ident, ".", 
                                    file_type)), device=file_type, dpi=300)
        }
      }
    } else if (labels == FALSE) {
      if (ident == "seurat_clusters") {
        resolution <- resolution * 100
        gg <- DimPlot(seurat_obj, reduction = "umap", group.by="seurat_clusters",
                      cols=colours[1:nrow(unique(seurat_obj[[ident]]))],
                      label = FALSE, pt.size = 0.4)
        ggsave(gg, filename=file.path(out_dir, paste0(project_name, "_umap_res_0_", resolution, 
                                  "_no_labels.", file_type)), device=file_type, dpi=300, width=11)
      } else { 
        if (is.character(seurat_obj[[ident]][,]) || is.factor(seurat_obj[[ident]][,]))  {
          gg <- DimPlot(seurat_obj, reduction = "umap", group.by=ident,
                        cols=colours[1:nrow(unique(seurat_obj[[ident]]))],
                        label = FALSE, pt.size = 0.4)
          ggsave(gg, filename=file.path(out_dir, paste0(project_name, "_umap_", ident, "_no_labels.", 
                                    file_type)), device=file_type, dpi=300, width=11)
        } else {
          gg <- FeaturePlot(seurat_obj, reduction = "umap", features=ident,
                            label = FALSE, pt.size = 0.4)
          ggsave(gg, filename=file.path(out_dir, paste0(project_name, "_umap_", ident, "_no_labels.", 
                                    file_type)), device=file_type, dpi=300, width=11)
        }
      }
    }
  }
  return(seurat_obj)
}



#' Identify and Plot Cluster Markers
#'
#' 1. Runs FindAllMarkers to identify upregulated genes for each cluster.
#' 2. Saves the top 'n_markers' per cluster to an Excel file.
#' 3. Generates a Heatmap of the top 5 markers per cluster.
#'
#' @param seurat_obj Seurat object.
#' @param idents Identity class to use (default: "seurat_clusters").
#' @param n_markers Number of top markers to save in the Excel report.
#' @return The Seurat object.
markers_seurat <- function(seurat_obj, out_dir, project_name, idents="seurat_clusters", 
                           n_markers=200) {
  
  Idents(seurat_obj) <- idents
  
  seurat_obj.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, 
                                       logfc.threshold = 0.25)
  
  # save the top markers of each cluster
  top_markers <- as.data.frame(seurat_obj.markers %>% group_by(cluster) %>% top_n(n = n_markers, wt = avg_log2FC))
  top_markers_df_list <- split(top_markers, f=top_markers$cluster)
  writexl::write_xlsx(top_markers_df_list, 
                      path=file.path(out_dir, paste0(project_name, "_markers.xlsx")),
                      col_names = TRUE,
                      format_headers = TRUE,
                      use_zip64 = FALSE)
  
  # draw heatmap of top marker expression
  top10 <- seurat_obj.markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5)
  
  gg <- DoHeatmap(seurat_obj, features = top10$gene) + NoLegend()
  ggsave(gg, filename=file.path(out_dir, paste0(project_name, "_heatmap_markers.pdf")), 
         device="pdf", dpi=300, useDingbats=FALSE)
  
  return(seurat_obj)
}



### Analysis
assay <- "SCT"

if (!exists("min_genes")) min_genes <- 200
if (!exists("max_genes")) max_genes <- 8000
if (!exists("max_perc_mt")) max_perc_mt <- 5
if (!exists("min_cells")) min_cells <- 3
if (!exists("n_pcs")) n_pcs <- 40
if (!exists("leiden_res")) leiden_res <- 0.4
if (!exists("drop_ribo_prot")) drop_ribo_prot <- "No"
min_genes <- as.numeric(min_genes)
max_genes <- as.numeric(max_genes)
max_perc_mt <- as.numeric(max_perc_mt)
min_cells <- as.numeric(min_cells)
n_pcs <- as.numeric(n_pcs)
leiden_res <- as.numeric(leiden_res)

if (!exists("blast_file")) stop("Variable 'blast_file' is required but missing.")
if (!exists("metadata_file")) stop("Variable 'metadata_file' is required but missing.")
if (!exists("cellcycle_genes_csv")) stop("Variable 'cellcycle_genes_csv' is required but missing.")
meta_data <- read.table(metadata_file, header=TRUE)
cellcycle_genes <- read.csv(cellcycle_genes_csv, header=TRUE)
s_genes <- cellcycle_genes$G1_S
g2m_genes <- cellcycle_genes$G2_M

if (!exists("gtf_file")) {
    gtf_file = NULL
}

if (!exists("genes_drop_csv") || file.size(genes_drop_csv) == 0) {
    genes_drop = NULL
} else {
    genes_drop <- read.csv(genes_drop_csv, header=FALSE)
    print("Reading in custom genes to drop from HVG.")
}


# ensure missing harmony_var is received as NULL in function cluster_seurat
if (!exists("harmony_var") || is.null(harmony_var) || !harmony_var %in% colnames(meta_data)) {
    harmony_var = NULL
}

# create output directory if it does not exist
if (!dir.exists(out_dir)) {
    dir.create(out_dir)
}

# check if ambient RNA corrected counts are available
if (exists("soupx_dir") && dir.exists(soupx_dir)) {
    data_dir <- soupx_dir
} else {
    if (exists("data_dir") && dir.exists(data_dir)) {
        data_dir <- data_dir
    } else {
      stop("Variable 'data_dir' or 'soupx_dir' must be set.")
    }
}

# check if input directory contains Cell Ranger output
matrix_file_list <- list.files(
  path = data_dir,
  pattern = "matrix\\.mtx\\.gz$", # Regular expression: matches files ending with 'matrix.mtx.gz'
  full.names = FALSE,
  recursive = FALSE
)

if (length(matrix_file_list) > 0) {
  if (!file.exists(file.path(out_dir, paste0(project_name, "_scRNAseq_with_doublets.rds")))) {
    # create Seurat object from Cell Ranger output
    seurat_obj <- create_seurat(data_dir=data_dir,
                                out_dir=out_dir, 
                                project_name=project_name,
                                min_features=min_genes, 
                                min_cells=min_cells,
                                metadata_file=metadata_file,
                                gtf_file=NULL)
    
    seurat_obj <- preprocess_seurat(seurat_obj=seurat_obj,
                                    s_genes=s_genes,
                                    g2m_genes=g2m_genes,
                                    assay=assay, 
                                    out_dir=out_dir, 
                                    project_name=project_name, 
                                    jack_straw=FALSE, 
                                    pct_mt_cutoff=max_perc_mt,
                                    drop_ribo_prot=drop_ribo_prot,
                                    genes_drop=genes_drop
                                    )
    if (file.size(blast_file) == 0) {
      print("Skipping gene symbol addition: Blast file is empty.")
    } else {
        blast_df <- read.table(blast_file, header=TRUE)
        seurat_obj <- seurat_add_gene_symbols(seurat_obj=seurat_obj, 
                                              blast_df=blast_df,
                                              assay=assay,
                                              out_dir=out_dir, 
                                              project_name=project_name
                                              )
    }

    seurat_obj <- filter_samples_by_count(seurat_obj = seurat_obj, 
                                          harmony_var = harmony_var, 
                                          out_dir = out_dir, 
                                          project_name = project_name, 
                                          threshold = 50
                                          )

    seurat_obj <- cluster_seurat(seurat_obj=seurat_obj,
                                 assay=assay,
                                 out_dir=out_dir, 
                                 project_name=project_name, 
                                 n_dims=n_pcs, 
                                 resolution=leiden_res, 
                                 harmony_var=NULL
                                 )

    saveRDS(seurat_obj, file.path(out_dir, paste0(project_name, "_scRNAseq_with_doublets.rds")))
  } else {
    seurat_obj <- readRDS(file.path(out_dir, paste0(project_name, "_scRNAseq_with_doublets.rds")))
  }
} else {
  # check if input directory contains h5 files
  h5_file_list <- list.files(
    path = data_dir,
    pattern = ".h5$", # Regular expression: matches files ending with '.h5'
    full.names = TRUE,
    recursive = FALSE
  )
  if (length(h5_file_list) > 0) {
    if (!file.exists(file.path(out_dir, paste0(project_name, "_scRNAseq_with_doublets.rds")))) {
      # create Seurat object from .h5 files
      seurat_obj <- create_seurat_h5(h5_file_paths=h5_file_list,
                                     out_dir=out_dir, 
                                     project_name=project_name,
                                     min_features=min_genes, 
                                     min_cells=min_cells, 
                                     metadata_file=metadata_file,
                                     gtf_file=NULL)

      seurat_obj <- preprocess_seurat(seurat_obj=seurat_obj,
                                      s_genes=s_genes,
                                      g2m_genes=g2m_genes,
                                      assay=assay, 
                                      out_dir=out_dir, 
                                      project_name=project_name, 
                                      jack_straw=FALSE, 
                                      pct_mt_cutoff=max_perc_mt,
                                      drop_ribo_prot=drop_ribo_prot,
                                      genes_drop=genes_drop
                                      )

      if (file.size(blast_file) == 0) {
        print("Skipping gene symbol addition: Blast file is empty.")
      } else {
        blast_df <- read.table(blast_file, header=TRUE)
        seurat_obj <- seurat_add_gene_symbols(seurat_obj=seurat_obj, 
                                              blast_df=blast_df,
                                              assay=assay,
                                              out_dir=out_dir, 
                                              project_name=project_name
                                              )
      }

      seurat_obj <- filter_samples_by_count(seurat_obj = seurat_obj, 
                                            harmony_var = harmony_var, 
                                            out_dir = out_dir, 
                                            project_name = project_name, 
                                            threshold = 50
                                            )

      seurat_obj <- cluster_seurat(seurat_obj=seurat_obj,
                                   assay=assay,
                                   out_dir=out_dir, 
                                   project_name=project_name, 
                                   n_dims=n_pcs, 
                                   resolution=leiden_res, 
                                   harmony_var=NULL
                                   )

      saveRDS(seurat_obj, file.path(out_dir, paste0(project_name, "_scRNAseq_with_doublets.rds")))
    } else {
      seurat_obj <- readRDS(file.path(out_dir, paste0(project_name, "_scRNAseq_with_doublets.rds")))
    }
  } else {
    stop("Please provide an input directory that contains either \
          a Cell Ranger 'matrix.mtx.gz' file or '.h5' files with mapped data.")
  }
}

print(seurat_obj)

if (!file.exists(file.path(out_dir, paste0(project_name, "_scRNAseq_analysed_no_doublets.rds")))) {

  # remove doublets
  doublets <- doublet_finder(seurat_obj=seurat_obj,
                             out_dir=out_dir, 
                             project_name=project_name, 
                             metadata_var=harmony_var, 
                             n_dims=30, 
                             resolution=leiden_res, 
                             pK=0.09, 
                             dfr=0.03
                             )
  write.csv(doublets, file.path(out_dir, paste0(project_name, "_doublets.csv")))

  seurat_obj[["cell_id"]] <- rownames(seurat_obj[[]])
  seurat_obj <- subset(x = seurat_obj, subset = cell_id %in% doublets$cell_id, invert=TRUE)
  
  # pre-process again after doublet removal
  seurat_obj <- preprocess_seurat(seurat_obj=seurat_obj,
                                  s_genes=s_genes,
                                  g2m_genes=g2m_genes,
                                  assay=assay, 
                                  out_dir=out_dir, 
                                  project_name=project_name, 
                                  jack_straw=FALSE, 
                                  pct_mt_cutoff=max_perc_mt,
                                  drop_ribo_prot=drop_ribo_prot,
                                  genes_drop=genes_drop
                                  )

  # cluster cells
  seurat_obj <- cluster_seurat(seurat_obj=seurat_obj,
                               assay=assay,
                               out_dir=out_dir, 
                               project_name=project_name, 
                               n_dims=n_pcs, 
                               resolution=leiden_res, 
                               harmony_var=harmony_var
                               )

  # For some reason, seurat cluster levels are 0-based, but values are 1-based
  # -> adjust levels to be 1-based
  seurat_obj$seurat_clusters <- as.factor(as.numeric(as.character(seurat_obj$seurat_clusters)) + 1)
  
  # plot UMAPs with metadata overlay
  idents_list = colnames(seurat_obj[[]])
  idents_list = idents_list[idents_list != "barcode"]
  seurat_obj <- umap_seurat(seurat_obj=seurat_obj, 
                            out_dir=out_dir, 
                            project_name=project_name, 
                            resolution=leiden_res, 
                            idents_list=idents_list, 
                            labels=FALSE, 
                            file_type="png"
                            )
  
  # Identify cluster markers
  seurat_obj <- markers_seurat(seurat_obj=seurat_obj, 
                               out_dir=out_dir, 
                               project_name=project_name, 
                               idents="seurat_clusters", 
                               n_markers=200
                               )

  # export metadata as csv
  meta_data_export <- seurat_obj[[]]
  write.csv(meta_data_export, file.path(out_dir, paste0(project_name, "_cell_metadata.csv")))
  write.csv(seurat_obj[[assay]][[]], file.path(out_dir, paste0(project_name, "_feature_metadata.csv")))
  
  # save Seurat object
  saveRDS(seurat_obj, file.path(out_dir, paste0(project_name, "_scRNAseq_analysed_no_doublets.rds")))
}