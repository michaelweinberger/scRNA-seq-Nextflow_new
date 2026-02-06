#!/usr/bin/env python



# -*- coding: utf-8 -*-
"""
Preprocess and cluster scRNA-seq data with scanpy

"""


### Modules
#import matplotlib
#from matplotlib import rcParams
#import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import anndata as ad
import matplotlib.pyplot as plt
import re
import os
import random
import argparse
from scipy.stats import median_abs_deviation

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.autosave = True	          # do not show plots, save them in figdir
sc.settings.set_figure_params(dpi=300, facecolor='white')

random.seed(0)



## unpack arguments imported from bash parent script
parser = argparse.ArgumentParser(description='Preprocess and cluster scRNA-seq data with scanpy')
parser.add_argument('-i','--input_dir', help='Directory containing cellranger output files \
                    (default: /outs/count/filtered_feature_bc_matrix)', required=True)
parser.add_argument('-m','--metadata', help='.tsv file containing "barcode" column with cell barcodes, \
                    "sample" column with sample names and additional metadata columns', required=True)
parser.add_argument('-gtf','--gtf_file', help='Filepath to a gtf file to subset features in adata', required=False)
parser.add_argument('-s','--soupx_dir', help='Path to a directory containing ambient RNA corrected, 10X formatted count matrix, feature + barcode files', required=False)
parser.add_argument('-b','--blast_file', help='Filepath to tab-delimited blast output file in -outfmt 6 format', required=False)
parser.add_argument('-ccg','--cellcycle_genes_csv', help='Filepath of a csv file containing cell cycle S and G2/M phase marker genes', required=True)
parser.add_argument('-o','--output_dir', help='Directory for output files', required=True)
parser.add_argument('-n','--output_name', help='Prefix for output file names', required=True)
parser.add_argument('-d','--drop_ribo_prot', help='Whether to exclude ribosomal protein genes from highly variable genes before clustering', required=True)
parser.add_argument('-g','--genes_drop_csv', help='Filepath of a csv file containing genes to drop from highly variable genes before clustering', required=True)
parser.add_argument('-mig','--min_genes', help='Minimum number of genes detected for a cell to be kept in the dataset', required=True)
parser.add_argument('-mag','--max_genes', help='Maximum number of genes detected for a cell to be kept in the dataset', required=True)
parser.add_argument('-mam','--max_perc_mt', help='Maximum percentage of mitochondrial gene counts detected for a cell to be kept in the dataset', required=True)
parser.add_argument('-mic','--min_cells', help='Minimum number of cells in which a gene needs to be detected for it to be kept in the dataset', required=True)
parser.add_argument('-npc','--n_pcs', help='Number of principal components to use for neighbourhood graph', required=True)
parser.add_argument('-hv','--harmony_var', help='Name of metadata column to use for data integration', required=True)
parser.add_argument('-r','--leiden_res', help='Resolution of Leiden clustering', required=False)

args = parser.parse_args()
in_dir = args.input_dir
metadata = args.metadata
gtf_file = args.gtf_file if args.gtf_file else None
soupx_dir = args.soupx_dir if args.soupx_dir else None
blast_file = args.blast_file
cellcycle_genes_csv = args.cellcycle_genes_csv
out_dir = args.output_dir
out_name = args.output_name
drop_ribo_prot = args.drop_ribo_prot
genes_drop_csv = args.genes_drop_csv
min_genes =  int(args.min_genes)
max_genes = int(args.max_genes)
max_perc_mt = int(args.max_perc_mt)
min_cells = int(args.min_cells)
n_pcs = int(args.n_pcs)
harmony_var = args.harmony_var
leiden_res = float(args.leiden_res) if args.leiden_res else 0.4

# set plot file suffix
sc.settings.plot_suffix = f"_scanpy_{out_name}"
sc.settings.figdir = out_dir



### Functions

def scanpy_create_10X(cellranger_dir, 
                      out_dir, 
                      out_name, 
                      metadata=None, 
                      min_genes=200, 
                      min_cells=3, 
                      gtf_file=None):
    
    """
    Function to create an Anndata object from 10X cellranger output (.mtx format).
    
    Processing steps:
    1. Detects matrix files in the directory.
    2. Handles single sample or multiple sample aggregation/concatenation.
    3. Optionally subsets genes based on a provided GTF file.
    4. Merges external metadata.
    5. Performs basic filtering (min genes per cell, min cells per gene).
    
    Input:
    cellranger_dir: Directory containing cellranger 'matrix.mtx.gz' output matrix.
    out_dir: Directory to write Anndata object and plots to.
    out_name: Name of dataset, used in plot/file names.
    metadata: Tab-separated file with 'barcode' column matching cell barcodes.
    gtf_file: Optional GTF filepath to subset features (by gene_name).
    """

    # --- Step 1: File Detection ---    
    print(f"Reading files in {cellranger_dir}")

    # Identify all 'matrix.mtx.gz' files in the target directory
    matrix_file_list = []
    with os.scandir(cellranger_dir) as it:
        for entry in it:
            if entry.name.endswith('matrix.mtx.gz') and entry.is_file():
                matrix_file_list.append(entry.name)

    # --- Step 2: Data Loading & Concatenation ---
    if len(matrix_file_list) == 0:
        print(f"Error: Did not find a 'matrix.mtx.gz' file in {cellranger_dir}")
    elif len(matrix_file_list) == 1:
        adata = sc.read_10x_mtx(
            cellranger_dir,
            var_names='gene_symbols',                
            cache=True)
    elif len(matrix_file_list) > 1:
        samples = []
        for entry in matrix_file_list:
            sample_name = entry.replace('_matrix.mtx.gz', '')
            s = sc.read(
                filename=f"{cellranger_dir}/{sample_name}_matrix.mtx.gz",
                cache=False
                ).T
            genes = pd.read_csv(f"{cellranger_dir}/{sample_name}_features.tsv.gz", header=None, sep='\t')
            s.var_names = genes[0]
            s.var['gene_symbols'] = genes[1].values
            s.obs_names = pd.read_csv(f"{cellranger_dir}/{sample_name}_barcodes.tsv.gz", 
                                      header=None)[0]
            s.obs['sample_id'] = sample_name
            samples.append(s)
        
        adata = ad.concat(samples, join='outer', index_unique="_", merge='first')
        adata.var['gene_ids'] = adata.var_names
        adata.var_names = adata.var['gene_symbols']
        adata.var.index.name = None
    
    adata.obs['barcode'] = adata.obs.index
    adata.var_names_make_unique()    
    
    # --- Step 3: GTF-based Gene Filtering (Optional) ---
    if gtf_file is not None:
        print(f"Filtering genes based on GTF file: {gtf_file}")
        allowed_genes = set()
        
        # Regex to extract gene_name value
        gene_name_pattern = re.compile(r'gene_name "([^"]+)"')
        
        try:
            with open(gtf_file, 'r') as f:
                for line in f:
                    # Skip header lines
                    if line.startswith('#'): continue
                    
                    parts = line.strip().split('\t')
                    
                    # Ensure line structure and check for 'gene' feature type
                    if len(parts) >= 9 and parts[2] == 'gene':
                        attributes = parts[8]
                        match = gene_name_pattern.search(attributes)
                        if match:
                            allowed_genes.add(match.group(1))
            
            print(f"Found {len(allowed_genes)} unique genes in GTF file.")
            
            # Subset the AnnData object
            # Identify genes in adata that are in the allowed set
            genes_to_keep = adata.var_names.isin(allowed_genes)
            n_dropped = adata.n_vars - genes_to_keep.sum()
            
            adata = adata[:, genes_to_keep].copy()
            print(f"Subsetted AnnData object. Kept {adata.n_vars} genes (dropped {n_dropped}).")
            
        except Exception as e:
            print(f"Warning: Failed to process GTF file. Skipping filtering. Error: {e}")

    # --- Step 4: Metadata Integration ---  
    if metadata!=None:
        metadata_df = pd.read_csv(metadata, sep='\t')
        # Ensure indices match for mapping
        metadata_df = metadata_df.set_index('barcode')
        # Map specific columns you need
        for col in metadata_df.columns:
            adata.obs[col] = adata.obs.index.map(metadata_df[col])
    
    # --- Step 5: Basic Filtering & Output ---
    sc.pl.highest_expr_genes(adata, n_top=20, save='.pdf')
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    # save final scanpy output file
    adata.write(f"{out_dir}/{out_name}.h5ad", compression='gzip')
    
    return(adata)



def scanpy_create_10X_h5(cellranger_dir, 
                         out_dir, 
                         out_name, 
                         metadata=None, 
                         min_genes=200, 
                         min_cells=3,
                         gtf_file=None):
    
    """
    Function to create an Anndata object from 10X-formatted hdf5 files (.h5).
    Similar logic to 'scanpy_create_10X' but optimized for HDF5 input.
    
    Input:
    cellranger_dir: Directory containing 10X hdf5 files.
    out_dir: Directory to write Anndata object and plots to.
    out_name: Name of dataset.
    metadata: Tab-separated file. 'sample_id' col must match h5 filenames (minus extension).
    gtf_file: Optional filepath to a GTF file for gene subsetting.
    """
    
    # --- Step 1: H5 File Detection ---
    h5_file_list = []
    with os.scandir(cellranger_dir) as it:
        for entry in it:
            if entry.name.endswith('.h5') and entry.is_file():
                h5_file_list.append(entry.name)
    
    # --- Step 2: Load and Concatenate ---
    samples = []
    for index,entry in enumerate(h5_file_list):
        s = sc.read_10x_h5(filename=f"{cellranger_dir}/{entry}")
        s.obs.index = s.obs.index + "_" + str(index)
        s.obs['sample_id'] = entry.replace('.h5', '')
        s.var_names_make_unique()
        samples.append(s)
        
    adata = ad.concat(samples, join='outer', merge='first')
    adata.var.index.name = None
    adata.obs['barcode'] = adata.obs.index
    adata.var_names_make_unique()    
    
    # --- Step 3: GTF-based Gene Filtering (Optional) ---
    if gtf_file is not None:
        print(f"Filtering genes based on GTF file: {gtf_file}")
        allowed_genes = set()
        
        # Regex to extract gene_name value
        gene_name_pattern = re.compile(r'gene_name "([^"]+)"')
        
        try:
            with open(gtf_file, 'r') as f:
                for line in f:
                    # Skip header lines
                    if line.startswith('#'): continue
                    
                    parts = line.strip().split('\t')
                    
                    # Ensure line structure and check for 'gene' feature type
                    if len(parts) >= 9 and parts[2] == 'gene':
                        attributes = parts[8]
                        match = gene_name_pattern.search(attributes)
                        if match:
                            allowed_genes.add(match.group(1))
            
            print(f"Found {len(allowed_genes)} unique genes in GTF file.")
            
            # Subset the AnnData object
            # Identify genes in adata that are in the allowed set
            genes_to_keep = adata.var_names.isin(allowed_genes)
            n_dropped = adata.n_vars - genes_to_keep.sum()
            
            adata = adata[:, genes_to_keep].copy()
            print(f"Subsetted AnnData object. Kept {adata.n_vars} genes (dropped {n_dropped}).")
            
        except Exception as e:
            print(f"Warning: Failed to process GTF file. Skipping filtering. Error: {e}")

    # --- Step 4: Metadata Integration ---   
    if metadata != None:
        metadata_df = pd.read_csv(metadata, sep='\t')
        tmp_df = adata.obs.merge(metadata_df, left_on='sample_id', right_on='sample_id', how='left', sort=False)
        tmp_df.index = tmp_df['barcode']
        check_order = adata.obs.index.equals(tmp_df.index)
        print(f".obs order after including metadata matches: {check_order}")
        adata.obs = tmp_df
        adata.obs.index.name = None
   
    # --- Step 5: Basic Filtering & Output ---
    sc.pl.highest_expr_genes(adata, n_top=20, save='.png')
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    # save final scanpy output file
    adata.write(f"{out_dir}/{out_name}.h5ad", compression='gzip')
    
    return(adata)



def anndata_add_gene_symbols(
        adata, 
        blast_df,
        out_dir, 
        out_name,
        ):
    
    """
    Function to annotate current species data with gene symbols from another species
    using BLAST results (outfmt 6).
    
    Logic:
    1. Calculates match quality (length - mismatch).
    2. Prioritizes exact name matches (e.g. gene A in species X = gene A in species Y).
    3. Sorts by match quality and e-value to find the single best hit.
    4. Merges the 'sseqid' (subject sequence ID) into the AnnData object.
    
    Input:
    adata: Anndata object.
    blast_df: DataFrame containing BLAST output.
    out_dir: Output directory.
    out_name: Dataset name prefix.
    """
    
    # compute the absolute number of matching basepairs in each row
    blast_df['matches'] = blast_df['length'] - blast_df['mismatch']
    
    # create column indicating whether uppercase query gene symbol matches uppercase search gene symbol
    blast_df['name_match_bool'] = np.where(blast_df['qseqid'].str.upper() == blast_df['sseqid'].str.upper(), 1, 0)
    
    # sort dataframe + drop duplicated query gene symbol rows
    blast_df = blast_df.sort_values(['qseqid', 'matches', 'evalue', 'name_match_bool'], 
                                    ascending=[True, False, True, False])
    blast_df = blast_df.drop_duplicates(subset='qseqid', keep='first')

    # sort dataframe + drop duplicated target gene symbol rows
    blast_df = blast_df.sort_values(['sseqid', 'matches', 'evalue', 'name_match_bool'], 
                                    ascending=[True, False, True, False])
    blast_df = blast_df.drop_duplicates(subset='sseqid', keep='first')

    # merge Anndata var and blast results
    adata.var = adata.var.join(blast_df[['qseqid', 'sseqid']].set_index('qseqid'), how='left')
    
    # save final scanpy output file
    adata.write(f"{out_dir}/{out_name}_added_gene_symbols_raw.h5ad", compression='gzip')
    
    """
    adata.var['mouse_id'] = adata.var.index
    mask = adata.var['sseqid'].notna()
    adata = adata[:, mask]
    adata.var = adata.var.set_index('sseqid', drop=False)
    adata.var.index.name = None
    """
    
    return(adata)



def doublet_detection(adata, 
                      out_dir, 
                      out_name, 
                      var=None
                      ):
    
    """
    Function to detect doublets (two cells captured as one) using Scrublet.
    
    Input:
    adata: Anndata object.
    out_dir: Output directory.
    out_name: Dataset name.
    var: (Optional) Batch key (e.g., 'sample_id'). If provided, Scrublet runs 
         independently on each batch to account for batch-specific doublet rates.
    """
    
    adata.obs['barcode'] = adata.obs.index
    adata.obs.index.name = None

    if var==None:
        sc.pp.scrublet(adata)
    else:
        sc.pp.scrublet(adata, batch_key=var)
         
    # save final scanpy output file
    adata.write(f"{out_dir}/{out_name}_doublets_detected.h5ad", compression='gzip')
    
    return(adata)



def is_outlier(adata, metric: str, nmads: int):
    """
    Helper function to identify outliers using Median Absolute Deviation (MAD).
    
    Input:
    adata: AnnData object.
    metric: Column name in adata.obs to test (e.g., 'total_counts').
    nmads: Number of MADs away from the median to consider an outlier.
    
    Returns:
    Boolean Series indicating if a cell is an outlier.
    """

    M = adata.obs[metric]
    lower_threshold = np.median(M) - nmads * median_abs_deviation(M)
    upper_threshold = np.median(M) + nmads * median_abs_deviation(M)
    outlier = (M < lower_threshold) | (M > upper_threshold)
    print(f"{metric}: Lower bound of included range is {lower_threshold}, upper bound is {upper_threshold}.")

    return outlier



def scanpy_analysis(adata,
                    s_genes,
                    g2m_genes,
                    out_dir,
                    out_name,
                    leiden_res=0.4,
                    pct_mt_cutoff=10,
                    neighbors=10,
                    pcs=30,
                    batch_var=None,
                    harmony_var=None,
                    drop_ribo_prot="No",
                    genes_drop=None
                    ):
    
    """
    Main pipeline for Scanpy analysis: QC, Normalization, Integration, Clustering, and Marker Genes.
    
    Input:
    adata: Anndata object.
    s_genes: List of S-phase cell cycle genes.
    g2m_genes: List of G2/M-phase cell cycle genes.
    out_dir, out_name: Output paths.
    leiden_res: Resolution for Leiden clustering (higher = more clusters).
    pct_mt_cutoff: Max percentage of mitochondrial reads allowed.
    neighbors, pcs: Parameters for graph construction and PCA.
    batch_var: Variable to calculate HVGs per-batch.
    harmony_var: Variable for Harmony batch correction/integration.
    drop_ribo_prot: "Yes"/"No" to exclude ribosomal proteins from clustering.
    genes_drop: List of specific genes to exclude from clustering (HVGs).
    """
    
    # --- Step 1: Initialization ---
    # Restore raw counts if they exist in .raw
    if not adata.raw is None:
        adata = adata.raw.to_adata()
    
    # --- Step 2: Quality Control (QC) ---
    # Identify mitochondrial genes (starting with MT-, Mt-, mt-)
    adata.var['mt'] = adata.var_names.str.startswith(('MT-', 'Mt-', 'mt-'))

    # Calculate QC metrics (total counts, n_genes, pct_mt)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=[20], log1p=True, inplace=True)
    
    sc.pl.violin(
        adata, 
        ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
        jitter=0.4, 
        multi_panel=True, 
        save='_qc_violin.png'
        )
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save='_qc_scatter_1.png')
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save='_qc_scatter_2.png')
    
    # --- Step 3: Outlier Removal ---
    # Use MAD-based outlier detection for counts and genes
    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", 5)
        | is_outlier(adata, "log1p_n_genes_by_counts", 5)
        | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
        )
    
    # Hard cutoff for Mitochondrial percentage
    adata.obs["mt_outlier"] = adata.obs["pct_counts_mt"] > pct_mt_cutoff
    
    n_outliers = adata.obs['outlier'].sum()
    n_mt_outliers = adata.obs['mt_outlier'].sum()
    print(f"Removed {n_outliers} cells due to outlier total counts or number of expressed genes")
    print(f"Removed {n_mt_outliers} cells due to outlier mitochondrial gene expression ratio")
    
    # Apply filter
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
    
    # Save raw counts layer
    adata.layers["counts"] = adata.X.copy()

    # --- Step 4: Normalization & Cell Cycle Scoring ---
    # Perform log-normalization, save into a new layer and store as .raw object
    adata_norm = sc.pp.normalize_total(adata, target_sum=1e4, copy=True)
    sc.pp.log1p(adata_norm)
    adata.layers["lognorm"] = adata_norm.X.copy()
    adata.raw = adata_norm

    # Score cell cycle
    s_genes = [x for x in s_genes if x in adata.var_names]
    g2m_genes = [x for x in g2m_genes if x in adata.var_names]
    sc.tl.score_genes_cell_cycle(adata, 
                                 s_genes=s_genes, 
                                 g2m_genes=g2m_genes,
                                 layer='lognorm',
                                 use_raw=False)
    
    # --- Step 5: Feature Selection (HVGs) ---
    # Select Highly Variable Genes using Pearson residuals (Flavor seurat_v3 / pearson)
    if batch_var and len(adata.obs[batch_var].unique()) > 1:
        sc.experimental.pp.highly_variable_genes(
            adata,
            flavor="pearson_residuals",
            n_top_genes=3000,
            batch_key=batch_var
            )
    else:
        sc.experimental.pp.highly_variable_genes(
            adata,
            flavor="pearson_residuals",
            n_top_genes=3000
            )
    
    # --- Step 6: Exclusion of unwanted genes from clustering ---
    # 6a. Exclude Ribosomal proteins (optional)
    if drop_ribo_prot == "Yes":
        ribo_genes = adata.var_names.str.contains('^RP[SL]', case=False, regex=True)
        n_ribo = ribo_genes.sum()
        
        if n_ribo > 0:
            print(f"Identified {n_ribo} ribosomal protein genes. Excluding them from HVG.")
            adata.var.loc[ribo_genes, 'highly_variable'] = False

    # 6b. Exclude Cell Cycle genes to prevent clustering by cycle phase
    all_cc_genes = s_genes + g2m_genes
    cc_genes_present = [g for g in all_cc_genes if g in adata.var_names]
    
    if len(cc_genes_present) > 0:
        print(f"Excluding {len(cc_genes_present)} cell cycle genes from highly variable genes.")
        adata.var.loc[cc_genes_present, 'highly_variable'] = False

    # 6c. Exclude Custom genes provided by user
    if genes_drop is not None:
    	genes_drop_present = [g for g in genes_drop if g in adata.var_names]
    	if len(genes_drop_present) > 0:
        	print(f"Excluding {len(genes_drop_present)} custom genes from highly variable genes.")
        	adata.var.loc[genes_drop_present, 'highly_variable'] = False

    # Subset dataset
    n_hvg = adata.var.highly_variable.sum()
    print(f"Subsetting Anndata object to {n_hvg} highly variable features")
    adata = adata[:, adata.var.highly_variable]
    
    # --- Step 7: Scaling & PCA ---
    # Compute pearson residuals (analytic normalization, replaces standard scale())
    sc.experimental.pp.normalize_pearson_residuals(
        adata,
        inplace=True
        )
    
    # compute PCs
    sc.pp.pca(
        adata, 
        svd_solver='arpack', 
        random_state=0
        )
    sc.pl.pca_variance_ratio(
        adata, 
        log=True, 
        save='.pdf'
        )
    sc.pl.pca_loadings(
        adata, 
        components = '1,2,3',
        save='.pdf'
        )
    
    # --- Step 8: Integration & Neighbor Graph ---
    # Optionally integrate data batches with Harmonypy
    if harmony_var==None:
        # compute neighbourhood graph
        sc.pp.neighbors(
            adata, 
            n_neighbors=neighbors, 
            n_pcs=pcs,
            random_state=0
            )
    else:
        sce.pp.harmony_integrate(adata=adata, key=harmony_var, max_iter_harmony=50)
        sc.pp.neighbors(
            adata, 
            n_neighbors=neighbors, 
            n_pcs=pcs, 
            use_rep='X_pca_harmony',
            random_state=0
            )
    
    # --- Step 9: Clustering & Embedding ---
    # Leiden Clustering
    sc.tl.leiden(
        adata, 
        resolution=leiden_res,
        random_state=0
        )
    
    # PAGA (Partition-based graph abstraction)
    sc.tl.paga(adata)
    sc.pl.paga(
        adata, 
        plot=True, 
        edge_width_scale=2, 
        node_size_scale=2, 
        save='.pdf'
        )
    
    # UMAP embedding (initialized with PAGA)
    sc.tl.umap(
        adata, 
        init_pos='paga', 
        random_state=0
        )
    sc.pl.umap(
        adata, 
        color=['leiden'], 
        legend_loc='right margin', 
        title='', 
        frameon=False, 
        ncols=1, 
        save=f"_leiden_{leiden_res}.png"
        )
    
    # --- Step 10: Marker Gene Identification ---
    # Find markers for each cluster using Wilcoxon rank-sum test on raw data
    sc.tl.rank_genes_groups(
        adata, 
        groupby='leiden', 
        method='wilcoxon',
        use_raw=True,
        n_genes=1000
        )
    sc.pl.rank_genes_groups(
        adata, 
        n_genes=25, 
        sharey=False, 
        save='.pdf'
        )
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    markers = pd.DataFrame({group + '_' + key[:1]: result[key][group] 
                            for group in groups for key in ['names', 'pvals']})
    markers.to_csv(f"{out_dir}/{out_name}_markers.csv")
    
    # Save final object
    adata.write(f"{out_dir}/{out_name}_analysed.h5ad", compression='gzip')
    
    return(adata)



def obs_bar(adata, 
            column_1, 
            column_2, 
            out_dir, 
            out_name,
            x_axis_levels=None, 
            fill_levels=None
            ):
    
    """
    Function to generate a stacked bar plot showing the distribution of labels.
    Typically used to visualize Cell Type (col2) proportions per Sample (col1).
    
    Input:
    adata: Anndata object.
    column_1: Variable for X-axis (e.g., 'sample_id').
    column_2: Variable for fill color/stack (e.g., 'predicted_doublet' or cluster).
    out_dir, out_name: Output settings.
    x_axis_levels: Custom order list for X-axis.
    fill_levels: Custom order list for the stacked categories.
    """

    # count cell type labels by leiden cluster
    plot_df = adata.obs.groupby([column_1, column_2]).agg({column_2: 'size'})
    plot_df = plot_df.rename(columns={column_2:'count'})
    plot_df = plot_df.reset_index()
    
    # generate dataframe with percentage counts for each entry in column_2 in separate column
    plot_df_1 = plot_df.pivot(index=column_1, columns=column_2, values='count')
    plot_df_1 = plot_df_1.fillna(0)
    plot_df_1.to_csv(f"{out_dir}/{out_name}_counts_{column_1}_vs_{column_2}.csv")
    plot_df_1 = plot_df_1.div(plot_df_1.sum(axis=1), axis=0)
    
    # adjust x-axis plotting order
    if x_axis_levels != None:
        plot_df_1 = plot_df_1.reindex(labels=x_axis_levels, axis=0)
        
    # adjust stacked bar plotting order
    if fill_levels != None:
        for fill_level in fill_levels:
            if fill_level not in plot_df_1.columns:
                plot_df_1[fill_level] = 0
        plot_df_1 = plot_df_1[reversed(fill_levels)]
    
    # generate stacked bar plot
    values = plot_df_1.to_dict(orient='list')
    labels = plot_df_1.index
    bottom = np.zeros(len(labels))
    width = 0.8
    #fig_width = len(labels) * 0.7
    
    fig, ax = plt.subplots(figsize=(14,9), frameon=False)
    
    for key, value in values.items():
        p = ax.bar(labels, value, width, label=key, bottom=bottom)
        bottom += value
        #ax.bar_label(p, label_type='center')
        
    ax.set_axisbelow(True)
    ax.set_title("")
    ax.set_xlabel(column_1, fontsize=13)
    ax.set_ylabel("Relative count", fontsize=13)
    ax.tick_params(axis='x', which='major', labelsize=13, labelrotation=30)
    ax.tick_params(axis='y', which='major', labelsize=13)
    ax.tick_params(axis='both', which='minor', labelsize=13)
    for tick in ax.xaxis.get_majorticklabels():
        tick.set_horizontalalignment("right")
    
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title=column_2, 
              title_fontsize=12, fontsize=13, frameon=False)
    
    ax.figure.savefig(f"{out_dir}/{out_name}_{column_1}_vs_{column_2}_stacked_bar.pdf",
                      dpi=300)
    
    plt.close()
    return(plot_df_1)





### Analysis
# create output directory
if not os.path.exists(out_dir):
    os.makedirs(out_dir, exist_ok=True)

# read in genes for cell cycle scoring
cellcycle_genes = pd.read_csv(cellcycle_genes_csv)
s_genes = cellcycle_genes['G1_S'].dropna().to_list()
g2m_genes = cellcycle_genes['G2_M'].dropna().to_list()

# read in list of genes to drop from HVG
file_size = os.path.getsize(genes_drop_csv)
if file_size == 0:
    genes_drop = None
else:
    print("Reading in custom genes to drop from HVG.")
    genes_drop = pd.read_csv(genes_drop_csv, header=None).iloc[:, 0].tolist()

# create Anndata object
if soupx_dir and os.path.exists(soupx_dir):
    cellranger_dir = soupx_dir
else:
    cellranger_dir = in_dir

# check how many 'matrix.mtx.gz' files there are in the input directory
matrix_file_list = []
with os.scandir(cellranger_dir) as it:
    for entry in it:
        if entry.name.endswith('matrix.mtx.gz') and entry.is_file():
            matrix_file_list.append(entry.name)
    
if len(matrix_file_list) > 0:
    if not os.path.isfile(f"{out_dir}/{out_name}.h5ad"):
        adata = scanpy_create_10X(cellranger_dir=cellranger_dir, 
                                  out_dir=out_dir, 
                                  out_name=out_name,
                                  metadata=metadata,
                                  min_genes=min_genes, 
                                  min_cells=min_cells,
                                  gtf_file=None)
    else:
        adata = sc.read_h5ad(f"{out_dir}/{out_name}.h5ad")
else:
    # check how many h5 files there are in the directory
    h5_file_list = []
    with os.scandir(cellranger_dir) as it:
        for entry in it:
            if entry.name.endswith(".h5") and entry.is_file():
                h5_file_list.append(entry.name)

    if len(h5_file_list) > 0:
        if not os.path.isfile(f"{out_dir}/{out_name}.h5ad"):
            adata = scanpy_create_10X_h5(cellranger_dir=cellranger_dir, 
                                         out_dir=out_dir, 
                                         out_name=out_name,
                                         metadata=metadata,
                                         min_genes=min_genes, 
                                         min_cells=min_cells,
                                         gtf_file=None)
        else:
            adata = sc.read_h5ad(f"{out_dir}/{out_name}.h5ad")
    else:
        print("Please provide an input directory that contains either \
              a Cell Ranger 'matrix.mtx.gz' file or '.h5' files with mapped data.")



# add human gene symbols
if blast_file and not os.path.isfile(f"{out_dir}/{out_name}_added_gene_symbols_raw.h5ad"):
    file_size = os.path.getsize(blast_file)

    if file_size == 0:
        print("Skipping gene symbol addition: Blast file is empty.")
    else:
        blast_df = pd.read_csv(blast_file, sep='\t')
        adata = anndata_add_gene_symbols(adata=adata, 
                                         blast_df=blast_df,
                                         out_dir=out_dir, 
                                         out_name=out_name)
elif os.path.isfile(f"{out_dir}/{out_name}_added_gene_symbols_raw.h5ad"):
    adata = sc.read_h5ad(f"{out_dir}/{out_name}_added_gene_symbols_raw.h5ad")
else:
    # Handle case where no blast file is provided and no pre-processed file exists
    print("Skipping gene symbol addition: Blast file not provided.")



# plot cell numbers per sample + exclude samples with number of cells below threshold
threshold = 50 

counts = adata.obs[harmony_var].value_counts().sort_values(ascending=True)
plt.figure(figsize=(12, 8))
counts.plot(kind='barh', color='skyblue', edgecolor='black')

# Add a vertical line for the threshold
plt.axvline(x=threshold, color='red', linestyle='--', label=f'Threshold ({threshold})')

plt.title('Number of Cells per Sample')
plt.xlabel('Cell Count')  # Now on the X-axis
plt.ylabel('Sample ID')   # Now on the Y-axis
plt.legend()
plt.tight_layout()
plt.savefig(f"{out_dir}/{out_name}_cell_counts_horizontal.pdf")

# filter the AnnData object
samples_to_remove = counts[counts < threshold].index.tolist()
samples_to_keep = counts[counts >= threshold].index.tolist()
adata = adata[adata.obs[harmony_var].isin(samples_to_keep)]

print("Removed samples: ")
print(samples_to_remove)

# clean up categorical metadata
if isinstance(adata.obs[harmony_var].dtype, pd.CategoricalDtype):
    adata.obs[harmony_var] = adata.obs[harmony_var].cat.remove_unused_categories()



# detect doublets
if not os.path.isfile(f"{out_dir}/{out_name}_doublets_detected.h5ad"):
    adata = doublet_detection(adata=adata, 
                              out_dir=out_dir, 
                              out_name=out_name, 
                              var=harmony_var)
    
    # analyse doublet frequency across samples
    plot_bar = obs_bar(adata=adata, 
                       column_1=harmony_var, 
                       column_2="predicted_doublet", 
                       out_dir=out_dir, 
                       out_name=out_name
                       )
else:
    adata = sc.read_h5ad(f"{out_dir}/{out_name}_doublets_detected.h5ad")



# scanpy analysis
if not os.path.isfile(f"{out_dir}/{out_name}_scRNAseq_analysed_no_doublets.h5ad"):
    
    # remove doublets
    print(f"Removing {sum(adata.obs['predicted_doublet'] == 1)} doublets.")
    adata = adata[adata.obs['predicted_doublet'] == 0].copy()
    
    if not harmony_var in adata.obs.columns:
        harmony_var=None
    
    adata = scanpy_analysis(adata=adata,
                            s_genes=s_genes,
                            g2m_genes=g2m_genes,
                            out_dir=out_dir, 
                            out_name=out_name,
                            batch_var=harmony_var,
                            harmony_var=harmony_var, 
                            pct_mt_cutoff=max_perc_mt,
                            neighbors=10, 
                            pcs=n_pcs, 
                            leiden_res=leiden_res,
                            drop_ribo_prot=drop_ribo_prot,
                            genes_drop=genes_drop
                            )
    
    # plot
    sc.pl.umap(adata, color=['sample_id'], 
               legend_loc='right margin', title='', frameon=False, ncols=1, 
               save='_sample_id.png')
    
    sc.pl.umap(adata, color=['leiden'], 
               legend_loc='right margin', title='', frameon=False, ncols=1, 
               save=f'_leiden_{leiden_res}.png')

    # save final scanpy output file
    adata.obs['dataset'] = out_name
    #adata.obs['sample_id'] = adata.obs['sample_id'].cat.rename_categories(lambda x: x + "_" + out_name)

    print(adata)
    print(adata.obs['sample_id'].value_counts())
    adata.obs.to_csv(f"{out_dir}/{out_name}_obs.csv")
    adata.var.to_csv(f"{out_dir}/{out_name}_var.csv")
    adata.write(f"{out_dir}/{out_name}_scRNAseq_analysed_no_doublets.h5ad", compression='gzip')

