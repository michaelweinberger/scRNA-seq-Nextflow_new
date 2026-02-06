#!/usr/bin/env python


# -*- coding: utf-8 -*-
"""
Integrate SCINA automated cell annotation results into Scanpy object

"""


### Modules
import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from adjustText import adjust_text
import anndata as ad
from scipy.sparse import csr_matrix
import argparse
import re

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.autosave = True	          # do not show plots, save them in figdir
sc.settings.set_figure_params(dpi=300, facecolor='white')


# --- ARGUMENT PARSING ---
## unpack arguments imported from bash parent script
parser = argparse.ArgumentParser(description='Integrate SCINA automated cell annotation results into Scanpy object')
parser.add_argument('-i','--scanpy_obj', help='Filepath to an Anndata object processed with Scanpy', required=True)
parser.add_argument('-j','--scanpy_raw_obj', help='Filepath to an Anndata object with full gene count matrix', required=True)
parser.add_argument('-s','--scina_results_csv', help='Filepath to a csv file containing SCINA results with columns \
                    cell_barcode, SCINA_predicted_celltype, and SCINA_confidence', required=True)
parser.add_argument('-m','--cell_type_markers_csv', help='Filepath to a csv file containing cell type markers', required=True)
parser.add_argument('-t','--cell_type_metadata_csv', help='Filepath to a csv file containing cell type metadata for ordering and colors, with columns: cell_type, ID, and color', required=True)           
parser.add_argument('-r','--leiden_res', help='Resolution of Leiden clustering', required=False)
parser.add_argument('-o','--output_dir', help='Directory for output files', required=True)
parser.add_argument('-n','--output_name', help='Prefix for output file names', required=True)

args = parser.parse_args()
adata_file = args.scanpy_obj
adata_raw_file = args.scanpy_raw_obj
scina_results_csv = args.scina_results_csv
cell_type_markers_csv = args.cell_type_markers_csv
cell_type_metadata_csv = args.cell_type_metadata_csv
leiden_res = float(args.leiden_res) if args.leiden_res else None 
out_dir = args.output_dir
out_name = args.output_name


### Functions

def plot_umap_with_annotations(adata, group_key, color_key, out_dir, out_name):
    """
    Generates an optimized UMAP plot where labels are automatically 
    repositioned to avoid overlapping cells and each other.

    Parameters
    ----------
    adata : anndata.AnnData
        The annotated data matrix of shape (n_obs, n_vars).
    group_key : str
        The column name in `adata.obs` containing the categories to label (e.g., 'cell_type_ID').
    color_key : str
        The key in `adata.uns` containing the color mapping for the categories.
    out_dir : str
        Path to the output directory where the plot will be saved.
    out_name : str
        Prefix for the output filename.

    Returns
    -------
    None
        Saves a PDF file to `out_dir` named `{out_name}_{group_key}_annotated.pdf`.

    Notes
    -----
    - Calculates the centroid and standard deviation of each cluster in UMAP space.
    - Uses `adjustText` to iteratively repel labels away from dense point clouds and other labels.
    """

    print(f"\n--- Generating Optimized UMAP plot for {group_key} ---")
    
    # 1. Validation
    if 'X_umap' not in adata.obsm:
        return print("Warning: 'X_umap' not found.")
    if group_key not in adata.obs.columns or color_key not in adata.uns:
        return print(f"Error: Missing {group_key} or {color_key}.")
    
    categories = adata.obs[group_key].cat.categories
    colors = adata.uns[color_key]
    color_map = dict(zip(categories, colors))
    
    # 2. Setup Plot
    fig, ax = plt.subplots(figsize=(9, 9))
    
    # Capture the scatter plot object to tell adjust_text what to avoid
    sc.pl.umap(
        adata,
        color=group_key,
        title=f"Cell Types Annotated ({group_key})",
        frameon=False,
        legend_loc='none',
        show=False,
        ax=ax
    )
    
    # 3. Calculate Coordinates
    umap_coords = pd.DataFrame(adata.obsm['X_umap'], index=adata.obs_names, columns=['UMAP1', 'UMAP2'])
    plot_df = pd.concat([umap_coords, adata.obs[group_key]], axis=1)
    group_stats = plot_df.groupby(group_key).agg(['mean', 'std'])
    global_center = umap_coords.mean()
    
    texts = []
    
    # 4. Create initial label positions
    for category in categories:
        if category in group_stats.index:
            cx, cy = group_stats.loc[category, ('UMAP1', 'mean')], group_stats.loc[category, ('UMAP2', 'mean')]
            sx, sy = group_stats.loc[category, ('UMAP1', 'std')], group_stats.loc[category, ('UMAP2', 'std')]
            
            # Use your vector logic to get the "starting" position outside the cluster
            dx, dy = cx - global_center[0], cy - global_center[1]
            norm = np.sqrt(dx**2 + dy**2)
            shift = 2.5 * max(sx, sy) if pd.notna(sx) else 0.5
            
            x_init = cx + (dx/norm * shift) if norm > 0 else cx + 0.5
            y_init = cy + (dy/norm * shift) if norm > 0 else cy + 0.5
            
            # Create the text object
            t = ax.text(
                x_init, y_init, category, 
                fontsize=9, weight='bold',
                color=color_map[category]
            )
            texts.append(t)

    # 5. Optimize label positions
    # adjust_text will move 'texts' to avoid overlapping with:
    # 1. Other texts in the 'texts' list
    # 2. The scatter points (ax.collections)
    scatter_obj = ax.collections[0]
    adjust_text(
        texts,
        add_objects=[scatter_obj],
        expand_points=(1.5, 1.5), # Buffer around points
        expand_text=(1.2, 1.2),   # Buffer around labels
        ax=ax
    )
    
    # Save
    plot_file = os.path.join(out_dir, f"{out_name}_{group_key}_annotated.pdf")
    fig.savefig(plot_file, bbox_inches='tight')
    plt.close(fig)
    print(f"Annotated UMAP saved: {plot_file}")



def plot_frequency_bar_chart(freq_df, out_dir, out_name, cell_type_colors=None):
    """
    Generates a multi-panel bar chart showing the frequency of SCINA predicted cell types 
    per cluster.

    Parameters
    ----------
    freq_df : pandas.DataFrame
        DataFrame containing columns 'cluster', 'SCINA_predicted_celltype', and 'Frequency'.
    out_dir : str
        Directory to save the plot.
    out_name : str
        Prefix for output filename.
    cell_type_colors : dict, optional
        Dictionary mapping cell type names to hex color codes. If None, uses 'viridis'.

    Returns
    -------
    None
        Saves a PDF file to `out_dir`.
    """
    
    # Replace underscores with spaces in SCINA_predicted_celltype labels
    freq_df['SCINA_predicted_celltype'] = freq_df['SCINA_predicted_celltype'].astype(str).str.replace('_', ' ')
    
    # Ensure cluster and cell type are categorical for consistent plotting
    freq_df['cluster'] = freq_df['cluster'].astype('category')
    freq_df['SCINA_predicted_celltype'] = freq_df['SCINA_predicted_celltype'].astype('category')
    
    # Create subplots for each cluster
    n_clusters = freq_df['cluster'].nunique()
    # Determine grid layout: 3 columns, auto rows
    n_cols = 3
    n_rows = int(np.ceil(n_clusters / n_cols))
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows), constrained_layout=True)
    axes = axes.flatten()
    
    # Prepare color mapping
    if cell_type_colors:
        # Map cell type to a specific color
        cell_type_to_color = lambda ctype: cell_type_colors.get(ctype, 'grey')
        # Use a consistent color map for all clusters
        cmap = None
    else:
        # Fallback to a default colormap if no colors are provided
        cmap = plt.get_cmap('viridis')
        
    for i, (cluster_id, group) in enumerate(freq_df.groupby('cluster')):
        ax = axes[i]
        
        # Sort by frequency descending
        group = group.sort_values(by='Frequency', ascending=False)
        
        # Set colors based on the map if provided, otherwise use default cmap
        if cell_type_colors:
            colors = [cell_type_to_color(ctype) for ctype in group['SCINA_predicted_celltype']]
        else:
            # Use the default viridis if no colors are provided
            colors = cmap(np.arange(len(group)))
            
        ax.bar(group['SCINA_predicted_celltype'], group['Frequency'], 
               color=colors, 
               edgecolor='black')        
        ax.set_title(f"Cluster {cluster_id}")
        ax.set_xlabel("SCINA Predicted Cell Type")
        ax.set_ylabel("Relative Frequency")
        ax.tick_params(axis='x', rotation=45)
        ax.grid(axis='y', linestyle='--')
        for tick in ax.xaxis.get_majorticklabels():
            tick.set_horizontalalignment("right")
        
    # Hide unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])
        
    fig.suptitle(f"Top SCINA Cell Type Frequencies per Cluster: {out_name}", fontsize=16, fontweight='bold')
    
    freq_plot_file = os.path.join(out_dir, f"{out_name}_scina_cluster_frequency_plot.pdf")
    plt.savefig(freq_plot_file, bbox_inches='tight')
    plt.close()
    print(f"Cluster frequency plot saved to: {freq_plot_file}")



def integrate_scina_results(adata, 
                            adata_raw, 
                            scina_results_df, 
                            marker_cols_df, 
                            cell_type_meta_df, 
                            out_dir, 
                            out_name,
                            n_markers=200, 
                            heatmap_n=5, 
                            violin_n=3, 
                            umap_marker_n=6, 
                            cluster_key='leiden'):

    """
    Integrates SCINA results into a Scanpy AnnData object, determines consensus
    cell type labels per cluster, performs differential expression analysis, 
    and generates various plots and exports top markers.

    Parameters
    ----------
    adata : anndata.AnnData
        The main processed Scanpy object (clustered, e.g., via Leiden).
    adata_raw : anndata.AnnData
        The raw count matrix object (used for expression plotting).
    scina_results_df : pandas.DataFrame
        Dataframe containing 'cell_barcode', 'SCINA_predicted_celltype', and 'SCINA_confidence'.
    marker_cols_df : pandas.DataFrame
        Dataframe where columns are cell types and rows are marker genes.
    cell_type_meta_df : pandas.DataFrame
        Metadata defining 'cell_type', 'ID', and 'color' for visualization.
    out_dir : str
        Directory for output files.
    out_name : str
        Prefix for output files.
    n_markers : int, optional (default=200)
        Number of top markers to export to Excel per cell type.
    heatmap_n : int, optional (default=5)
        Number of markers per cell type to include in the dotplot heatmap.
    violin_n : int, optional (default=3)
        Number of markers per cell type to include in the stacked violin plot.
    umap_marker_n : int, optional (default=6)
        Number of markers to extract for subsequent UMAP feature plotting.
    cluster_key : str, optional (default='leiden')
        The column in `adata.obs` containing cluster assignments.

    Returns
    -------
    tuple
        (adata_annotated, adata_full, umap_markers_list, cell_type_colors)
        - adata_annotated: The input adata with new 'cell_type' and 'cell_type_ID' columns.
        - adata_full: A subset of `adata_raw` normalized for plotting.
        - umap_markers_list: Dictionary of top marker genes for the consensus cell types.
        - cell_type_colors: Dictionary mapping cell types to colors.
    """
    
    print("\n--- Integrating SCINA labels and calculating cluster frequencies ---\n")
    
    # --- 1. Integrate Results back into AnnData Metadata (adata.obs) ---
    # Ensure scina_results_df is indexed by cell barcode
    if 'cell_barcode' in scina_results_df.columns:
        scina_results_df = scina_results_df.set_index('cell_barcode')
        
    # Check and ensure the order of labels matches the order of cells
    if not all(adata.obs_names.isin(scina_results_df.index)):
        raise ValueError("Fatal error: Predicted cell labels indices do not match AnnData object cell barcodes.")
    
    # Align and add the predictions to adata.obs
    predicted_labels = scina_results_df.loc[adata.obs_names, 'SCINA_predicted_celltype']
    adata.obs['SCINA_predicted_celltype'] = predicted_labels.astype('category')
    adata.obs['SCINA_confidence'] = scina_results_df.loc[adata.obs_names, 'SCINA_confidence']
    
    print("SCINA labels successfully added to AnnData object metadata ('SCINA_predicted_celltype').")
    
    # --- 2. Check for Clusters and Prepare Data ---
    if cluster_key not in adata.obs.columns:
        print(f"Warning: Metadata column '{cluster_key}' not found. Skipping frequency analysis and consensus labeling.")
        return (adata, {}, {})
    
    meta_data_df = adata.obs[[cluster_key, 'SCINA_predicted_celltype']].copy()
    meta_data_df.rename(columns={cluster_key: 'cluster'}, inplace=True) # Temporarily rename cluster column
    
    # Calculate absolute counts per cluster/cell type
    counts_df = meta_data_df.groupby(['cluster', 'SCINA_predicted_celltype']).size().reset_index(name='Count')
    
    # Calculate relative frequency (proportion) within each cluster
    counts_df['Frequency'] = counts_df['Count'] / counts_df.groupby('cluster')['Count'].transform('sum')
    freq_df = counts_df
    
    # --- 3. DETERMINE AND INTEGRATE CONSENSUS LABEL ---
    # Find the single most frequent cell type for each cluster (the consensus annotation)
    consensus_labels = freq_df.loc[freq_df.groupby('cluster')['Frequency'].idxmax()]
    
    consensus_labels = consensus_labels[['cluster', 'SCINA_predicted_celltype']].rename(
        columns={'SCINA_predicted_celltype': 'SCINA_cluster_label'}
    )
    
    # Map the consensus label back to all cells in adata.obs
    cluster_to_label = consensus_labels.set_index('cluster')['SCINA_cluster_label'].to_dict()
    
    # Add the new consensus column to adata.obs
    # Initial assignment uses the raw SCINA label
    adata.obs['cell_type'] = adata.obs[cluster_key].astype(str).map(cluster_to_label).astype('category')
    
    print("Consensus SCINA label ('cell_type') added to AnnData object metadata.")
    
    # --- Re-leveling 'cell_type' based on cell_type_meta_df ---
    print("Re-leveling 'cell_type' categories based on 'cell_type' column order in metadata file...")
    
    # 3.1 Standardize labels in adata.obs
    adata.obs['cell_type'] = adata.obs['cell_type'].astype(str).str.replace('_', ' ').astype('category')
    print("Underscores replaced with spaces in 'cell_type' labels for all Scanpy plots.")
    
    # 3.2 Get the desired order from the metadata file, standardizing to match the consensus labels
    # Use .unique() to handle potential duplicates in the metadata file, keeping order
    new_cell_type_order = cell_type_meta_df['cell_type'].astype(str).str.replace('_', ' ').unique().tolist()
    
    # 3.3 Filter the order to only include cell types present in the current adata object
    present_categories = set(adata.obs['cell_type'].unique())
    
    # Prepare the final order by filtering
    final_order = [
        ctype for ctype in new_cell_type_order if ctype in present_categories
    ]
    
    # 3.4 If consensus cell types have been correctly assigned, re-level the column
    if final_order:
        adata.obs['cell_type'] = pd.Categorical(values=adata.obs.cell_type, categories=final_order, ordered=True)
        print(f"'cell_type' categories re-leveled to: {final_order[:5]}... ({len(final_order)} categories total)")
    else:
        print("Warning: Could not re-level 'cell_type'. No matching cell types found between consensus labels and 'cell_type' column in metadata file.")
        
    # --- Extract and Set Color Map from cell_type_meta_df ---
    cell_type_colors = {}
        
    # Map cell type names (standardized) to colors and IDs
    meta_df_standardized = cell_type_meta_df.copy()
    meta_df_standardized['cell_type'] = meta_df_standardized['cell_type'].astype(str).str.replace('_', ' ')
    meta_df_standardized = meta_df_standardized.set_index('cell_type')
    
    cell_type_colors_all = dict(zip(meta_df_standardized.index, meta_df_standardized['color']))
    
    # Create color dictionary based on the final order
    colors_for_plot = [
        meta_df_standardized.loc[ctype, 'color'] for ctype in final_order
    ]
    
    if colors_for_plot and all(pd.notna(colors_for_plot)):
        # Set the color palette for full cell type names
        adata.uns['cell_type_colors'] = colors_for_plot
        cell_type_colors = dict(zip(final_order, colors_for_plot))
        print("Set 'cell_type' color map in adata.uns based on provided metadata colors.")
    else:
        # Fallback if no colors are provided or are missing
        print("Warning: Could not set 'cell_type_colors' from metadata. Default colors will be used.")
        
    # Map shortened IDs to a new column 'cell_type_ID'
    try:
        id_map = meta_df_standardized['ID'].to_dict()
        # Map only the present categories
        id_map_present = {k: v for k, v in id_map.items() if k in present_categories}
        
        # Use the ID map for the new column
        adata.obs['cell_type_ID'] = adata.obs['cell_type'].map(id_map_present)
        
        # Filter the order of IDs based on the final order of cell types
        id_order = [id_map_present[ctype] for ctype in final_order]
        adata.obs['cell_type_ID'] = pd.Categorical(values=adata.obs.cell_type_ID, categories=id_order, ordered=True)
        
        # Also set the color map for the ID column for plotting
        adata.uns['cell_type_ID_colors'] = colors_for_plot
        
        print("Shortened cell type labels ('cell_type_ID') and corresponding colors added to AnnData object metadata.")
        
        # Use new custom plotting function for 'cell_type_ID'
        plot_umap_with_annotations(
            adata, 
            group_key='cell_type_ID', 
            color_key='cell_type_ID_colors', 
            out_dir=out_dir, 
            out_name=out_name
        )
        
    except Exception as e:
        print(f"Warning: Could not map 'cell_type' to 'ID' from metadata. Skipping 'cell_type_ID' creation. Error: {e}")
        # Plot UMAP Feature plot with cell type overlay (using full name) as fallback
        sc.pl.umap(adata, color=['cell_type'], 
                   legend_loc='right margin', title='', frameon=False, ncols=1, 
                   save=f"_{out_name}_cell_type.png")
        pass # Continue without the ID column if mapping fails
    
    # Plot UMAP Feature plot with cell type overlay (using full name, regular legend)
    sc.pl.umap(adata, color=['cell_type'], 
               legend_loc='right margin', title='', frameon=False, ncols=1, 
               save=f"_{out_name}_cell_type.png")
    
    # --- 4. Filter for Top Cell Labels per Cluster (for plotting) ---
    # Filter to show only the top most frequent cell types in each cluster
    top_freq_df = freq_df.sort_values(by=['cluster', 'Frequency'], ascending=[True, False])
    top_freq_df.to_csv(f"{out_dir}/{out_name}_top_freq_scina_celltypes_per_cluster.csv")
        
    # --- 5. Plot Cluster Frequency Bar Chart ---
    top_freq_df = top_freq_df.groupby('cluster').head(4).reset_index(drop=True)
    # Pass the extracted/set color dictionary to the plotting function
    plot_frequency_bar_chart(top_freq_df, out_dir, out_name, cell_type_colors=cell_type_colors_all)
    
    # --- 6. FIND MARKERS, EXPORT, AND GENERATE PLOTS ---
    print("Finding markers for SCINA consensus labels...")
    
    # Ensure 'cell_type' is used as the grouping variable
    group_by_key = 'cell_type'
    
    # Find markers
    sc.tl.rank_genes_groups(
        adata, 
        groupby=group_by_key, 
        method='wilcoxon',
        use_raw=False,
        layer='lognorm',
        key_added='scina_markers'
    )
    
    # Convert results to a pandas DataFrame (long format)
    markers_raw = sc.get.rank_genes_groups_df(adata, key='scina_markers', group=None)
    
    # Filter for positive log fold changes (equivalent to only.pos=TRUE)
    markers = markers_raw[markers_raw['logfoldchanges'] > 0]
    markers = markers.rename(columns={'logfoldchanges': 'avg_log2FC'})
    markers = markers.sort_values(by=['group', 'avg_log2FC'], ascending=[True, False])
    
    if markers.empty:
        print("Warning: No significant positive markers found for consensus cell types. Skipping remaining steps.")
        return (adata, {}, cell_type_colors)
    
    # --- 6.1 EXCEL EXPORT (Top n_markers) ---
    print(f"Saving top {n_markers} markers to Excel file...")
    
    top_excel_markers = markers.groupby('group').head(n_markers)
    
    # Prepare for multi-sheet Excel export (one sheet per cell type/group)
    top_markers_df_list = {
        name: group for name, group in top_excel_markers.groupby('group')
    }
    
    marker_excel_file = os.path.join(out_dir, f"{out_name}_scina_cell_type_top_markers.xlsx")
    
    with pd.ExcelWriter(marker_excel_file) as writer:
        for sheet_name, df in top_markers_df_list.items():
            df.to_excel(writer, sheet_name=str(sheet_name), index=False)
            
    print(f"Top markers saved to: {marker_excel_file}")
    
    # subset raw object to barcodes contained in clustered object
    adata_full = adata_raw[adata.obs_names, :].copy()
    adata_full.obs = adata.obs.copy()
    adata_full.obsm = adata.obsm.copy()
    adata_full.uns = adata.uns.copy()

    # normalize + scale data
    sc.pp.normalize_total(adata_full, target_sum=1e4)
    sc.pp.log1p(adata_full)

    # --- 6.2 HEATMAP DOTPLOT (Using Top heatmap_n Markers from marker_cols_df) ---
    print(f"\n--- Preparing Top {heatmap_n} markers from marker_cols_df for Heatmap Dotplot ---")
    
    # Standardize marker_cols_df column names for matching with adata.obs['cell_type'] groups
    marker_cols_df.columns = marker_cols_df.columns.str.replace('_', ' ')
    
    # Filter markers from marker_cols_df based on consensus cell types and presence in adata.var
    heatmap_genes = []
    # Get the cell types present in the annotated object (already standardized)
    present_cell_types = adata_full.obs['cell_type'].cat.categories
    all_genes_in_adata = set(adata_full.var_names)
    
    for cell_type in present_cell_types:
        if cell_type in marker_cols_df.columns:
            # Get the top N markers for this cell type from the original marker list
            valid_markers = [
                m for m in marker_cols_df[cell_type] if m in all_genes_in_adata
            ]
            # Take the top 'heatmap_n' valid markers
            heatmap_genes.extend(valid_markers[:heatmap_n])
    
    # Remove duplicates while maintaining order
    heatmap_genes = list(dict.fromkeys(heatmap_genes))
    
    if not heatmap_genes:
         print("Warning: No significant positive markers found for SCINA consensus clusters. Skipping heatmap plot.")
    else:
        heatmap_plot_file = os.path.join(out_dir, f"{out_name}_scina_cell_type_heatmap_top{heatmap_n}_markers.pdf")
        print(f"Saving marker heatmap (Top {heatmap_n} markers) to: {heatmap_plot_file}")
        
        # Scanpy equivalent of DoHeatmap (using DotPlot as proxy)
        dp = sc.pl.dotplot(
            adata_full, 
            var_names=heatmap_genes, 
            groupby=group_by_key, 
            dendrogram=False, 
            cmap='viridis', 
            standard_scale='var', 
            figsize=(25, 9),
            show=False,
            return_fig=True
        )
        dp.savefig(heatmap_plot_file)
        print("Marker dotplot saved.")
        
    # --- 6.3 STACKED VIOLIN PLOT (Using Top violin_n Markers from marker_cols_df) ---
    print(f"\n--- Preparing Top {violin_n} markers from marker_cols_df for Stacked Violin Plot ---")
    
    # Filter markers from marker_cols_df based on consensus cell types and presence in adata.var
    plot_genes = []
    
    for cell_type in present_cell_types:
        if cell_type in marker_cols_df.columns:
            # Get the top N markers for this cell type from the original marker list
            valid_markers = [
                m for m in marker_cols_df[cell_type] if m in all_genes_in_adata
            ]
            # Take the top 'violin_n' valid markers
            plot_genes.extend(valid_markers[:violin_n])
    
    # Remove duplicates while maintaining order
    plot_genes = list(dict.fromkeys(plot_genes))
    
    if not plot_genes:
        print("Warning: No significant positive markers found. Skipping stacked violin plot.")
    else:
        vln_plot_file = os.path.join(out_dir, f"{out_name}_scina_cell_type_stacked_violin_top{violin_n}_markers.pdf")
        print(f"Saving stacked violin plot for top {violin_n} markers to: {vln_plot_file}")
        
        # Scanpy equivalent of Stacked VlnPlot
        vp = sc.pl.stacked_violin(
            adata_full, 
            var_names=plot_genes, 
            groupby=group_by_key, 
            dendrogram=False, 
            figsize=(15, 10),
            row_palette=adata_full.uns['cell_type_colors'] if cell_type_colors else None,
            show=False,
            return_fig=True
        )
        vp.savefig(vln_plot_file)
        print("Stacked violin plot saved.")
        
    # --- 6.4 Consensus Markers ---
    # Filter for the top UMAP_MARKERS (user-defined) per cluster
    top_umap_markers = markers.groupby('group').head(umap_marker_n)
    
    # Convert the filtered DataFrame into the requested list structure: {cell_type: [marker_genes]}
    umap_markers_list = top_umap_markers.groupby('group')['names'].apply(list).to_dict()
    
    # --- 7. Return the AnnData Object and the UMAP Markers List ---
    return (adata, adata_full, umap_markers_list, cell_type_colors)



def plot_umap_feature_plots(adata, 
                            markers_dict, 
                            out_dir, 
                            out_name):

    """
    Generates UMAP plots colored by the expression of cell-type-specific marker genes
    (equivalent to Seurat's FeaturePlot).

    Parameters
    ----------
    adata : anndata.AnnData
        The AnnData object containing expression data and 'X_umap' coordinates.
    markers_dict : dict
        A dictionary where keys are cell type names and values are lists of marker gene strings.
    out_dir : str
        Directory to save the plots.
    out_name : str
        Prefix for the output filenames.

    Returns
    -------
    None
        Saves one or more PNG files (depending on number of genes) to `out_dir`.
    """

    print("\n--- Generating UMAP FeaturePlots for cell-type-specific marker genes ---")
    
    # 1. Check for UMAP reduction
    if 'X_umap' not in adata.obsm:
        print("Warning: 'X_umap' reduction not found in the AnnData object (.obsm). Skipping UMAP FeaturePlot generation. Please run 'sc.tl.umap()' on your AnnData object before calling this function.")
        return
    
    # 2. Get all genes present in the AnnData object
    all_genes_in_adata = set(adata.var_names)
    # Print only the first few genes for brevity/check
    print(f"First 10 genes in AnnData object: {list(all_genes_in_adata)[:10]}")
    
    for cell_type, markers_for_type in markers_dict.items():
        
        # Filter markers: keep only those present in the AnnData object
        valid_markers = [marker for marker in markers_for_type if marker in all_genes_in_adata]
        
        if valid_markers:
            
            # Sanitize cell_type name for file saving
            sanitized_cell_type = re.sub(r'[^A-Za-z0-9_]', '_', cell_type)
            plot_file = f"_{out_name}_Markers_{sanitized_cell_type}.png"
            
            print(f"Saving UMAP FeaturePlots for {cell_type} to: {plot_file}")
            
            # Generate the plot using scanpy
            sc.pl.umap(
                adata,
                color=valid_markers,
                frameon=False,
                save=plot_file,
                show=False,
                ncols=2
            )
        else:
            print(f"Warning: No valid marker genes found for cell type: {cell_type}. Skipping UMAP FeaturePlot.")
                
    print("Finished generating UMAP FeaturePlots for cell type markers.")



### Analysis
sc.settings.figdir = out_dir

if not os.path.exists(out_dir):
    os.makedirs(out_dir, exist_ok=True)

# 1. Read input files
adata = sc.read_h5ad(adata_file)
adata_raw = sc.read_h5ad(adata_raw_file)
scina_results_df = pd.read_csv(scina_results_csv)
marker_cols_df = pd.read_csv(cell_type_markers_csv)
cell_type_meta_df = pd.read_csv(cell_type_metadata_csv)

markers_dict = {col: marker_cols_df[col].dropna().tolist() for col in marker_cols_df.columns}


# 2. Integrate SCINA results into Scanpy object
adata_annotated, adata_full, umap_markers_list, cell_type_colors = integrate_scina_results(
    adata=adata,
    adata_raw=adata_raw,
    scina_results_df=scina_results_df,
    marker_cols_df=marker_cols_df,
    cell_type_meta_df=cell_type_meta_df,
    out_dir=out_dir,
    out_name=out_name,
    n_markers=200,          # For Excel
    heatmap_n=8,            # For Heatmap/DotPlot
    violin_n=2,             # For Stacked Violin Plot
    umap_marker_n=6         # For UMAP FeaturePlot
)


# 3. Plot UMAP feature plots for pre-defined marker genes
plot_umap_feature_plots(
    adata=adata_full, 
    markers_dict=markers_dict, 
    out_dir=out_dir, 
    out_name=out_name
)


# 4. Plot UMAP feature plots for annotated cell type markers
if (len(umap_markers_list) > 0):
    # Call the dedicated function with the new list of consensus markers
    plot_umap_feature_plots(
      adata = adata_full, 
      markers_dict = umap_markers_list, 
      out_dir = out_dir, 
      out_name = f"{out_name}_scina_cell_type"
    )
else:
    print("No significant positive markers found for consensus clusters. Skipping UMAP FeaturePlot generation based on consensus labels.")


# 5. Plot umap with metadata overlay
# drop unnamed columns from metadata
adata_annotated.obs = adata_annotated.obs.loc[:, ~adata_annotated.obs.columns.astype(str).str.contains('^Unnamed')]

# Use a set of columns to plot UMAP for
metadata_cols_to_plot = [
    col for col in adata_annotated.obs.columns 
    if col not in ['barcode', 'cluster', 'cell_type_ID'] # Exclude 'cell_type_ID' since it's custom plotted
]

for col in metadata_cols_to_plot:
    # Use standard scanpy plotting for all other metadata columns
    sc.pl.umap(adata_annotated, color=[col], 
               legend_loc='right margin', title='', frameon=False, ncols=1, 
               save=f"_{out_name}_{col}.png")


# 6. Generate blue-print csv file for manual cell type annotation
annotation_df = adata_annotated.obs.drop_duplicates(subset=["leiden"])[["leiden", "cell_type"]].copy()
annotation_df.columns = ["cluster", "cell_type"]
annotation_df['cluster'] = pd.to_numeric(annotation_df['cluster'])
annotation_df = annotation_df.sort_values(by="cluster")
annotation_df['order'] = annotation_df['cell_type'].cat.codes
annotation_df.to_csv(os.path.join(out_dir, f"{out_name}_annotation.csv"), index=False)

# 7. Export metadata as csv
adata_annotated.obs.to_csv(os.path.join(out_dir, f"{out_name}_obs.csv"))

# 8. Save final object
out_file = os.path.join(out_dir, f"{out_name}_scRNAseq_no_doublets_annotated.h5ad")
adata_annotated.write(out_file, compression='gzip')

print(f"\nFinal object saved to {out_file}")