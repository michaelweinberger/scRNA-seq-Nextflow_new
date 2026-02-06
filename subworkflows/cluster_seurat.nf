#!/usr/bin/env nextflow



// Enable DSL2 for more powerful functionality
nextflow.enable.dsl = 2



// Processes

/*
 * run Seurat clustering
 */
process CLUSTER_SEURAT_PR {
    debug true
    publishDir "${outdir}/seurat", pattern: "", mode: "copy", saveAs: { filename -> "${filename}" }

    container { ( params.docker_enabled == "true" || params.docker_enabled == true ) ? "michaelweinberger/seurat-soupx-scina:v1" : "" }

    input:
    path ( cellranger_out_dir )
    path ( metadata )
    path ( cellranger_gtf_file )
    path ( transcriptome_blast_txt )
    path ( cellcycle_genes_csv )
    val  ( out_name )
    val  ( drop_ribo_prot )
    path ( genes_drop_csv )
    val  ( min_genes )
    val  ( max_genes )
    val  ( max_perc_mt )
    val  ( min_cells )
    val  ( n_pcs )
    val  ( harmony_var )
    val  ( leiden_res )
    val  ( outdir )
    val  ( docker_enabled )
    val  ( r_module )
    
    output:
    path ( "*${out_name}*.pdf"                             ), emit: plots_pdf
    path ( "*${out_name}*.png"                             ), emit: plots_png
    path ( "${out_name}_doublets.csv"                      ), emit: doublets_csv
    path ( "${out_name}_markers.xlsx"                      ), emit: markers_csv
    path ( "${out_name}_cell_metadata.csv"                 ), emit: cell_metadata_csv
    path ( "${out_name}_feature_metadata.csv"              ), emit: feature_metadata_csv
    path ( "${out_name}_scRNAseq_analysed_no_doublets.rds" ), emit: scrnaseq_object
    path ( "versions.txt"                                  ), emit: versions

    script:
    """
    if [ "$docker_enabled" == "false" ] ; then
        module load "$r_module"
    fi

    4_3_seurat.R --args data_dir="$cellranger_out_dir" \
                        metadata_file="$metadata" \
                        gtf_file="$cellranger_gtf_file" \
                        blast_file="$transcriptome_blast_txt" \
                        cellcycle_genes_csv="$cellcycle_genes_csv" \
                        out_dir="\$PWD" \
                        project_name="$out_name" \
                        drop_ribo_prot="$drop_ribo_prot" \
                        genes_drop_csv="$genes_drop_csv" \
                        min_genes="$min_genes" \
                        max_genes="$max_genes" \
                        max_perc_mt="$max_perc_mt" \
                        min_cells="$min_cells" \
                        n_pcs="$n_pcs" \
                        harmony_var="$harmony_var" \
                        leiden_res="$leiden_res"

    echo "${task.process}:" > versions.txt
        R --version | head -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(dplyr); print(paste('dplyr', packageVersion('dplyr')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(ggplot2); print(paste('ggplot2', packageVersion('ggplot2')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(patchwork); print(paste('patchwork', packageVersion('patchwork')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(pheatmap); print(paste('pheatmap', packageVersion('pheatmap')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(Seurat); print(paste('Seurat', packageVersion('Seurat')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(viridis); print(paste('viridis', packageVersion('viridis')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(writexl); print(paste('writexl', packageVersion('writexl')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
    """
}



// Workflow

workflow CLUSTER_SEURAT_WF {

    take:
        cellranger_out_dir
        metadata
        cellranger_gtf_file
        transcriptome_blast_txt
        cellcycle_genes_csv
        out_name
        drop_ribo_prot
        genes_drop_csv
        min_genes
        max_genes
        max_perc_mt
        min_cells
        n_pcs
        harmony_var
        leiden_res
        outdir
        docker_enabled
        r_module

    main:
        ch_versions = Channel.empty()

        CLUSTER_SEURAT_PR ( 
            cellranger_out_dir,
            metadata,
            cellranger_gtf_file,
            transcriptome_blast_txt,
            cellcycle_genes_csv,
            out_name,
            drop_ribo_prot,
            genes_drop_csv,
            min_genes,
            max_genes,
            max_perc_mt,
            min_cells,
            n_pcs,
            harmony_var,
            leiden_res,
            outdir,
            docker_enabled,
            r_module,
        )
        ch_versions = ch_versions.mix(CLUSTER_SEURAT_PR.out.versions)

    emit:
        versions          = ch_versions
        scrnaseq_object   = CLUSTER_SEURAT_PR.out.scrnaseq_object
}
