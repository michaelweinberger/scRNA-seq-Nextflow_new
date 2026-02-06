#!/usr/bin/env nextflow



// Enable DSL2 for more powerful functionality
nextflow.enable.dsl = 2



// Processes

/*
 * run Scanpy clustering
 */
process CLUSTER_SCANPY_PR {
    debug true
    publishDir "${outdir}/scanpy", mode: "copy"

    container { ( params.docker_enabled == "true" || params.docker_enabled == true ) ? "michaelweinberger/python-3.11.14-scanpy:v1" : "" }

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
    val  ( python_module )
    
    output:
    path ( "${out_name}_obs.csv"                            ), emit: cell_metadata_csv
    path ( "${out_name}.h5ad"                               ), emit: scrnaseq_object_raw
    path ( "versions.txt"                                   ), emit: versions

    script:
    """
    if [ "$docker_enabled" == "false" ] ; then
        module load "$python_module"
    fi

    4_3_scanpy.py -i "$cellranger_out_dir" \
                  -m "$metadata" \
                  -gtf "$cellranger_gtf_file" \
                  -b "$transcriptome_blast_txt" \
                  -ccg "$cellcycle_genes_csv" \
                  -o "\$PWD" \
                  -n "$out_name" \
                  -d "$drop_ribo_prot" \
                  -g "$genes_drop_csv" \
                  -mig "$min_genes" \
                  -mag "$max_genes" \
                  -mam "$max_perc_mt" \
                  -mic "$min_cells" \
                  -npc "$n_pcs" \
                  -hv "$harmony_var" \
                  -r "$leiden_res"

    echo "${task.process}:" > versions.txt
        python --version | sed -e \$'s/^/\t/' >> versions.txt
        python -c "import scanpy; print(f'scanpy,{scanpy.__version__}')" | sed -e \$'s/^/\t/' >> versions.txt
        python -c "import anndata; print(f'anndata,{anndata.__version__}')" | sed -e \$'s/^/\t/' >> versions.txt
        python -c "import numpy; print(f'numpy,{numpy.__version__}')" | sed -e \$'s/^/\t/' >> versions.txt
        python -c "import pandas; print(f'pandas,{pandas.__version__}')" | sed -e \$'s/^/\t/' >> versions.txt

    # add delay to allow for HPC sync
    sleep 10
    """
}



/*
 * convert Scanpy object into Seurat object
 */
process CONVERT_SCANPY_SEURAT_PR {
    debug true
    publishDir "${outdir}/scanpy", mode: "copy"

    container { ( params.docker_enabled == "true" || params.docker_enabled == true ) ? "michaelweinberger/seurat-soupx-scina:v1" : "" }

    input:
    path ( scrnaseq_object_raw )
    path ( cell_metadata_csv )
    val  ( outdir )
    val  ( out_name )
    val  ( docker_enabled )
    val  ( r_module )

    output:
    path ( "${out_name}.rds" ), emit: seurat_object
    path ( "versions.txt"    ), emit: versions

    script:
    """
    if [ "$docker_enabled" == "false" ] ; then
        module load "$r_module"
    fi

    4_4_scanpy_to_seurat.R --args scanpy_obj="$scrnaseq_object_raw" \
                                  cell_metadata_csv="$cell_metadata_csv" \
                                  out_dir="\$PWD" \
                                  out_name="$out_name"

    echo "${task.process}:" > versions.txt
        R --version | head -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(Seurat); print(paste('Seurat', packageVersion('Seurat')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(anndataR); print(paste('anndataR', packageVersion('anndataR')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(hdf5r); print(paste('hdf5r', packageVersion('hdf5r')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
    """
}



/*
 * correct for ambient RNA contamination
 */
process AMBIENT_RNA_SOUPX_PR {
    debug true
    publishDir "${outdir}/soupx", mode: "copy"

    container { ( params.docker_enabled == "true" || params.docker_enabled == true ) ? "michaelweinberger/seurat-soupx-scina:v1" : "" }

    input:
    val  ( log_mapping_mode )
    path ( seurat_object_file )
    path ( cellranger_count_info )
    val  ( harmony_var )
    val  ( outdir )
    val  ( out_name )
    val  ( docker_enabled )
    val  ( r_module )

    output:
    path ( "*${out_name}*.pdf"             ), emit: plots_pdf
    path ( "*${out_name}*.csv"             ), emit: csv_files
    path ( "*${out_name}*.xlsx"            ), emit: xlsx_files
    path ( "${out_name}_cleaned_10X"       ), emit: soupx_dir
    path ( "versions.txt"                  ), emit: versions

    script:
    """
    if [ "$docker_enabled" == "false" ] ; then
        module load "$r_module"
    fi

    4_soupx.R --args mapping_mode="$log_mapping_mode" \
                     seurat_obj_file="$seurat_object_file" \
                     cellranger_count_info="$cellranger_count_info" \
                     harmony_var="$harmony_var" \
                     out_dir="\$PWD" \
                     out_name="${out_name}"

    echo "${task.process}:" > versions.txt
        R --version | head -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(SoupX); print(paste('SoupX', packageVersion('SoupX')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(DropletUtils); print(paste('DropletUtils', packageVersion('DropletUtils')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(Seurat); print(paste('Seurat', packageVersion('Seurat')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(Matrix); print(paste('Matrix', packageVersion('Matrix')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(ggplot2); print(paste('ggplot2', packageVersion('ggplot2')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(dplyr); print(paste('dplyr', packageVersion('dplyr')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
    """
}



/*
 * run Scanpy clustering following Soupx ambient RNA correction
 */
process CLUSTER_SCANPY_SOUPX_PR {
    debug true
    publishDir "${outdir}/scanpy", mode: "copy"

    container { ( params.docker_enabled == "true" || params.docker_enabled == true ) ? "michaelweinberger/python-3.11.14-scanpy:v1" : "" }

    input:
    path ( cellranger_out_dir )
    path ( metadata )
    path ( cellranger_gtf_file )
    path ( soupx_dir )
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
    val  ( python_module )
    
    output:
    path ( "*${out_name}*.pdf"                                    ), emit: plots_pdf
    path ( "*${out_name}*.png"                                    ), emit: plots_png
    path ( "${out_name}_soupx_markers.csv"                        ), emit: markers_csv
    path ( "${out_name}_soupx_obs.csv"                            ), emit: cell_metadata_csv
    path ( "${out_name}_soupx_var.csv"                            ), emit: feature_metadata_csv
    path ( "${out_name}_soupx_scRNAseq_analysed_no_doublets.h5ad" ), emit: scrnaseq_object
    path ( "${out_name}_soupx_doublets_detected.h5ad"             ), emit: scrnaseq_object_raw_1
    path ( "${out_name}_soupx.h5ad"                               ), emit: scrnaseq_object_raw
    path ( "versions.txt"                                         ), emit: versions

    script:
    """
    if [ "$docker_enabled" == "false" ] ; then
        module load "$python_module"
    fi

    4_3_scanpy.py -i "$cellranger_out_dir" \
                  -m "$metadata" \
                  -gtf "$cellranger_gtf_file" \
                  -s "$soupx_dir" \
                  -b "$transcriptome_blast_txt" \
                  -ccg "$cellcycle_genes_csv" \
                  -o "\$PWD" \
                  -n "${out_name}_soupx" \
                  -d "$drop_ribo_prot" \
                  -g "$genes_drop_csv" \
                  -mig "$min_genes" \
                  -mag "$max_genes" \
                  -mam "$max_perc_mt" \
                  -mic "$min_cells" \
                  -npc "$n_pcs" \
                  -hv "$harmony_var" \
                  -r "$leiden_res"

    echo "${task.process}:" > versions.txt
        python --version | sed -e \$'s/^/\t/' >> versions.txt
        python -c "import scanpy; print(f'scanpy,{scanpy.__version__}')" | sed -e \$'s/^/\t/' >> versions.txt
        python -c "import anndata; print(f'anndata,{anndata.__version__}')" | sed -e \$'s/^/\t/' >> versions.txt
        python -c "import numpy; print(f'numpy,{numpy.__version__}')" | sed -e \$'s/^/\t/' >> versions.txt
        python -c "import pandas; print(f'pandas,{pandas.__version__}')" | sed -e \$'s/^/\t/' >> versions.txt

    # add delay to allow for HPC sync
    sleep 10
    """
}



// Workflow

workflow CLUSTER_SCANPY_SOUPX_WF {

    take:
        log_mapping_mode
        cellranger_out_dir
        metadata
        cellranger_gtf_file
        cellranger_count_info
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
        python_module
        r_module

    main:
        ch_versions = Channel.empty()

        CLUSTER_SCANPY_PR ( 
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
            python_module,
        )
        ch_versions = ch_versions.mix(CLUSTER_SCANPY_PR.out.versions)

        CONVERT_SCANPY_SEURAT_PR (
            CLUSTER_SCANPY_PR.out.scrnaseq_object_raw,
            CLUSTER_SCANPY_PR.out.cell_metadata_csv,
            outdir,
            out_name,
            docker_enabled,
            r_module,
        )
        ch_versions = ch_versions.mix(CONVERT_SCANPY_SEURAT_PR.out.versions)

        AMBIENT_RNA_SOUPX_PR (
            log_mapping_mode,
            CONVERT_SCANPY_SEURAT_PR.out.seurat_object,
            cellranger_count_info,
            harmony_var,
            outdir,
            out_name,
            docker_enabled,
            r_module,
        )
        ch_versions = ch_versions.mix(AMBIENT_RNA_SOUPX_PR.out.versions)

        CLUSTER_SCANPY_SOUPX_PR ( 
            cellranger_out_dir,
            metadata,
            cellranger_gtf_file,
            AMBIENT_RNA_SOUPX_PR.out.soupx_dir,
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
            python_module,
        )
        ch_versions = ch_versions.mix(CLUSTER_SCANPY_SOUPX_PR.out.versions)

    emit:
        versions             = ch_versions
        cell_metadata_csv    = CLUSTER_SCANPY_SOUPX_PR.out.cell_metadata_csv
        scrnaseq_object_raw  = CLUSTER_SCANPY_SOUPX_PR.out.scrnaseq_object_raw
        scrnaseq_object      = CLUSTER_SCANPY_SOUPX_PR.out.scrnaseq_object
}
