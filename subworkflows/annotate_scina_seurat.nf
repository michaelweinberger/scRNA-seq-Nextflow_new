#!/usr/bin/env nextflow



// Enable DSL2 for more powerful functionality
nextflow.enable.dsl = 2



// Processes

/*
 * run SCINA annotation + integrate results into Seurat object
 */
process SCINA_ANNOTATION_INTEGRATION_PR {
    debug true
    publishDir "${outdir}/seurat/scina", pattern: "", mode: "copy", saveAs: { filename -> "${filename}" }

    container { ( params.docker_enabled == "true" || params.docker_enabled == true ) ? "michaelweinberger/seurat-soupx-scina:v1" : "" }

    input:
    path ( scrnaseq_object )
    path ( cell_type_markers_csv )
    path ( cell_type_metadata_csv )
    val  ( leiden_res )
    val  ( outdir )
    val  ( out_name )
    val  ( docker_enabled )
    val  ( r_module )

    output:
    path ( "*.pdf"                                          ), emit: plots_pdf
    path ( "*.png"                                          ), emit: plots_png
    path ( "*.csv"                                          ), emit: csv_files
    path ( "*.xlsx"                                         ), emit: excel_files
    path ( "${out_name}_metadata.csv"                       ), emit: cell_metadata_csv
    path ( "${out_name}_scRNAseq_no_doublets_annotated.rds" ), emit: scrnaseq_object
    path ( "versions.txt"                                   ), emit: versions

    script:
    """
    if [ "$docker_enabled" == "false" ] ; then
        module load "$r_module"
    fi

    #[ ! -d "${outdir}/seurat/scina" ] && mkdir -p "${outdir}/seurat/scina"

    4_4_scina_annotation_integration.R --args seurat_obj_file="$scrnaseq_object" \
                                              cell_type_markers_csv="$cell_type_markers_csv" \
                                              cell_type_metadata_csv="$cell_type_metadata_csv" \
                                              leiden_res="$leiden_res" \
                                              out_dir="\$PWD" \
                                              out_name="$out_name"

    echo "${task.process}:" > versions.txt
        R --version | head -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(Seurat); print(paste('Seurat', packageVersion('Seurat')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(readr); print(paste('readr', packageVersion('readr')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(SCINA); print(paste('SCINA', packageVersion('SCINA')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(pheatmap); print(paste('pheatmap', packageVersion('pheatmap')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(RColorBrewer); print(paste('RColorBrewer', packageVersion('RColorBrewer')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(dplyr); print(paste('dplyr', packageVersion('dplyr')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(ggplot2); print(paste('ggplot2', packageVersion('ggplot2')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(writexl); print(paste('writexl', packageVersion('writexl')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
    """
}



// Workflow

workflow ANNOTATE_SCINA_SEURAT_WF {

    take:
        scrnaseq_object
        cell_type_markers_csv
        cell_type_metadata_csv
        leiden_res
        outdir
        out_name
        docker_enabled
        r_module

    main:
        ch_versions = Channel.empty()

        SCINA_ANNOTATION_INTEGRATION_PR (
            scrnaseq_object,
            cell_type_markers_csv,
            cell_type_metadata_csv,
            leiden_res,
            outdir,
            out_name,
            docker_enabled,
            r_module
        )
        ch_versions = ch_versions.mix(SCINA_ANNOTATION_INTEGRATION_PR.out.versions)

    emit:
        versions        = ch_versions
        scrnaseq_object = SCINA_ANNOTATION_INTEGRATION_PR.out.scrnaseq_object
}