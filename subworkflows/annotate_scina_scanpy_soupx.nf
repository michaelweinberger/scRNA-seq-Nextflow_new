#!/usr/bin/env nextflow



// Enable DSL2 for more powerful functionality
nextflow.enable.dsl = 2



// Processes

/*
 * convert Scanpy object into Seurat object
 */
process CONVERT_SCANPY_SEURAT_PR {
    debug true
    publishDir "${outdir}/scanpy/scina", pattern: "", mode: "copy", saveAs: { filename -> "${filename}" }

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

    #[ ! -d "${outdir}/scanpy/scina" ] && mkdir -p "${outdir}/scanpy/scina"

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
 * run SCINA annotation
 */
process SCINA_ANNOTATION_PR {
    debug true
    publishDir "${outdir}/scanpy/scina", pattern: "", mode: "copy", saveAs: { filename -> "${filename}" }

    container { ( params.docker_enabled == "true" || params.docker_enabled == true ) ? "michaelweinberger/seurat-soupx-scina:v1" : "" }

    input:
    path ( seurat_object )
    path ( cell_type_markers_csv )
    val  ( outdir )
    val  ( out_name )
    val  ( docker_enabled )
    val  ( r_module )

    output:
    path ( "*.pdf"                                      ), emit: plots_pdf
    path ( "${out_name}_scina_celltype_correlation.csv" ), emit: scina_correlations_csv
    path ( "${out_name}_scina_predictions.csv"          ), emit: scina_results_csv
    path ( "versions.txt"                               ), emit: versions

    script:
    """
    if [ "$docker_enabled" == "false" ] ; then
        module load "$r_module"
    fi

    #[ ! -d "${outdir}/scanpy/scina" ] && mkdir -p "${outdir}/scanpy/scina"

    4_4_scina_annotation.R --args seurat_obj_file="$seurat_object" \
                                  cell_type_markers_csv="$cell_type_markers_csv" \
                                  out_dir="\$PWD" \
                                  out_name="$out_name"

    echo "${task.process}:" > versions.txt
        R --version | head -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(Seurat); print(paste('Seurat', packageVersion('Seurat')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(readr); print(paste('readr', packageVersion('readr')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(SCINA); print(paste('SCINA', packageVersion('SCINA')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(pheatmap); print(paste('pheatmap', packageVersion('pheatmap')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(dplyr); print(paste('dplyr', packageVersion('dplyr')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(ggplot2); print(paste('ggplot2', packageVersion('ggplot2')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
    """
}



/*
 * integrate SCINA results into Scanpy object
 */
process SCINA_SCANPY_INTEGRATION_PR {
    debug true
    publishDir "${outdir}/scanpy/scina", pattern: "", mode: "copy", saveAs: { filename -> "${filename}" }

    container { ( params.docker_enabled == "true" || params.docker_enabled == true ) ? "michaelweinberger/python-3.11.14-scanpy:v1" : "" }

    input:
    path ( scrnaseq_object )
    path ( scrnaseq_object_raw )
    path ( scina_results_csv )
    path ( cell_type_markers_csv )
    path ( cell_type_metadata_csv )
    val  ( leiden_res )
    val  ( outdir )
    val  ( out_name )
    val  ( docker_enabled )
    val  ( python_module )

    output:
    path ( "*.pdf"                                           ), emit: plots_pdf
    path ( "*.png"                                           ), emit: plots_png
    path ( "*.csv"                                           ), emit: csv_files
    path ( "*.xlsx"                                          ), emit: excel_files
    path ( "${out_name}_scRNAseq_no_doublets_annotated.h5ad" ), emit: scrnaseq_object
    path ( "versions.txt"                                    ), emit: versions

    script:
    """
    if [ "$docker_enabled" == "false" ] ; then
        module load "$python_module"
    fi

    #[ ! -d "${outdir}/scanpy/scina" ] && mkdir -p "${outdir}/scanpy/scina"

    4_4_scina_scanpy_integration.py -i "$scrnaseq_object" \
                                    -j "$scrnaseq_object_raw" \
                                    -s "$scina_results_csv" \
                                    -m "$cell_type_markers_csv" \
                                    -t "$cell_type_metadata_csv" \
                                    -r "$leiden_res" \
                                    -o "\$PWD" \
                                    -n "$out_name"

    echo "${task.process}:" > versions.txt
        python --version | sed -e \$'s/^/\t/' >> versions.txt
        python -c "import scanpy; print(f'scanpy,{scanpy.__version__}')" | sed -e \$'s/^/\t/' >> versions.txt
        python -c "import pandas; print(f'pandas,{pandas.__version__}')" | sed -e \$'s/^/\t/' >> versions.txt
        python -c "import numpy; print(f'numpy,{numpy.__version__}')" | sed -e \$'s/^/\t/' >> versions.txt
        python -c "import matplotlib; print(f'matplotlib,{matplotlib.__version__}')" | sed -e \$'s/^/\t/' >> versions.txt
        python -c "import anndata; print(f'anndata,{anndata.__version__}')" | sed -e \$'s/^/\t/' >> versions.txt
        python -c "import scipy; print(f'scipy,{scipy.__version__}')" | sed -e \$'s/^/\t/' >> versions.txt
    """
}



// Workflow

workflow ANNOTATE_SCINA_SCANPY_SOUPX_WF {

    take:
        scrnaseq_object_raw
        cell_metadata_csv
        scrnaseq_object
        cell_type_markers_csv
        cell_type_metadata_csv
        leiden_res
        outdir
        out_name
        docker_enabled
        r_module
        python_module

    main:
        ch_versions = Channel.empty()

        CONVERT_SCANPY_SEURAT_PR (
            scrnaseq_object_raw,
            cell_metadata_csv,
            outdir,
            out_name,
            docker_enabled,
            r_module
        )
        ch_versions = ch_versions.mix(CONVERT_SCANPY_SEURAT_PR.out.versions)

        SCINA_ANNOTATION_PR (
            CONVERT_SCANPY_SEURAT_PR.out.seurat_object,
            cell_type_markers_csv,
            outdir,
            out_name,
            docker_enabled,
            r_module
        )
        ch_versions = ch_versions.mix(SCINA_ANNOTATION_PR.out.versions)

        SCINA_SCANPY_INTEGRATION_PR (
            scrnaseq_object,
            scrnaseq_object_raw,
            SCINA_ANNOTATION_PR.out.scina_results_csv,
            cell_type_markers_csv,
            cell_type_metadata_csv,
            leiden_res,
            outdir,
            out_name,
            docker_enabled,
            python_module
        )
        ch_versions = ch_versions.mix(SCINA_SCANPY_INTEGRATION_PR.out.versions)

    emit:
        versions        = ch_versions
        scrnaseq_object = SCINA_SCANPY_INTEGRATION_PR.out.scrnaseq_object
}
