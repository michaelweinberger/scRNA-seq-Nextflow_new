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
    path ( "*${out_name}*.pdf"                              ), emit: plots_pdf
    path ( "*${out_name}*.png"                              ), emit: plots_png
    path ( "${out_name}_markers.csv"                        ), emit: markers_csv
    path ( "${out_name}_obs.csv"                            ), emit: cell_metadata_csv
    path ( "${out_name}_var.csv"                            ), emit: feature_metadata_csv
    path ( "${out_name}_scRNAseq_analysed_no_doublets.h5ad" ), emit: scrnaseq_object
    path ( "${out_name}_doublets_detected.h5ad"             ), emit: scrnaseq_object_raw_1
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



// Workflow

workflow CLUSTER_SCANPY_WF {

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
        python_module

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

    emit:
        versions            = ch_versions
        cell_metadata_csv   = CLUSTER_SCANPY_PR.out.cell_metadata_csv
        scrnaseq_object_raw = CLUSTER_SCANPY_PR.out.scrnaseq_object_raw
        scrnaseq_object     = CLUSTER_SCANPY_PR.out.scrnaseq_object
}
