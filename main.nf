#!/usr/bin/env nextflow



// Enable DSL2 for more powerful functionality
nextflow.enable.dsl = 2



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PIPELINE_INIT_WF                } from "./subworkflows/pipeline_init.nf"
include { CELLRANGER_REF_WF               } from "./subworkflows/cellranger_ref.nf"
include { CELLRANGER_COUNT_AGGR_WF        } from "./subworkflows/cellranger_count_aggr.nf"
include { CELLRANGER_MULTI_AGGR_WF        } from "./subworkflows/cellranger_multi_aggr.nf"
include { TRANSCRIPTOME_PROTEOME_BLAST_WF } from "./subworkflows/transcriptome_proteome_blast.nf"
include { CLUSTER_SCANPY_WF               } from "./subworkflows/cluster_scanpy.nf"
include { CLUSTER_SCANPY_SOUPX_WF         } from "./subworkflows/cluster_scanpy_soupx.nf"
include { CLUSTER_SEURAT_WF               } from "./subworkflows/cluster_seurat.nf"
include { CLUSTER_SEURAT_SOUPX_WF         } from "./subworkflows/cluster_seurat_soupx.nf"
include { ANNOTATE_SCINA_SCANPY_WF        } from "./subworkflows/annotate_scina_scanpy.nf"
include { ANNOTATE_SCINA_SCANPY_SOUPX_WF  } from "./subworkflows/annotate_scina_scanpy_soupx.nf"
include { ANNOTATE_SCINA_SEURAT_WF        } from "./subworkflows/annotate_scina_seurat.nf"
include { ANNOTATE_SCINA_SEURAT_SOUPX_WF  } from "./subworkflows/annotate_scina_seurat_soupx.nf"



// Workflow

workflow {

    main:

        // Create output directory if it does not exist
        if ( params.outdir ) {
            outdir = params.outdir
        } else {
            outdir = workflow.projectDir / "out"
        }
        new File ( outdir ).mkdirs()

        // Set output directory for clustering with ambient RNA correction
        if ( params.run_soupx == "Yes" ) {
            outdir_soupx = outdir + "/soupx"
            out_name_soupx = params.project + "_soupx"
        }

        // Set optional parameters
        if ( params.genes_drop_csv ) {
            genes_drop_csv = params.genes_drop_csv
        } else {
            def empty_genes_drop_csv = file("${outdir}/genes_drop_csv_empty.csv")
            if( !empty_genes_drop_csv.exists() ) { empty_genes_drop_csv.text = "" }
            genes_drop_csv = Channel.value(empty_genes_drop_csv)
        }

        // print logging info
        if ( params.input_count && params.input_multi ) {
            error "Parameters 'input_count' and 'input_multi' cannot both be set."
        }

        if ( params.input_count ) {
            log_mapping_mode = "cell ranger count"
            log_input = params.input_count
        } else if ( params.input_multi ) {
            log_mapping_mode = "cell ranger multi"
            log_input = params.input_multi
        } else {
            log_mapping_mode = params.mapping_mode
            log_input = "Mapping skipped"
        }

        log.info """\
            ${params.manifest.name} v${params.manifest.version}
            ============================================================
            project name                : ${params.project}
            input from                  : ${log_input}
            output to                   : ${outdir}
            mapping via                 : ${log_mapping_mode}
            clustering via              : ${params.clustering_mode}
            species                     : ${params.species}
            10X chemistry               : ${params.chem}
            ------------------------------------------------------------
            run as                      : ${workflow.commandLine}
            started at                  : ${workflow.start}
            config files                : ${workflow.configFiles}
            container                   : ${workflow.containerEngine}:${workflow.container}
            ============================================================
            """
         .stripIndent()


        // Create channel to collect software versions
        ch_versions = Channel.empty()


        //
        // SUBWORKFLOW: Generate Cell Ranger reference genome index
        //
        CELLRANGER_REF_WF ( 
            params.species,
            params.genome,
            params.genome_ucsc,
            params.species_latin,
            params.ensembl_version,
            outdir,
            params.docker_enabled,
            params.cellranger_module,
        )
        ch_versions = ch_versions.unique().mix(CELLRANGER_REF_WF.out.versions)


        //
        // SUBWORKFLOW: Process input sample sheet
        //
        if ( params.input_count ) {
            PIPELINE_INIT_WF (
                params.input_count,
            )
        } else if ( params.input_multi ) {
            PIPELINE_INIT_WF (
                params.input_multi,
            )
        }


        //
        // SUBWORKFLOW: Align to reference genome with Cell Ranger
        //
        if ( params.input_count ) {
            CELLRANGER_COUNT_AGGR_WF (
                CELLRANGER_REF_WF.out.cellranger_index,
                params.input_count,
                PIPELINE_INIT_WF.out.input,
                params.project,
                outdir,
                params.docker_enabled,
                params.cellranger_module,
                params.chem,
            )
            ch_versions = ch_versions.unique().mix(CELLRANGER_COUNT_AGGR_WF.out.versions)
        } else if ( params.input_multi ) {
            CELLRANGER_MULTI_AGGR_WF (
                CELLRANGER_REF_WF.out.cellranger_ref_dir,
                params.input_multi,
                PIPELINE_INIT_WF.out.input,
                params.project,
                outdir,
                params.docker_enabled,
                params.cellranger_module,
            )
            ch_versions = ch_versions.unique().mix(CELLRANGER_MULTI_AGGR_WF.out.versions)
        } else {
            println "No input found for Cell Ranger Count or Multi: Skipping Cell Ranger alignment."
        }

        // Set parameters when skipping mapping
        if ( !params.input_count && !params.input_multi && !params.cellranger_out_dir ) {
            error "File not found - parameter 'cellranger_out_dir' needs to be set if both 'input_count' and 'input_multi' are not set."
        } else if ( !params.input_count && !params.input_multi && params.cellranger_out_dir ) {
            cellranger_out_dir = params.cellranger_out_dir

            if ( params.metadata ) {
                cellranger_metadata = params.metadata
            } else {
                error "File not found - parameter 'metadata' needs to be set if both 'input_count' and 'input_multi' are not set."
            }

            if ( params.cellranger_info_tsv ) {
                cellranger_count_info = params.cellranger_info_tsv
                println "Pipeline runs ambient RNA correction."
            } else {
                println "Pipeline does not run ambient RNA correction."
            }
        } else if ( ( params.input_count || params.input_multi ) && params.cellranger_out_dir ) {
            error "'input_count' or 'input_multi' cannot be set at the same time as 'cellranger_out_dir'."
        } else if ( params.input_count && !params.input_multi && !params.cellranger_out_dir ) {
            cellranger_out_dir    = CELLRANGER_COUNT_AGGR_WF.out.cellranger_aggr_bc_matrix
            cellranger_metadata   = CELLRANGER_COUNT_AGGR_WF.out.cellranger_metadata
            cellranger_count_info = CELLRANGER_COUNT_AGGR_WF.out.cellranger_count_info
        } else if ( !params.input_count && params.input_multi && !params.cellranger_out_dir ) {
            cellranger_out_dir    = CELLRANGER_MULTI_AGGR_WF.out.cellranger_aggr_bc_matrix
            cellranger_metadata   = CELLRANGER_MULTI_AGGR_WF.out.cellranger_metadata
            cellranger_count_info = CELLRANGER_MULTI_AGGR_WF.out.cellranger_count_info
        }


        //
        // SUBWORKFLOW: Blast transcriptome/proteome of species set to human transcriptome/proteome
        //
        transcriptome_blast_txt = Channel.empty()

        if ( params.species != "human" ) {
            TRANSCRIPTOME_PROTEOME_BLAST_WF (
                params.species,
                params.biomart_dataset,
                params.biomart_genome,
                params.chrom_number,
                params.human_biomart_dataset,
                params.human_biomart_genome,
                params.uniprot_id,
                params.ncbi_id,
                outdir,
                params.docker_enabled,
                params.r_module,
                params.blast_module,
            )
            ch_versions = ch_versions.unique().mix(TRANSCRIPTOME_PROTEOME_BLAST_WF.out.versions)
            transcriptome_blast_txt = TRANSCRIPTOME_PROTEOME_BLAST_WF.out.transcriptome_blast_txt
        } else {
            def empty_file = file("${outdir}/blast_empty.txt")
            if( !empty_file.exists() ) { empty_file.text = "" }
            transcriptome_blast_txt = Channel.value(empty_file)
        }

        //
        // SUBWORKFLOW: Cluster cells with Scanpy or Seurat
        //
        if ( params.clustering_mode == "scanpy" ) {
            CLUSTER_SCANPY_WF (
                cellranger_out_dir,
                cellranger_metadata,
                CELLRANGER_REF_WF.out.cellranger_gtf_file,
                transcriptome_blast_txt,
                params.cellcycle_genes_csv,
                params.project,
                params.drop_ribo_prot,
                genes_drop_csv,
                params.min_genes,
                params.max_genes,
                params.max_perc_mt,
                params.min_cells,
                params.n_pcs,
                params.harmony_var,
                params.leiden_res,
                outdir,
                params.docker_enabled,
                params.python_module,
            )
            ch_versions = ch_versions.unique().mix(CLUSTER_SCANPY_WF.out.versions)
        } else if ( params.clustering_mode == "seurat" ) {
            CLUSTER_SEURAT_WF (
                cellranger_out_dir,
                cellranger_metadata,
                CELLRANGER_REF_WF.out.cellranger_gtf_file,
                transcriptome_blast_txt,
                params.cellcycle_genes_csv,
                params.project,
                params.drop_ribo_prot,
                genes_drop_csv,
                params.min_genes,
                params.max_genes,
                params.max_perc_mt,
                params.min_cells,
                params.n_pcs,
                params.harmony_var,
                params.leiden_res,
                outdir,
                params.docker_enabled,
                params.r_module,
            )
            ch_versions = ch_versions.unique().mix(CLUSTER_SEURAT_WF.out.versions)
        } else {
            error "'clustering_mode' parameter needs to be: 'scanpy' or 'seurat'"
        }

        // Ambient RNA correction + clustering
        if ( params.run_soupx == "Yes" && cellranger_count_info ) {
            if ( params.clustering_mode == "scanpy" ) {
                CLUSTER_SCANPY_SOUPX_WF (
                    log_mapping_mode,
                    cellranger_out_dir,
                    cellranger_metadata,
                    CELLRANGER_REF_WF.out.cellranger_gtf_file,
                    cellranger_count_info,
                    transcriptome_blast_txt,
                    params.cellcycle_genes_csv,
                    out_name_soupx,
                    params.drop_ribo_prot,
                    genes_drop_csv,
                    params.min_genes,
                    params.max_genes,
                    params.max_perc_mt,
                    params.min_cells,
                    params.n_pcs,
                    params.harmony_var,
                    params.leiden_res,
                    outdir_soupx,
                    params.docker_enabled,
                    params.python_module,
                    params.r_module,
                )
                ch_versions = ch_versions.unique().mix(CLUSTER_SCANPY_SOUPX_WF.out.versions)
            } else if ( params.clustering_mode == "seurat" ) {
                CLUSTER_SEURAT_SOUPX_WF (
                    log_mapping_mode,
                    cellranger_out_dir,
                    cellranger_metadata,
                    CELLRANGER_REF_WF.out.cellranger_gtf_file,
                    cellranger_count_info,
                    transcriptome_blast_txt,
                    params.cellcycle_genes_csv,
                    out_name_soupx,
                    params.drop_ribo_prot,
                    genes_drop_csv,
                    params.min_genes,
                    params.max_genes,
                    params.max_perc_mt,
                    params.min_cells,
                    params.n_pcs,
                    params.harmony_var,
                    params.leiden_res,
                    outdir_soupx,
                    params.docker_enabled,
                    params.r_module,
                )
                ch_versions = ch_versions.unique().mix(CLUSTER_SEURAT_SOUPX_WF.out.versions)
            } else {
                error "'clustering_mode' parameter needs to be: 'scanpy' or 'seurat'"
            }
        }


        //
        // SUBWORKFLOW: Automatically annotate cell types in scRNA-seq object
        //
        if ( params.clustering_mode == "scanpy" ) {
            ANNOTATE_SCINA_SCANPY_WF (
                CLUSTER_SCANPY_WF.out.scrnaseq_object_raw,
                CLUSTER_SCANPY_WF.out.cell_metadata_csv,
                CLUSTER_SCANPY_WF.out.scrnaseq_object,
                params.cell_type_markers_csv,
                params.cell_type_metadata_csv,
                params.leiden_res,
                outdir,
                params.project,
                params.docker_enabled,
                params.r_module,
                params.python_module,
            )
            ch_versions = ch_versions.unique().mix(ANNOTATE_SCINA_SCANPY_WF.out.versions)
        } else if ( params.clustering_mode == "seurat" ) {
            ANNOTATE_SCINA_SEURAT_WF (
                CLUSTER_SEURAT_WF.out.scrnaseq_object,
                params.cell_type_markers_csv,
                params.cell_type_metadata_csv,
                params.leiden_res,
                outdir,
                params.project,
                params.docker_enabled,
                params.r_module,
            )
            ch_versions = ch_versions.unique().mix(ANNOTATE_SCINA_SEURAT_WF.out.versions)
        } else {
            error "'clustering_mode' parameter needs to be: 'scanpy' or 'seurat'"
        }

        // Automatic annotation following ambient RNA correction + clustering
        if ( params.run_soupx == "Yes" && cellranger_count_info ) {
            if ( params.clustering_mode == "scanpy" ) {
                ANNOTATE_SCINA_SCANPY_SOUPX_WF (
                    CLUSTER_SCANPY_SOUPX_WF.out.scrnaseq_object_raw,
                    CLUSTER_SCANPY_SOUPX_WF.out.cell_metadata_csv,
                    CLUSTER_SCANPY_SOUPX_WF.out.scrnaseq_object,
                    params.cell_type_markers_csv,
                    params.cell_type_metadata_csv,
                    params.leiden_res,
                    outdir_soupx,
                    out_name_soupx,
                    params.docker_enabled,
                    params.r_module,
                    params.python_module,
                )
            } else if ( params.clustering_mode == "seurat" ) {
                ANNOTATE_SCINA_SEURAT_SOUPX_WF (
                    CLUSTER_SEURAT_SOUPX_WF.out.scrnaseq_object,
                    params.cell_type_markers_csv,
                    params.cell_type_metadata_csv,
                    params.leiden_res,
                    outdir_soupx,
                    out_name_soupx,
                    params.docker_enabled,
                    params.r_module,
                )
            } else {
                error "'clustering_mode' parameter needs to be: 'scanpy' or 'seurat'"
            }
        }


        // Write final versions file
        ch_versions
        .collectFile (
            name: "versions.txt", 
            storeDir: "${outdir}", 
            keepHeader: false
        )
}
