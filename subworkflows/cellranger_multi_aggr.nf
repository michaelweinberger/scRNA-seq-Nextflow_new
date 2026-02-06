#!/usr/bin/env nextflow



// Enable DSL2 for more powerful functionality
nextflow.enable.dsl = 2



// Processes

/*
 * generate cell ranger multi config .csv file
 */
process CELLRANGER_MULTI_CONFIG_PR {
    cpus 1
    debug false
    publishDir (
        "${outdir}/cellranger", 
        pattern: "", 
        mode: "copy", 
        saveAs: { filename -> "${filename}" }
    )

    input:
    path  ( cellranger_ref_dir )                // Cell Ranger reference genome absolute path
    tuple val( sample_id ), path( fastq_dir )
    path  ( sample_sheet )          		// sample IDs & metadata
    val   ( outdir )
    
    output:
    tuple val( sample_id ), path( "cellranger_multi_config.csv" ), emit: cellranger_multi_config   // emitted for cellranger multi

    script:
    """
    cellranger_ref_dir_1="\$(cat ${cellranger_ref_dir})"

    2_cellranger_multi_input.sh -i "${sample_id}" -s "${sample_sheet}" -r "\${cellranger_ref_dir_1}"
    """
}



/*
 * run cell ranger multi
 */
process CELLRANGER_MULTI_PR {
    debug       false
    tag         "$sample_id"
    publishDir (
        "${outdir}/cellranger", 
        pattern: "versions.txt", 
        mode: "copy", 
        saveAs: { filename -> "${filename}" }
    )
    publishDir (
        "${outdir}/cellranger", 
        pattern: "cellranger_count_info.tsv", 
        mode: "copy", 
        saveAs: { filename -> "${filename}" }
    )
    publishDir (
        "${outdir}/cellranger", 
        pattern: "cellranger_count_dir/outs", 
        mode: "copy", 
       saveAs: { fn -> "${sample_id.replace(',', '_')}_count/outs" }
    )

    container { ( params.docker_enabled == "true" || params.docker_enabled == true ) ? "litd/docker-cellranger:v9.0.1" : "" }

    input:
    tuple val( sample_id ), path( cellranger_multi_config_csv )
    val   ( outdir )
    val   ( docker_enabled )
    val   ( cellranger_module )
    
    output:
    path ( "cellranger_count_info.tsv" ), emit: cellranger_count_info 
    path ( "cellranger_count_dir/outs" ), emit: cellranger_count_dir  // emitted for publishing
    path ( "versions.txt"              ), emit: versions

    script:
    def sample_id_clean = sample_id.replace(',', '_')

    """
    if [ "$docker_enabled" == "false" ] ; then
        module load "$cellranger_module"
    fi

    # run cell ranger multi
    cellranger multi \
	--id="cellranger_count_dir" \
	--csv="$cellranger_multi_config_csv"

    printf "sample_id\tcellranger_dir\n" > cellranger_count_info.tsv
    cellranger_dir="\$(echo \${PWD}/cellranger_count_dir)"
    printf "${sample_id_clean}\t\${cellranger_dir}\n" >> cellranger_count_info.tsv

    echo "${task.process}:" > versions.txt
        echo cellranger: "\$(cellranger --version 2>&1 | awk '{print \$(NF)}' )" | sed -e \$'s/^/\t/' >> versions.txt
    """
}



/*
 * generate cell ranger aggr input .csv file
 */
process CELLRANGER_AGGR_INPUT_PR {
    cpus 1
    debug false
    publishDir (
        "${outdir}/cellranger", 
        pattern: "", 
        mode: "copy", 
        saveAs: { filename -> "${filename}" }
    )

    input:
    path  ( sample_sheet )          // sample IDs & metadata
    path  ( cellranger_count_info ) // sample IDs & full filepaths cellranger count output dirs
    val   ( outdir )
    
    output:
    path  ( "cellranger_aggr_input.csv" ), emit: cellranger_aggr_input // emitted for cellranger aggr

    script:
    """
    2_cellranger_multi_aggr_input.sh -s "${sample_sheet}" -i "${cellranger_count_info}"
    """
}



/*
 * generate barcode metadata file
 */
process CELLRANGER_METADATA_PR {
    cpus 1
    debug false
    publishDir (
        "${outdir}/cellranger", 
        pattern: "", 
        mode: "copy", 
        saveAs: { filename -> "${filename}" }
    )

    input:
    path ( cellranger_aggr_input_csv )
    val  ( outdir )
    
    output:
    path ( "cellranger_aggr_cell_metadata.tsv" ), emit: metadata   // emitted for scanpy, seurat and loom merge

    script:
    """
    3_cellranger_metadata.sh -i "${cellranger_aggr_input_csv}"
    """
}



/*
 * run cell ranger aggr
 */
process CELLRANGER_AGGR_PR {
    debug false
    publishDir (
        "${outdir}/cellranger", 
        pattern: "versions.txt", 
        mode: "copy", 
        saveAs: { filename -> "${filename}" }
    )
    publishDir (
        "${outdir}/cellranger", 
        pattern: "${out_name}/outs", 
        mode: "copy", 
        saveAs: { fn -> "${out_name}/outs" }
    )

    container { ( params.docker_enabled == "true" || params.docker_enabled == true ) ? "litd/docker-cellranger:v9.0.1" : "" }

    input:
    path ( cellranger_aggr_input_csv )
    val  ( out_name )
    val  ( outdir )
    val  ( docker_enabled )
    val  ( cellranger_module )
    
    output:
    path ( "${out_name}/outs/count/filtered_feature_bc_matrix" ), emit: cellranger_aggr_bc_matrix // emitted for scanpy,seurat
    path ( "${out_name}/outs"                                  ), emit: cellranger_aggr_outdir    // emitted for publishing
    path ( "versions.txt"                                      ), emit: versions

    script:
    """
    if [ "$docker_enabled" == "false" ] ; then
        module load "$cellranger_module"
    fi

    # Run cellranger aggr
    cellranger aggr --id="$out_name" \
	--csv="$cellranger_aggr_input_csv"

    echo "${task.process}:" > versions.txt
        echo cellranger: "\$(cellranger --version 2>&1 | awk '{print \$(NF)}' )" | sed -e \$'s/^/\t/' >> versions.txt
    """
}



// Workflow

workflow CELLRANGER_MULTI_AGGR_WF {

    take:
        cellranger_ref_dir
        sample_sheet
        input                     // tuple sample ID + fastq dir
        cellranger_aggr_out_name
        outdir
        docker_enabled
        cellranger_module

    main:
        ch_versions = Channel.empty()

        // generate cellranger multi config .csv file
        CELLRANGER_MULTI_CONFIG_PR (
            cellranger_ref_dir,
            input,
            sample_sheet,
            outdir,
        )

        // run cellranger multi
        CELLRANGER_MULTI_PR (
            CELLRANGER_MULTI_CONFIG_PR.out.cellranger_multi_config,
            outdir,
            docker_enabled,
            cellranger_module,
        )
        ch_versions = ch_versions.mix(CELLRANGER_MULTI_PR.out.versions)

        // collect sample IDs and full paths of cellranger output dirs
        CELLRANGER_MULTI_PR.out.cellranger_count_info
        .collectFile (
            name: "cellranger_count_info.tsv", 
            storeDir: "${outdir}/cellranger", 
            keepHeader: true
        )
        .set { cellranger_count_info }

        // generate tuple of sample IDs and full paths of cellranger count output dirs for velocyto
        cellranger_count_info
        .splitCsv( header: true , sep: "\t" )
        .map { line ->
            tuple( line.sample_id, line.cellranger_dir )
        }
        .set { cellranger_count_out }

        // generate cellranger aggr input .csv file
        CELLRANGER_AGGR_INPUT_PR (
            sample_sheet, 
            cellranger_count_info,
            outdir,
        )

        // generate cell barcode level metadata file
        CELLRANGER_METADATA_PR (
            CELLRANGER_AGGR_INPUT_PR.out.cellranger_aggr_input,
            outdir,
        )

        // run cellranger aggr
        CELLRANGER_AGGR_PR ( 
            CELLRANGER_AGGR_INPUT_PR.out.cellranger_aggr_input,
            cellranger_aggr_out_name,
            outdir,
            docker_enabled,
            cellranger_module,
        )
        ch_versions = ch_versions.mix(CELLRANGER_AGGR_PR.out.versions)

    emit:
        versions                   = ch_versions
        cellranger_count_info      = cellranger_count_info
        cellranger_count_out       = cellranger_count_out
        cellranger_metadata        = CELLRANGER_METADATA_PR.out.metadata
        cellranger_aggr_bc_matrix  = CELLRANGER_AGGR_PR.out.cellranger_aggr_bc_matrix
}
