#!/usr/bin/env nextflow



// Enable DSL2 for more powerful functionality
nextflow.enable.dsl = 2



// Workflow

/*
 * prepare input for pipeline
 */

workflow PIPELINE_INIT_WF {

    take:
        sample_sheet

    main:

        // process input sample sheet
        ch_input = Channel.empty()
        file_ext = file( sample_sheet, checkIfExists: true ).Extension

        // extract sample IDs and fastq directories from sample sheet
        if ( file_ext = "txt" ) {
        
            Channel
                .fromPath( sample_sheet )
                .splitCsv( header: true, sep: "\t" )
                .map { line ->
                    tuple( line.sample_id, file( line.fastq_dir, type: "dir", checkIfExists: true ) )
                }
		.unique()
                .set { ch_txt }
            ch_input = ch_input.mix( ch_txt )

        } else if ( file_ext = "csv" ) {

            Channel
                .fromPath( sample_sheet )
                .splitCsv( header: true )
                .map { line ->
                    tuple( line.sample_id, file( line.fastq_dir, type: "dir", checkIfExists: true ) )
                }
                .set { ch_csv }
            ch_input = ch_input.mix( ch_csv )

        } else {

            println "Error: Input sample sheet should be a tab-delimited '.txt' file. Set sample sheet file path in the nextflow.config file."

        }

    emit:
        input = ch_input
}