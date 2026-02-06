#!/usr/bin/env nextflow



// Enable DSL2 for more powerful functionality
nextflow.enable.dsl = 2



// Processes

/*
 * Download genome files
 */
process GENOME_FILES_PR {
    debug false
    publishDir "${outdir}/genomes", pattern: "", mode: "copy", overwrite: true, saveAs: { filename -> "${filename}" }

    container { ( params.docker_enabled == "true" || params.docker_enabled == true ) ? "michaelweinberger/ubuntu-22.04:v1" : "" }
    
    input:
    val  ( species )
    val  ( genome )
    val  ( genome_ucsc )
    val  ( species_latin )
    val  ( ensembl_version )
    val  ( outdir )
    val  ( docker_enabled )
    
    output:
    path ( "${genome}.fa"                     ), emit: genome_fasta, optional: true
    path ( "${genome}.${ensembl_version}.gtf" ), emit: genome_gtf, optional: true
    path ( "${genome}_rmsk.txt"               ), emit: genome_repeat_mask, optional: true

    script:
    """
    # capitalise first letter
    species_latin_1="${species_latin}"
    species_latin_2="\${species_latin_1^}"

    # download fasta and gtf files
    rsync -avzP "rsync://ftp.ebi.ac.uk/ensemblorg/pub/release-${ensembl_version}/fasta/${species_latin}/dna/\${species_latin_2}.${genome}.dna.toplevel.fa.gz" .
    mv "\${species_latin_2}.${genome}.dna.toplevel.fa.gz" "${genome}.fa.gz"
    gunzip "${genome}.fa.gz"

    rsync -avzP "rsync://ftp.ebi.ac.uk/ensemblorg/pub/release-${ensembl_version}/gtf/${species_latin}/\${species_latin_2}.${genome}.${ensembl_version}.gtf.gz" .
    mv "\${species_latin_2}.${genome}.${ensembl_version}.gtf.gz" "${genome}.${ensembl_version}.gtf.gz"
    gunzip "${genome}.${ensembl_version}.gtf.gz"

    # download file containing positions of repetitive elements in genome
    #wget -L "http://hgdownload.soe.ucsc.edu/goldenPath/${genome_ucsc}/database/rmsk.txt.gz" .
    #gunzip "rmsk.txt.gz"
    #mv "rmsk.txt" "${genome}_rmsk.txt"
    """
}



/*
 * Prepare cellranger genome index
 */
process CELLRANGER_REF_PR {
    debug false
    publishDir "${outdir}/genomes", pattern: "", mode: "copy", overwrite: true, saveAs: { filename -> "${filename}" }

    container { ( params.docker_enabled == "true" || params.docker_enabled == true ) ? "litd/docker-cellranger:v9.0.1" : "" }
    
    input:
    val  ( species )
    val  ( genome )
    path ( genome_fasta )
    path ( genome_gtf )
    val  ( outdir )
    val  ( docker_enabled )
    val  ( cellranger_module )
    
    output:
    path ( "refdata-cellranger-${genome}"                 ), emit: cellranger_index
    path ( "cellranger_ref_dir.txt"                       ), emit: cellranger_ref_dir
    path ( "refdata-cellranger-${genome}/genes/genes.gtf" ), emit: cellranger_gtf_file
    path ( "versions.txt"                                 ), emit: versions

    script:
    """
    if [ "$docker_enabled" = "false" ] ; then
        module load "$cellranger_module"
    fi

    # prepare cellranger reference data
    if [ "$species" = "human" ] ; then
	cellranger_ref="refdata-gex-GRCh38-2024-A"
    elif [ "$species" = "mouse" ] ; then
	cellranger_ref="refdata-gex-GRCm39-2024-A"
    fi

    if [ "$species" = "human" ] || [ "$species" = "mouse" ] ; then
        wget "https://cf.10xgenomics.com/supp/cell-exp/\${cellranger_ref}.tar.gz"
        tar -zxvf "\${cellranger_ref}.tar.gz"
        rm "\${cellranger_ref}.tar.gz"
        mv "\${cellranger_ref}" "refdata-cellranger-${genome}"
    else 
        cellranger mkgtf \
            "${genome_gtf}" \
            "${genome}.filtered.gtf" \
            --attribute=gene_biotype:protein_coding \
    	    --attribute=gene_biotype:lncRNA \
    	    --attribute=gene_biotype:antisense \
    	    --attribute=gene_biotype:IG_LV_gene \
    	    --attribute=gene_biotype:IG_V_gene \
    	    --attribute=gene_biotype:IG_V_pseudogene \
    	    --attribute=gene_biotype:IG_D_gene \
    	    --attribute=gene_biotype:IG_J_gene \
    	    --attribute=gene_biotype:IG_J_pseudogene \
    	    --attribute=gene_biotype:IG_C_gene \
    	    --attribute=gene_biotype:IG_C_pseudogene \
    	    --attribute=gene_biotype:TR_V_gene \
    	    --attribute=gene_biotype:TR_V_pseudogene \
    	    --attribute=gene_biotype:TR_D_gene \
    	    --attribute=gene_biotype:TR_J_gene \
    	    --attribute=gene_biotype:TR_J_pseudogene \
    	    --attribute=gene_biotype:TR_C_gene

        cellranger mkref \
            --genome="refdata-cellranger-${genome}" \
            --fasta="${genome_fasta}" \
            --genes="${genome}.filtered.gtf" \
            --output-dir="\${PWD}/refdata-cellranger-${genome}"
    fi

    gunzip "refdata-cellranger-${genome}/genes/genes.gtf.gz"

    # write absolute path to Cell Ranger reference genome
    echo "\${PWD}/refdata-cellranger-${genome}" > cellranger_ref_dir.txt

    echo "${task.process}:" > versions.txt
        echo cellranger: "\$(cellranger --version 2>&1 | awk '{print \$(NF)}' )" | sed -e \$'s/^/\t/' >> versions.txt
    """
}



// Workflow

workflow CELLRANGER_REF_WF {

    take:
        species
        genome
        genome_ucsc
        species_latin
        ensembl_version
        outdir
        docker_enabled
        cellranger_module
	
    main:
        ch_versions = Channel.empty()

        GENOME_FILES_PR (
            species,
            genome,
            genome_ucsc,
            species_latin,
            ensembl_version,
            outdir,
            docker_enabled,
        )

        CELLRANGER_REF_PR (
            species,
            genome,
            GENOME_FILES_PR.out.genome_fasta,
            GENOME_FILES_PR.out.genome_gtf,
            outdir,
            docker_enabled,
            cellranger_module,
        )
        ch_versions         = ch_versions.mix(CELLRANGER_REF_PR.out.versions)
        cellranger_index    = CELLRANGER_REF_PR.out.cellranger_index
        cellranger_ref_dir  = CELLRANGER_REF_PR.out.cellranger_ref_dir
        cellranger_gtf_file = CELLRANGER_REF_PR.out.cellranger_gtf_file

    emit:
        versions            = ch_versions
        cellranger_index    = cellranger_index
        cellranger_ref_dir  = cellranger_ref_dir
        cellranger_gtf_file = cellranger_gtf_file
        //genome_repeat_mask  = GENOME_FILES_PR.out.genome_repeat_mask
}
