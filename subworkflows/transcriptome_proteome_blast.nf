#!/usr/bin/env nextflow



// Enable DSL2 for more powerful functionality
nextflow.enable.dsl = 2



// Processes

/*
 * Make transcriptome fasta files
 */
process MAKE_TRANSCRIPTOME_FASTA_PR {
    debug false
    publishDir "${outdir}/genomes", pattern: "", mode: "copy", overwrite: true, saveAs: { filename -> "${filename}" }

    container { ( params.docker_enabled == "true" || params.docker_enabled == true ) ? "michaelweinberger/ubuntu-22.04:v1" : "" }
    
    input:
    val  ( species )
    val  ( biomart_dataset )
    val  ( biomart_genome )
    val  ( chrom_number )
    val  ( human_biomart_dataset )
    val  ( human_biomart_genome )
    val  ( outdir )
    val  ( docker_enabled )
    val  ( r_module )
    
    output:
    path ( "${biomart_genome}_cds.fa"       ), emit: transcriptome_fasta
    path ( "${human_biomart_genome}_cds.fa" ), emit: human_transcriptome_fasta
    path ( "versions.txt"                   ), emit: versions

    script:
    """
    if [ "$docker_enabled" == "false" ] ; then
        module load "$r_module"
    fi

    4_1_make_transcriptome_fasta.R --args species="$species" \
        out_dir="\$PWD" \
        biomart_dataset="$biomart_dataset" \
        biomart_genome="$biomart_genome" \
        chrom_number="$chrom_number" \
        human_biomart_dataset="$human_biomart_dataset" \
        human_biomart_genome="$human_biomart_genome"

    echo "${task.process}:" > versions.txt
        R --version | head -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(BiocManager); print(paste('BiocManager', packageVersion('BiocManager')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(dplyr); print(paste('dplyr', packageVersion('dplyr')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(tidyr); print(paste('tidyr', packageVersion('tidyr')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(GenomicFeatures); print(paste('GenomicFeatures', packageVersion('GenomicFeatures')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(rtracklayer); print(paste('rtracklayer', packageVersion('rtracklayer')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(BSgenome); print(paste('BSgenome', packageVersion('BSgenome')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(gprofiler2); print(paste('gprofiler2', packageVersion('gprofiler2')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(biomaRt); print(paste('biomaRt', packageVersion('biomaRt')))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(eval(as.name($biomart_genome))); print(paste($biomart_genome, packageVersion($biomart_genome)))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
        R -e "library(eval(as.name($human_biomart_genome))); print(paste($human_biomart_genome, packageVersion($human_biomart_genome)))" | grep '[1]' | tail -n 1 | sed -e \$'s/^/\t/' >> versions.txt
    """
}



/*
 * Make proteome fasta files
 */
process MAKE_PROTEOME_FASTA_PR {
    debug false
    publishDir "${outdir}/genomes", pattern: "", mode: "copy", overwrite: true, saveAs: { filename -> "${filename}" }

    container { ( params.docker_enabled == "true" || params.docker_enabled == true ) ? "michaelweinberger/ubuntu-22.04:v1" : "" }
    
    input:
    val  ( species )
    val  ( uniprot_id )
    val  ( ncbi_id )
    val  ( outdir )
    val  ( docker_enabled )
    
    output:
    path ( "${species}_proteome.fa" ), emit: proteome_fasta
    path ( "human_proteome.fa"      ), emit: human_proteome_fasta

    script:
    """
    # download proteome fasta file (not "_DNA.fasta" file)
    wget "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/${uniprot_id}/${uniprot_id}_${ncbi_id}.fasta.gz"
    gunzip "${uniprot_id}_${ncbi_id}.fasta.gz"
    mv "${uniprot_id}_${ncbi_id}.fasta" "${species}_proteome_tmp.fa"

    # download human proteome fasta file
    wget "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz"
    gunzip "UP000005640_9606.fasta.gz"
    mv "UP000005640_9606.fasta" "human_proteome_tmp.fa"

    # adjust fasta file headers to only contain gene names (in GN field)
    awk '/^>/ { for(i=1; i<=NF; i++) { if(\$i ~ /^GN=/) { sub(/^GN=/, ">", \$i); print \$i; next; } } } /^[^>]/ { print \$0 }' "${species}_proteome_tmp.fa" > "${species}_proteome.fa"
    grep -c "^>" "${species}_proteome.fa"

    awk '/^>/ { for(i=1; i<=NF; i++) { if(\$i ~ /^GN=/) { sub(/^GN=/, ">", \$i); print \$i; next; } } } /^[^>]/ { print \$0 }' "human_proteome_tmp.fa" > "human_proteome.fa"
    grep -c "^>" "human_proteome.fa"
    """
}



/*
 * Blast-compare non-human to human genes, using transcriptome or proteome fasta files
 */
process BLAST_FASTA_PR {
    debug false
    publishDir "${outdir}/genomes/transcriptome_maps", pattern: "", mode: "copy", overwrite: true, saveAs: { filename -> "${filename}" }

    container { ( params.docker_enabled == "true" || params.docker_enabled == true ) ? "michaelweinberger/ubuntu-22.04:v1" : "" }
    
    input:
//    path ( transcriptome_fasta )
//    path ( human_transcriptome_fasta )
    path ( proteome_fasta )
    path ( human_proteome_fasta )
    val  ( species )
    val  ( outdir )
    val  ( docker_enabled )
    val  ( blast_module )
    
    output:
    path ( "${species}_to_human.txt" ), emit: transcriptome_blast_txt

    script:
    """
    if [ "$docker_enabled" == "false" ] ; then
        module load "$blast_module"
    fi

    # script can be run with transcriptome fasta files and --t1 + --t2 set to "nucl", or with proteome fasta files and --t1 + --t2 set to "prot"
    # map_genes.sh downloaded from https://github.com/atarashansky/SAMap/blob/main/
    4_2_map_genes.sh --tr1 "${proteome_fasta}" \
                     --t1 "prot" \
                     --n1 "${species}" \
                     --tr2 "${human_proteome_fasta}" \
                     --t2 "prot" \
                     --n2 "human" \
                     --threads "${task.cpus}"

    # add header
    echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" > standard_blast_fields.tsv
    cat standard_blast_fields.tsv "${species}_to_human.txt" > t && mv t "${species}_to_human.txt"
    """
}



// Workflow

workflow TRANSCRIPTOME_PROTEOME_BLAST_WF {

    take:
        species
        biomart_dataset
        biomart_genome
        chrom_number
        human_biomart_dataset
        human_biomart_genome
        uniprot_id
        ncbi_id
        outdir
        docker_enabled
        r_module
        blast_module

    main:
        ch_versions = Channel.empty()

//        MAKE_TRANSCRIPTOME_FASTA_PR (
//            species,
//            biomart_dataset,
//            biomart_genome,
//            chrom_number,
//            human_biomart_dataset,
//            human_biomart_genome,
//            outdir,
//            docker_enabled,
//            r_module,
//        )
//        ch_versions = ch_versions.mix(MAKE_TRANSCRIPTOME_FASTA_PR.out.versions)

        MAKE_PROTEOME_FASTA_PR (
            species,
            uniprot_id,
            ncbi_id,
            outdir,
            docker_enabled,
        )

        BLAST_FASTA_PR (
//            MAKE_TRANSCRIPTOME_FASTA_PR.out.transcriptome_fasta,
//            MAKE_TRANSCRIPTOME_FASTA_PR.out.human_transcriptome_fasta,
            MAKE_PROTEOME_FASTA_PR.out.proteome_fasta,
            MAKE_PROTEOME_FASTA_PR.out.human_proteome_fasta,
            species,
            outdir,
            docker_enabled,
            blast_module,
        )

    emit:
        versions                = ch_versions
        transcriptome_blast_txt = BLAST_FASTA_PR.out.transcriptome_blast_txt
}
