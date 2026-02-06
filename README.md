# scRNA-seq Nextflow Pipeline

A comprehensive Nextflow pipeline for the analysis of single-cell RNA sequencing (scRNA-seq) data. This workflow automates genome reference generation, alignment, quality control, clustering, ambient RNA correction, and cell type annotation.

---

## Pipeline Overview

This pipeline supports both **Scanpy** and **Seurat** based analysis paths and includes the following steps:

1.  **Reference Generation**
    * Automated download of genome FASTA and GTF files from Ensembl.
    * Generation of Cell Ranger genome indices.
2.  **Alignment & Counting**
    * Supports two modes via **[Cellranger](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-ct)**:
        * **Standard Count**: Standard gene expression quantification (via `--input_count`).
        * **Multi**: Multiplexing/Feature Barcode analysis (via `--input_multi`).
    * Performs sample aggregation via `cellranger aggr`.
3.  **Quality Control**
    * Transcriptome/Proteome BLAST checks.
4.  **Clustering & Integration**
    * User-selectable analysis backend: **[Scanpy](https://scanpy.readthedocs.io/en/stable/)<sup>1</sup>** or **[Seurat](https://satijalab.org/seurat/)<sup>2</sup>**.
    * Includes QC filtering.
    * Includes doublet detection via [DoubletDetection](https://github.com/JonathanShor/DoubletDetection?tab=readme-ov-file)<sup>3</sup> (as part of Scanpy analysis) or [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder)<sup>4</sup> (as part of Seurat analysis)
    * Dimensionality reduction (PCA) and Visualization (UMAP).
    * Sample integration via **[Harmony](https://github.com/immunogenomics/harmony)<sup>5</sup>**.
5.  **Ambient RNA Correction (Optional)**
    * Integration with **[SoupX] (https://github.com/constantAmateur/SoupX)<sup>6</sup>** to remove ambient RNA contamination \
    (enabled via `--run_soupx "Yes"`).
6.  **Cell Type Annotation**
    * Automated cell type assignment using **[SCINA] (https://github.com/jcao89757/SCINA)<sup>7</sup>** based on provided marker gene lists.

---

## Requirements

* **[Nextflow](https://www.nextflow.io/docs/latest/index.html)<sup>8</sup>** (>= 25.04.7)
* **Docker** (if running with `docker_enabled: true`) or appropriate software modules loaded in your environment (R, Python, Cell Ranger).

---

## Usage

#### 1. Clone the Repository

```bash
git clone [https://github.com/michaelweinberger/scRNA-seq-Nextflow.git](https://github.com/michaelweinberger/scRNA-seq-Nextflow.git)
cd scRNA-seq-Nextflow
```

#### 2. Adjust the `nextflow.config` file

#### 3.  Start the pipeline via `nextflow run` like
```
nextflow run main.nf \
-profile HPC_docker,mouse \
-resume
```

In this example, the "HPC_no_docker" profile directs Nextflow to run the pipeline in a high-performance computing environment using the slurm scheduler, and the "mouse" profile sets some genome parameters used to generate cellranger genome index files. \
When setting the `-resume` flag, the pipeline will resume from a previous run.

<br/>



## Usage scenarios

### Starting from fastq files

If the `input_count` or `input_multi` parameters are set in the nextflow.config file, the workflow will start from fastq files. See below under "Parameters" for more information about `input` specifications.

The workflow first generates Cell Ranger genome index files. Input fastq files are then aligned to the reference genome via the `cellranger count` or `cellranger multi` commands and the output from individual samples is combined via `cellranger aggr`. This is followed by doublet removal, clustering, data integration and cell type annotation.

<br/>

### Starting from Cell Ranger outputs

To start the workflow from previously generated Cell Ranger output files instead of from fastq files, leave the `input` parameters in the nextflow.config file empty. Instead, specify the `cellranger_out_dir` and `metadata` parameters. These specify a file path to a directory containing cell Ranger outputs, and a tab delimited file with cell barcodes and metadata, respectively. Please see below under "Parameters" for more details.

<br/>



## Parameters

All parameters can be set on the command line with `--parameter_name` or in the `params` scope of the nextflow.config file located in the main pipeline directory.

- `project`   The name of the analysis project, defaults to "nf_scRNAseq_analysis".

- `outdir`    The name of the directory to save pipeline outputs to, defaults to "/out" in the main pipeline directory.\
            Within this directory, outputs from individual parts of the pipeline are written to different subdirectories:\
            - "genomes" for Cellranger genome index files\
            - "cellranger" for Cellranger mapping outputs\
            - "scanpy" for Scanpy clustering outputs (without SoupX)\
                - "scina" for SCINA cell type annotation outputs\
            - "seurat" for Seurat clustering outputs (without SoupX)\
                - "scina" for SCINA cell type annotation outputs\
            - (If enabled) "soupx" for SoupX outputs:\
                "soupx", "scanpy" or "seurat"

<br/>

### Parameters specific to starting from fastq files (set only one):

- `input_count`     The file path to an input ".txt" (tab delimited) file containing sample information to run Cell Ranger Count.\
            - The first column of the sample sheet must be named "sample_id" and contain sample-specific identifiers that also are prefixes in the corresponding fastq file names.\
            For example: Put "sample_x" as sample ID if your fastq files are named "sample_x_S2_L001_I1_001.fastq.gz", "sample_x_S2_L001_R1_001.fastq.gz", "sample_x_S2_L001_R2_001.fastq.gz" etc.\
            - The second column must be named "fastq_dir" and contain file paths to directories with fastq files to be analysed. Fastq files of multiple samples may be located within the same directory. \
	    - Further metadata columns can be added: e.g. "sample_name", "tissue", "condition" etc.

- `input_multi`     The file path to an input ".txt" (tab delimited) file containing sample information to run Cell Ranger Multi.\
            - The first column must be named "parent_id" and contain the sample-specific prefixes of the fastq files, \
              IMPORTANT: use identical parent IDs to name gene expression and multiplex capture (CMO) fastq files originating from the same sample \
              -> add suffix "GEX" to gene expression parent ID/fastq file name, add suffix "CMO" to multiplex capture parent ID/fastq file name. \
              Include only GEX parent IDs (not CMO parent IDs) in the sample sheet. \
              If samples have been re-sequenced, supply all prefixes as a comma-separated list. \
            - The second column must be named "fastq_dir" and contain file paths to directories containing GEX and CMO fastq files specified by "parent_id". \
              If samples have been re-sequenced and are stored in multiple directories, supply all directories as a comma-separated list. \
            - The third column must be named "sample_id" and contain the names of the de-multiplexed samples present in the fastq files, \
              put one sample ID per line and duplicate "parent_id" and "fastq_dir" column entries as necessary. \
            - The fourth column must be named "cmo_id" and contain the CMO tags corresponding to de-multiplexed sample IDs. \
	    - Further metadata columns can be added: e.g. "sample_name", "tissue", "condition" etc.
<br/>

### Parameters specific to starting from Cellranger output files:

- `cellranger_out_dir`  The file path to a directory containing cellranger output "barcodes.tsv.gz", "features.tsv.gz" and "matrix.mtx.gz" files, typically a "/outs/count/filtered_feature_bc_matrix" directory. \
ALTERNATIVELY, this directory can contain '.h5' files with Cell Ranger mapped data.

- `metadata`      The file path to a tab delimited file containing \
                - a column named "barcode" of cell barcodes, \
                - a column named "sample_id" of sample identifiers, \
                - optional metadata columns \
                When  '.h5' , the metadata file does not need to contain a "barcode" column. \
                Doublet detection is performed using the "sample_id" column.

- `mapping_mode`  The mode that Cell Ranger was run in to generate the mapped data ("cell ranger count" or "cell ranger multi")

- `cellranger_info_tsv` (Optional) If you would like to run SoupX to correct ambient RNA contamination, supply a file path to a \
                tab-delimited text file with columns "sample_id", containing sample IDs for Cell Ranger Count or \
                non de-multiplexed parent IDs for Cell Ranger Multi, \
                and "cellranger_dir", containing paths to Cell Ranger output directories (parent directory of "outs/")
<br/>

### Parameters related to scRNA-seq clustering:

- `clustering_mode`   Clustering analysis tool to use. Values: "scanpy" or "seurat"
- `drop_ribo_prot`    Should ribosomal protein genes be excluded from highly variable genes before clustering? Values: "Yes" or "No"
- `genes_drop_csv`    Optional: Path to a csv file containing genes to be dropped from highly variable genes before clustering
- `min_genes`         Minimum number of genes expressed for a cell to be kept in the dataset, defaults to 200
- `max_genes`         Maximum number of genes expressed for a cell to be kept in the dataset, defaults to 2500
- `max_perc_mt`       Maximum percentage of mitochondrial reads for a cell to be kept in the dataset, defaults to 5
- `min_cells`         Minimum number of cells for a gene to be expressed in to be kept in the dataset, defaults to 3
- `n_pcs`             Number of principal components to be computed, defaults to 30
- `harmony_var`       Name of the metadata column to use for Harmony data integration, defaults to ""
- `leiden_res`        Resolution of cell clustering, defaults to 0.4
- `run_soupx`         Should ambient RNA contamination be corrected? Values: "Yes" or "No"

<br/>




## Profiles

Multiple parameters can be bundled into profiles. These are defined in the `profiles` scope in the nextflow.config file and can be invoked on the command line via the `-profile` flag. Additional executor or genome profiles can be added in the nextflow.config file.

<br/>

### Pre-defined executor profiles

#### HPC_no_docker
Use for pipeline execution on a high performance cluster without using Docker.

#### HPC_docker
Use for pipeline execution on a high performance cluster that allows the use of Docker.

#### standard
Use for local pipeline execution.

<br/>

### Pre-defined genome profiles

#### human

#### mouse

#### zebrafish

<br/>



## References
1.	Wolf, F.A., Angerer, P., and Theis, F.J. (2018). SCANPY: large-scale single-cell gene expression data analysis. Genome Biol 19, 15. [10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0)
2.	Hao, Y., Hao, S., Andersen-Nissen, E., Mauck, W.M., 3rd, Zheng, S., Butler, A., Lee, M.J., Wilk, A.J., Darby, C., Zager, M., et al. (2021). Integrated analysis of multimodal single-cell data. Cell 184, 3573-3587 e3529. [10.1016/j.cell.2021.04.048](https://doi.org/10.1016/j.cell.2021.04.048)
3.	Gayoso, Adam, Shor, Jonathan, Carr, Ambrose J., Sharma, Roshan, Pe'er, Dana (2020, December 18). DoubletDetection (Version v3.0). Zenodo. [10.5281/zenodo.2678041](https://doi.org/10.5281/zenodo.2678041)
4.	McGinnis, C.S., Murrow, L.M., and Gartner, Z.J. (2019). DoubletFinder: Doublet Detection in Single-Cell RNA Sequencing Data Using Artificial Nearest Neighbors. Cell Syst 8, 329-337 e324. [10.1016/j.cels.2019.03.003](https://doi.org/10.1016/j.cels.2019.03.003)
5.	Korsunsky, I., Millard, N., Fan, J., Slowikowski, K., Zhang, F., Wei, K., Baglaenko, Y., Brenner, M., Loh, P.R., and Raychaudhuri, S. (2019). Fast, sensitive and accurate integration of single-cell data with Harmony. Nat Methods 16, 1289-1296. [10.1038/s41592-019-0619-0](https://doi.org/10.1038/s41592-019-0619-0)
6.	Young, M.D., Behjati, S. (2020). SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data. GigaScience, Volume 9, Issue 12, December 2020. [10.1093/gigascience/giaa151](https://doi.org/10.1093/gigascience/giaa151)
7.	Zhang Z., Luo D., et al. (2018). SCINA: A Semi-Supervised Subtyping Algorithm of Single Cells and Bulk Samples. Genes, Volume 10, Issue 7. [10.3390/genes10070531](https://doi.org/10.3390/genes10070531)
8.  P. Di Tommaso, et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology 35, 316–319. [10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820)