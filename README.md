# VSEARCH nextflow pipeline

[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.10076630-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.10076630)

[![Nextflow](https://img.shields.io/badge/Nextflow%20-%20%3E%3D25.10.2%20-%20green?logo=Nextflow)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**vsearchpipeline** is a bioinformatics pipeline that uses [VSEARCH](https://github.com/torognes/vsearch) to infer ASVs and make a count table from 16S sequencing reads. The input is a samplesheet with sample names and file paths to the fastq files, and a sheet with primer sequences if primer trimming is necessary. The pipeline uses DADA2 for taxonomic assignment using the SILVA reference database. The resulting count table, taxonomic table and phylogenetic tree resulting from the pipeline are stored in a phyloseq object.
![vsearch](https://github.com/user-attachments/assets/be411868-e12e-4a54-9b6d-60dad859814c)

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Trim primers ([`Seqtk`](https://github.com/lh3/seqtk))
3. Infer ASVs and make count table ([`VSEARCH`](https://github.com/torognes/vsearch))
4. Multiple sequence alignment ([`MAFFT`]()) to make phylogenetic tree ([`Fasttree`](https://www.microbesonline.org/fasttree/))
5. Taxonomy assignment ([`DADA2`](https://benjjneb.github.io/dada2/)) using SILVA 138.1 database for DADA2
6. Phyloseq object with count table, taxonomic table and phylogenetic tree ([`Phyloseq`](https://joey711.github.io/phyloseq/))
7. MultiQC report ([`MultiQC`](http://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
> to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
> with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a pair of fastq files (paired-end).

Then, prepare a sheet with the forward and reverse primers. This sheet should look as follows:

`primers.csv`

```console
forward_primer, reverse_primer
CCTACGGGAGGCAGCAG,TACNVGGGTATCTAAKCC
```

If there are no primers to be trimmed, simply add the `--skip_primers` flag to the nextflow run command.

Now, you can run the pipeline using:

```bash
nextflow run barbarahelena/vsearchpipeline \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --primers primers.csv \
   --outdir <OUTDIR>
```

For more details and further functionality, please refer to the [usage documentation](https://github.com/barbarahelena/vsearchpipeline/blob/master/docs/usage.md) and the [parameter documentation](https://github.com/barbarahelena/vsearchpipeline/blob/master/docs/parameters.md).

## Pipeline output

All output of the different parts of the pipeline are stored in subdirectories of the output directory. These directories are named after the tools that were used ('vsearch', 'dada2', etc.). In the phyloseq folder, you can find the end result of the pipeline, which is the phyloseq object. Other important outputs are the multiqc report in the multiqc folder and the execution html report in the pipeline_info folder.

For more details on the pipeline output, please refer to the [output documentation](https://github.com/barbarahelena/vsearchpipeline/blob/master/docs/output.md).

## Citations

If you use this VSEARCH workflow for your analysis, please cite it using the following doi: [10.5281/zenodo.10076629](https://doi.org/10.5281/zenodo.10076629). An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
