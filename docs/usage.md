# vsearchpipeline: Usage

## Introduction

This bioinformatic pipeline uses VSEARCH to generate a count table (samples+ASVs) from 16S sequencing reads. The input is a samplesheet with sample names and file paths to the fastq files, and a sheet with primer sequences if primer trimming is necessary. The pipeline uses DADA2 for taxonomic assignment (as opposed to VSEARCH sintax) using the SILVA v.138.1 database. The resulting count table, taxonomic table and phylogenetic tree resulting from the pipeline are stored in a phyloseq object.

As an experimental feature, there is a rarefaction/pruning process and a process that assembles 'nice' assembled taxonomic names ('Roseburia spp.' if species level is unknown, Lachnospiraceae spp., if genus level is unknown, etcetera). It is possible to skip these processes. Either way, both the unrarefied and the rarefied phyloseq objects are stored in the output directory.

There is a separate [paragraph](#hpc-settings) on use of the pipeline on computing clusters in this document. I had some struggles getting Nextflow to work on an HPC with a lot of strange settings, and hope this will help others in similar situations.

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Full samplesheet

The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 3 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 6 samples, where `TREATMENT_REP3` has been sequenced twice.

```console
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,AEG588A4_S4_L003_R2_001.fastq.gz
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,AEG588A5_S5_L003_R2_001.fastq.gz
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,AEG588A6_S6_L003_R2_001.fastq.gz
TREATMENT_REP3,AEG588A6_S6_L004_R1_001.fastq.gz,AEG588A6_S6_L004_R2_001.fastq.gz
```

| Column           | Description                                                                                        |
| ---------------- | -------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `nucl_acid_conc` | _(Optional)_ DNA/RNA concentration of the sample (e.g. measured by Qubit). Required for frequency-based decontam. Leave empty for negative controls and samples without a measurement. |

An [example samplesheet](../testdata/Samplesheet.csv) has been provided with the pipeline.

## Decontam

If your experiment includes **negative (extraction/PCR) controls**, you can run [`decontam`](https://benjjneb.github.io/decontam/) to identify and remove contaminant ASVs before downstream analyses.

### Enabling decontam

```bash
--run_decontam
```

### Requirements

| Requirement | Details |
|---|---|
| **Negative control samples** | Sample names must contain `NEGCON` (e.g. `NEGCON_1`, `Sample_NEGCON`). These are used for prevalence-based decontam. |
| **DNA concentration** _(optional)_ | Add a `nucl_acid_conc` column to your samplesheet (e.g. Qubit ng/µL). This enables frequency-based decontam in addition to prevalence-based. The original samplesheet file (`--input`) is passed directly to the decontam process and the column is read from there. Rows with missing or zero concentration are excluded from frequency-based decontam only — they still contribute to prevalence-based decontam. If the column is absent entirely, only prevalence-based decontam runs. |

Both methods are run if both are available; each is run individually if only one data source is present. If neither negative controls nor concentration data are found, the process saves the unmodified phyloseq and writes a warning to the report.

> [!NOTE]
> Prevalence-based decontam (using `NEGCON` samples) runs independently of concentration data. If you have negative controls but no DNA concentrations — or only partial concentration data — decontam will still run using the prevalence method.

### Decontam methods

1. **Frequency-based** (`isContaminant(method = "frequency", conc = "nucl_acid_conc")`): Identifies ASVs whose abundance is inversely correlated with DNA concentration — a hallmark of contamination.
2. **Prevalence-based** (`isContaminant(method = "prevalence", neg = "Ctrl")`): Identifies ASVs that are more prevalent in negative controls than in true samples.

The union of both methods is removed.

### Control sample removal

After contaminant ASVs are removed, all samples whose names match `NEGCON` or start with `SB` are removed from the phyloseq object. This also happens automatically in the **rarefaction step** as a safeguard, even when `--run_decontam` is not set.

### Output

Results are saved to `<outdir>/phyloseq/decontam/`:

| File | Description |
|---|---|
| `phyloseq_decontam.RDS` | Cleaned phyloseq object (contaminants and control samples removed) |
| `decontam_contaminants.csv` | List of ASVs identified as contaminants |
| `decontam_prev_plot.pdf` | Prevalence plot (negative controls vs. true samples) |
| `decontam_report.txt` | Summary report of the decontam run |

### Effect on rarefaction

When `--run_decontam` is enabled, the **rarefaction step uses the decontam output** (`phyloseq_decontam.RDS`) as its input instead of the complete phyloseq object. The complete (pre-decontam) phyloseq object is always saved separately.

## Primersheet input

If you want the primers to be trimmed, you need to specify a file path to a csv with the forward and reverse primer sequences that were used.

```bash
--primers '[path to primersheet file]'
```

This csv should look as follows:

```console
forward_primer, reverse_primer
CCTACGGGAGGCAGCAG,TACNVGGGTATCTAAKCC
```

| Column           | Description                                                                                        |
| ---------------- | -------------------------------------------------------------------------------------------------- |
| `forward_primer` | Sequence of the forward primer. All nucleotides and degenerate nucleotides are allowed characters. |
| `reverse_primer` | Sequence of the reverse primer. All nucleotides and degenerate nucleotides are allowed characters. |

An [example primersheet](../testdata/Primers.csv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run vsearchpipeline --input ./Samplesheet.csv --outdir ./results --primers ./Primers.csv -profile singularity
```

This will launch the pipeline with the `singularity` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run barbarahelena/vsearchpipeline -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
primers: './primers.csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch) or by writing them manually.

### HPC settings

I had some problems getting Nextflow pipelines to work on an HPC. This is a short list of issues I encountered with solutions. These solutions are all implemented, although for slurm use and the job name replacement you have to use `-profile snellius`.

#### Use singularity

Docker usually does not work on HPCs for security reasons; the alternative is Singularity (or Apptainer). Run the pipeline with `-profile singularity`. If you use `-profile snellius`, singularity is also enabled.

#### Retry failed jobs

Because slurm jobs might occasionally fail, some processes have the `error_retry` label. With this label, jobs are retried 5 times before the pipeline fails. This is already implemented in the pipeline.

#### Use slurm

If you're using a scheduler on a computing cluster, you should define that in the process scope configuration file. The queue size and rate limit of job submissions are also defined. The configs below are used in the `-profile snellius` but with a specified partition name. Make sure you replace the partition name when using this profile on other computing clusters.

```
process{
    executor = 'slurm'
    clusterOptions = '-p partitionname'
    queueSize = 50
	  submitRateLimit = '10 sec'
}
```

#### Job names with special characters

On some HPCs, certain characters such as parentheses and columns are not allowed in the slurm job names. The configs below are used in the `-profile snellius`.

```
executor{
    jobName = {
        n = "${task.name}" -> n.replaceAll( '[():\\s]', '_')
    }
}
```

### Updating the pipeline

To make sure that you're running the latest version of the pipeline, make sure that you regularly clone the git repository of the pipeline:

```bash
git clone git@github.com:barbarahelena/vsearchpipeline.git
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since. You can specify the version number when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!NOTE]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations when it runs. Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.
- `snellius`
  - A configuration specific for the ['Snellius HPC'](). It uses the settings described in the paragraph about ['HPC settings'](#hpc-settings). In short: singularity, slurm configuration and job name configuration. You can use this profile for other HPCs but you should replace the name of the partition where jobs are submitted to.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, edit the labels in `conf/base.config` directly, or supply a custom config file with `-c`.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. Containers and software are sourced from [BioContainers](https://biocontainers.pro/) or [Bioconda](https://bioconda.github.io/).

To use a different container, override `container` for the relevant process label or name in a custom config file passed with `-c`.

### Custom Tool Arguments

To pass additional arguments to a particular tool, set `ext.args` for the relevant process in a custom config file. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for details.

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
