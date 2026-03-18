# vsearchpipeline: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0 - Mar 18 2026

### `Added`

- New `VSEARCH_MAPPINGRATE` process that calculates per-sample mapping rates from existing `usearch_global` outputs without re-running vsearch
- `VSEARCH_FASTQFILTER` now captures stderr (`*.filter_stats.txt`) to record post-filter read counts per sample, used as the denominator for mapping rate
- `VSEARCH_USEARCHGLOBAL` now captures stderr (`mapping_stats.txt`) to record overall unique-sequence mapping rate
- Mapping rate output: `mapping_rate_summary.tsv` (per-sample) and `mapping_rate_overall.txt` (global)
- SILVA database version now configurable via pipeline parameter (`silva_db`)
- PiCrust2 stratified option added
- Decontam module added with negative control detection and nucleic acid concentration input
- Picrust2 module with extArgs support
- `PIPELINE_INITIALISATION` and `PIPELINE_COMPLETION` subworkflows following modern nf-core pattern
- Samplesheet reading and primer validation moved into `PIPELINE_INITIALISATION` (utils subworkflow)
- Named `VSEARCHFLOW` wrapper workflow in `main.nf` following modern nf-core pattern
- MultiQC skipped automatically with a warning when input exceeds 6000 FastQC files

### `Fixed`

- Fixed rarefaction level auto-calculation: level is now capped at `max_counts` so it can never exceed the deepest sample and wipe all samples; `min_counts` and `max_counts` are computed once and reused throughout
- Fixed rarefaction user-level guard: fallback to auto-calculated level now only triggers when `user_rarelevel > max_counts` (i.e. no samples would survive), instead of incorrectly triggering whenever the level exceeded the shallowest sample
- Fixed decontam bugs including column name for nucleic acid concentration (`Nucl_Acid_Conc`)
- Fixed overall mapping rate not being found due to incorrect search string (`"Matching unique query sequences"` vs `"Matching query sequences"`)
- Fixed version collection for R-based processes
- Fixed test file paths and testdata layout
- Fixed schema for metamap input
- Fixed small bugs in vsearch scripts
- Fixed broken link in docs
- Removed nf-core template boilerplate, logos and references
- Removed unused `check_max` function
- Cleaned repo structure
- Fixed `MissingMethodException` crash in MultiQC caused by FastQC zip `LinkedList` being passed as closure arguments
- Fixed MultiQC channel construction: zip files were incorrectly passed as `--config` and `--logo` CLI arguments
- Removed unused `ch_multiqc_config`, `ch_multiqc_logo`, and `ch_multiqc_custom_config` channels
- Removed deprecated script-level `def` variable declarations; all logic now inside `workflow {}` blocks
- Removed `NfcoreTemplate`, `WorkflowMain`, and `WorkflowVsearchpipeline` Groovy helper class usage
- Inlined `INPUT_CHECK` and `PRIMERS_CHECK` subworkflows; subworkflow files are now unused
- Updated `assets/multiqc_config.yml` section IDs and links to match `barbarahelena/vsearchpipeline`
- Updated `assets/methods_description_template.yml` with correct pipeline references and full tool citation list
- Replaced metaboflow ASCII art in pipeline logo with VSEARCH ASCII art

### `Dependencies`

- Updated FastQC and MultiQC nf-core modules
- Updated DADA2 to v1.38.0
- Updated SILVA database to 138.2
- Updated vsearch scripts
- Updated mafft version

### `Deprecated`

- `subworkflows/local/input_check.nf` — replaced by `PIPELINE_INITIALISATION`
- `subworkflows/local/primers_check.nf` — replaced by `PIPELINE_INITIALISATION`
- `lib/NfcoreTemplate.groovy`, `lib/WorkflowMain.groovy`, `lib/WorkflowVsearchpipeline.groovy` — replaced by utils subworkflows

---

## v0.7.1-beta - Feb 18 2024

### `Fixed`

- Updated FastQC and MultiQC nf-core modules to latest versions

---

## v0.7-beta - Jan 23 2024

### `Added`

- New module to sort and remove singletons
- `skip_tree` parameter to make phylogenetic tree construction optional
- `maxdiffpct` parameter for vsearch merge pairs
- New HPC profile with wider CPU/memory/time limits
- `error_retry` label for failed jobs
- Zenodo DOI badge

### `Fixed`

- Solved bugs in phyloseq modules
- Ensured phyloseq outputs are published to outdir
- Vertical line in rarefaction plot now placed at the actual rarefaction level
- Fixed memory and resource labels
- Corrected typos in documentation

### `Dependencies`

- Replaced IQ-TREE with FastTree for phylogenetic tree construction
- Updated vsearch and DADA2 containers
- Updated nf-core modules
- Returned to SILVA 138.1 (latest stable)
- Moved container and conda settings to `modules.config`

### `Deprecated`

- IQ-TREE option removed

---

## v0.6-beta - Nov 6 2023

### `Added`

- Check for presence of primers before trimming (primer check is now optional)
- Testdata FASTQ files without primers (`*_noprimers.fastq.gz`)
- `skip_primers` parameter

### `Fixed`

- Resolved issues in metrics and rarefaction modules
- Fixed rarefaction level fallback when all sample counts are very low
- Restructured modules for reusability with postfix (complete / rarefied)

### `Dependencies`

- Updated phyloseq container

---

## v0.5-beta - Nov 6 2023

Initial release of vsearchpipeline, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Full VSEARCH-based amplicon pipeline: merge pairs, quality filter, dereplication, chimera removal, clustering (UNOISE3), taxonomy assignment (DADA2 / SILVA), phyloseq object construction
- Optional phylogenetic tree (FastTree / IQ-TREE) and MSA (MAFFT)
- Optional PiCrust2 functional prediction
- Rarefaction and diversity/composition metrics modules
- Decontam module for negative-control-based contaminant removal
- MultiQC and FastQC quality control
- nf-core-compatible samplesheet and primer sheet input
- Test profile with small testdata
- HPC (Snellius) profile

### `Dependencies`

- VSEARCH
- DADA2
- phyloseq
- SILVA 138.1
- FastTree / MAFFT
- PiCrust2
- FastQC / MultiQC
