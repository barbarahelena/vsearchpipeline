# vsearchpipeline: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.0 - Oct 19 2026

Strict syntax and static typing throughout. **Requires Nextflow 26.10.0 or later.**

### Breaking changes

- Nextflow 26.10.0 or later is required (was 25.10.3)
- Misspelled or unknown parameters now stop the run instead of being ignored; values must have the declared type
- `--max_cpus`, `--max_memory`, `--max_time` removed: set `process.resourceLimits` in a config file instead (see the usage docs)
- `--hook_url` (Slack/Teams notifications) removed
- Custom configs that set `publishDir` no longer have an effect: results are published through workflow outputs
- Process names lost the `VSEARCHFLOW:` prefix (`VSEARCHPIPELINE:...`): update `withName:` selectors in custom configs, and expect `-resume` from a v1 run to rerun all tasks
- Results on macOS without containers change, because primer trimming now actually runs there

### `Changed`

- Strict syntax everywhere (`nextflow lint` reports no errors); typed processes for all local modules; typed workflows for `main.nf` and `VSEARCHPIPELINE`
- Samples are `Sample` records (`id`, `single_end`, `fastq`) instead of `meta` maps; per-sample steps take only the fields they need and return a record with the sample id and their files
- All analysis parameters declared once, with types, in a `params` block in `main.nf`. Wrong types, misspelled parameters and missing input files now stop the run at launch
- Parameter validation uses `nf-schema` 3.0.0 instead of the deprecated `nf-validation`; JSON schemas converted to draft 2020-12
- nf-core utility subworkflows replaced by current versions, tracked in `modules.json`; nf-core FastQC (0.13.0) and MultiQC (1.35) modules updated
- VSEARCH updated to 2.32.0 (bug fixes; thread count now respects container/Slurm CPU limits; quality scores above 41 accepted). Results on the test data are unchanged
- Results are published through workflow outputs (`output` block in `main.nf`) instead of `publishDir`; the results folder layout is unchanged
- `--max_cpus`, `--max_memory`, `--max_time` replaced by `process.resourceLimits` (they were not enforced)
- Tool versions collected from the `versions` topic channel; MultiQC now lists all tools (it previously showed only FastQC)
- `git_commit_hash.txt` written by Nextflow itself; no longer needs git inside a container
- `task.ext.args` is now passed to every tool
- The `arm` profile also builds conda environments for `osx-64`, so `-profile mamba,arm` works on Apple Silicon
- nf-test suite (`tests/default.nf.test`): the test profile with checksums of the scientific results, and a stub run with all optional steps
- CI: `nextflow lint` job, then nf-test on Nextflow 26.10.0 and latest stable

### `Fixed`

- Primer trimming was silently skipped on macOS without containers (`zcat` there only reads `.Z` files)
- Boolean parameters set to `false` on the command line (e.g. `--skip_picrust false`) were treated as `true`
- ASV numbering could differ between runs of the same data (per-sample files were concatenated in task finish order)
- `--help` showed `cluster_minsize` default 8; the default used is 10
- `*.merged.fastq.gz` files were not gzip-compressed
- `ext.args` set in a config was ignored
- Fasttree, seqtk and vsearch versions were recorded incorrectly; PICRUSt2 and dereplication stub blocks were broken

### `Removed`

- `--hook_url` (Slack/Teams notifications), dropped from the nf-core utilities

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
- GitHub Actions CI workflow (`ci.yml`) with matrix jobs for `test` and `test_noprimers` profiles
- GitHub Actions benchmark workflow (`benchmark.yml`) triggered on release or manually; validates ASV recovery from 22-species mock community
- GitHub Actions `check_updates.yml` workflow: checks Bioconda/conda-forge weekly for new tool versions and sends email notification via Gmail
- `assets/check_updates.py` — Anaconda API-based version checker with bioconda → conda-forge fallback
- `assets/tools.txt` — list of pipeline tool packages to monitor for updates
- `assets/versions_seen.json` — tracks last known version per tool; auto-committed by the workflow
- Stub SILVA databases for CI (`assets/testdata/small/SILVA_asv_db_test.fa.gz`, `SILVA_species_db_test.fa.gz`) filtered to taxa present in testdata
- Stub SILVA databases for benchmark (`assets/testdata/benchmark/SILVA_benchmark_asv_db.fa.gz`, `SILVA_benchmark_species_db.fa.gz`) filtered to the 19 mock community genera
- `conf/test_noprimers.config` — CI test profile without primer trimming
- `conf/test_benchmark.config` — benchmark profile using mock community data with known 22-species composition
- `assets/testdata/benchmark/mock_sequences_V4.fasta` — 22 reference ASV sequences for benchmark validation
- `assets/testdata/benchmark/validate_benchmark.py` — validates recovered ASV count (21–23) and per-reference identity ≥ 0.97
- `docs/check_updates.md` — setup instructions for the version monitoring workflow

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
