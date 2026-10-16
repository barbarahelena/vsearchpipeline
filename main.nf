#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    vsearchpipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/barbarahelena/vsearchpipeline
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { VSEARCHPIPELINE         } from './workflows/vsearchpipeline'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_vsearchpipeline_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_vsearchpipeline_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMETERS
    Typed: command-line values are converted to these types at launch, and a wrong type,
    a missing file or an undeclared (misspelled) parameter stops the run before it starts.
    Infrastructure params read by config files (outdir, multiqc_*, config_profile_*, ...)
    are defined in nextflow.config. Keep defaults in sync with nextflow_schema.json.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params {
    // Input (required for a run; optional here so --help works without them)
    input: Path?                        = null
    primers: Path?                      = null   // required unless --skip_primers

    // Merge read pairs (vsearch --fastq_mergepairs)
    merge_allowmergestagger: Boolean    = false
    merge_maxdiffs: Integer             = 30
    merge_minlen: Integer               = 0      // 0 = no minimum
    merge_maxlen: Integer               = 0      // 0 = no maximum
    merge_maxdiffpct: Integer           = 100

    // Quality filter (vsearch --fastq_filter)
    filter_maxee: Float                 = 1.0
    filter_maxns: Integer               = 0
    filter_minlen: Integer              = 0
    filter_maxlen: Integer              = 0

    // Dereplication
    derep_strand: String                = 'plus'
    derep_all_strand: String            = 'plus'
    derep_all_fastawidth: Integer       = 0
    derep_all_minunique: Integer        = 2

    // Denoising (vsearch --cluster_unoise)
    cluster_minsize: Integer            = 10     // vsearch default is 8; 10 was more accurate on mock communities
    cluster_alpha: Float                = 2.0

    // Sort and remove singletons
    sort_fastawidth: Integer            = 0
    sort_minsize: Integer               = 2

    // Chimera removal and count table
    uchime_label: String                = 'ASV_'
    usearch_id: Float                   = 0.97

    // Taxonomy (DADA2 with SILVA)
    silva_asv_db: Path?                 = null   // downloaded when not given
    silva_species_db: Path?             = null
    save_silva_db: Boolean              = false
    dada2_minboot: Integer              = 80
    dada2_allowmultiple: Integer        = 3
    dada2_tryrevcompl: Boolean          = true

    // phyloseq
    rarelevel: Integer                  = 0      // 0 = choose automatically
    run_decontam: Boolean               = false

    // Skip steps
    skip_primers: Boolean               = true
    skip_tree: Boolean                  = false
    skip_rarefaction: Boolean           = false
    skip_fixtaxonomy: Boolean           = false
    skip_metrics: Boolean               = false
    skip_picrust: Boolean               = true

    // Generic
    help_full: Boolean                  = false
    show_hidden: Boolean                = false
    version: Boolean                    = false
    validate_params: Boolean            = true
    plaintext_email: Boolean            = false
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        args,
        params.outdir,
        params.input,
        params.primers,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    VSEARCHPIPELINE (
        PIPELINE_INITIALISATION.out.reads,
        PIPELINE_INITIALISATION.out.primers
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs
    )

    publish:
    fastqc            = VSEARCHPIPELINE.out.fastqc
    seqtk             = VSEARCHPIPELINE.out.seqtk
    vsearch           = VSEARCHPIPELINE.out.vsearch
    mafft             = VSEARCHPIPELINE.out.mafft
    fasttree          = VSEARCHPIPELINE.out.fasttree
    silva_db          = VSEARCHPIPELINE.out.silva_db
    dada2             = VSEARCHPIPELINE.out.dada2
    picrust2          = VSEARCHPIPELINE.out.picrust2
    phyloseq_complete = VSEARCHPIPELINE.out.phyloseq_complete
    phyloseq_decontam = VSEARCHPIPELINE.out.phyloseq_decontam
    phyloseq_rarefied = VSEARCHPIPELINE.out.phyloseq_rarefied
    multiqc           = VSEARCHPIPELINE.out.multiqc
    pipeline_info     = VSEARCHPIPELINE.out.pipeline_info
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW OUTPUTS: where each published channel goes, relative to --outdir
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

output {
    fastqc            { path 'fastqc' }
    seqtk             { path 'seqtk' }
    vsearch           { path 'vsearch' }
    mafft             { path 'mafft' }
    fasttree          { path 'fasttree' }
    silva_db          { path 'silva_db' }
    dada2             { path 'dada2' }
    picrust2          { path 'picrust2' }
    phyloseq_complete { path 'phyloseq/complete' }
    phyloseq_decontam { path 'phyloseq/decontam' }
    phyloseq_rarefied { path 'phyloseq/rarefied' }
    multiqc           { path 'multiqc' }
    pipeline_info     { path 'pipeline_info' }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
