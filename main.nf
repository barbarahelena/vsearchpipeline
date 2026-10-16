#!/usr/bin/env nextflow

nextflow.enable.types = true
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

    // Output and notifications (outdir and monochrome_logs are also set in nextflow.config, which config files read; config wins)
    outdir: String?                     = null
    email: String?                      = null
    email_on_fail: String?              = null
    monochrome_logs: Boolean            = false
    multiqc_methods_description: String? = null

    // Generic
    help: String?                       = null   // --help, or --help <parameter> for one parameter
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
    def init = PIPELINE_INITIALISATION (
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
    // The options VSEARCHPIPELINE uses, as a typed record (see VsearchOptions in workflows/vsearchpipeline.nf)
    def opts = record(
        input                      : params.input,
        outdir                     : params.outdir,
        skip_primers               : params.skip_primers,
        merge_allowmergestagger    : params.merge_allowmergestagger,
        merge_maxdiffs             : params.merge_maxdiffs,
        merge_minlen               : params.merge_minlen,
        merge_maxlen               : params.merge_maxlen,
        merge_maxdiffpct           : params.merge_maxdiffpct,
        filter_maxee               : params.filter_maxee,
        filter_minlen              : params.filter_minlen,
        filter_maxlen              : params.filter_maxlen,
        filter_maxns               : params.filter_maxns,
        derep_strand               : params.derep_strand,
        derep_all_strand           : params.derep_all_strand,
        derep_all_fastawidth       : params.derep_all_fastawidth,
        derep_all_minunique        : params.derep_all_minunique,
        cluster_minsize            : params.cluster_minsize,
        cluster_alpha              : params.cluster_alpha,
        sort_fastawidth            : params.sort_fastawidth,
        sort_minsize               : params.sort_minsize,
        uchime_label               : params.uchime_label,
        usearch_id                 : params.usearch_id,
        skip_tree                  : params.skip_tree,
        silva_asv_db               : params.silva_asv_db,
        silva_species_db           : params.silva_species_db,
        save_silva_db              : params.save_silva_db,
        dada2_minboot              : params.dada2_minboot,
        dada2_allowmultiple        : params.dada2_allowmultiple,
        dada2_tryrevcompl          : params.dada2_tryrevcompl,
        skip_picrust               : params.skip_picrust,
        run_decontam               : params.run_decontam,
        skip_fixtaxonomy           : params.skip_fixtaxonomy,
        skip_metrics               : params.skip_metrics,
        skip_rarefaction           : params.skip_rarefaction,
        rarelevel                  : params.rarelevel,
        multiqc_methods_description: params.multiqc_methods_description
    )

    def pipeline = VSEARCHPIPELINE (
        init.reads,
        init.primers,
        opts
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
    fastqc            = pipeline.fastqc
    seqtk             = pipeline.seqtk
    vsearch           = pipeline.vsearch
    mafft             = pipeline.mafft
    fasttree          = pipeline.fasttree
    silva_db          = pipeline.silva_db
    dada2             = pipeline.dada2
    picrust2          = pipeline.picrust2
    phyloseq_complete = pipeline.phyloseq_complete
    phyloseq_decontam = pipeline.phyloseq_decontam
    phyloseq_rarefied = pipeline.phyloseq_rarefied
    multiqc           = pipeline.multiqc
    pipeline_info     = pipeline.pipeline_info
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW OUTPUTS: where each published channel goes, relative to --outdir
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

output {
    fastqc: Channel<Path>            { path 'fastqc' }
    seqtk: Channel<Path>             { path 'seqtk' }
    vsearch: Channel<Path>           { path 'vsearch' }
    mafft: Channel<Path>             { path 'mafft' }
    fasttree: Channel<Path>          { path 'fasttree' }
    silva_db: Channel<Path>          { path 'silva_db' }
    dada2: Channel<Path>             { path 'dada2' }
    picrust2: Channel<Path>          { path 'picrust2' }
    phyloseq_complete: Channel<Path> { path 'phyloseq/complete' }
    phyloseq_decontam: Channel<Path> { path 'phyloseq/decontam' }
    phyloseq_rarefied: Channel<Path> { path 'phyloseq/rarefied' }
    multiqc: Channel<Path>           { path 'multiqc' }
    pipeline_info: Channel<Path>     { path 'pipeline_info' }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
