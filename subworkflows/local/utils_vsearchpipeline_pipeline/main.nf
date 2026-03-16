//
// Subworkflow with functionality specific to the vsearchpipeline pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFVALIDATION_PLUGIN } from '../../nf-core/utils_nfvalidation_plugin'
include { paramsSummaryMap          } from 'plugin/nf-validation'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { dashedLine                } from '../../nf-core/utils_nfcore_pipeline'
include { nfCoreLogo                } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { workflowCitation          } from '../../nf-core/utils_nfcore_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved

    main:

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    pre_help_text  = nfCoreLogo(monochrome_logs)
    post_help_text = '\n' + workflowCitation() + '\n' + dashedLine(monochrome_logs)
    def String workflow_command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"
    UTILS_NFVALIDATION_PLUGIN (
        help,
        workflow_command,
        pre_help_text,
        post_help_text,
        validate_params,
        "nextflow_schema.json"
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Create channel from input samplesheet
    //
    channel.fromSamplesheet('input')
        .map { meta, fastq_1, fastq_2 ->
            def new_meta = meta + [ single_end: !fastq_2 ]
            def reads    = fastq_2 ? [ fastq_1, fastq_2 ] : [ fastq_1 ]
            [ new_meta, reads ]
        }
        .set { ch_reads }

    //
    // Create and validate primer channel (if primers are not skipped)
    //
    if (!params.skip_primers) {
        if (!params.primers) {
            error "[vsearchpipeline] params.primers is not set. Please provide a primers file with --primers or set --skip_primers true."
        }
        channel.fromPath(params.primers, checkIfExists: true)
            .map { csv -> parsePrimersheet(csv) }
            .first()
            .set { ch_primers }
    } else {
        ch_primers = channel.empty()
    }

    emit:
    reads   = ch_reads
    primers = ch_primers
}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications

    main:

    summary_params = paramsSummaryMap(workflow)

    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs)
        }
        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://github.com/barbarahelena/vsearchpipeline"
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

//
// Parse and validate a one-row primer CSV file
//
def parsePrimersheet(csv) {
    def rows = csv.splitCsv(header: true, strip: true)
    if (rows.size() != 1) {
        error "[vsearchpipeline] Primersheet must contain exactly one row of primers, found ${rows.size()} rows."
    }
    def row = rows[0]
    def fwd = row['forward_primer'] ?: row['fwd_primer']
    def rev = row['reverse_primer'] ?: row['rev_primer']
    if (!fwd || !rev) {
        error "[vsearchpipeline] Primersheet must have columns 'forward_primer' and 'reverse_primer' (or 'fwd_primer'/'rev_primer')."
    }
    def iupac = '[ACGTURYSWKMBDHVNacgturyswkmbdhvn]+'
    if (!fwd.matches(iupac)) {
        error "[vsearchpipeline] Forward primer '${fwd}' is not a valid IUPAC DNA sequence."
    }
    if (!rev.matches(iupac)) {
        error "[vsearchpipeline] Reverse primer '${rev}' is not a valid IUPAC DNA sequence."
    }
    return [ forward: fwd.trim(), reverse: rev.trim() ]
}

//
// Generate methods description for MultiQC
//
def methodsDescriptionText(mqc_methods_yaml) {
    def meta          = [:]
    meta.workflow     = workflow.toMap()
    meta.manifest_map = workflow.manifest.toMap()

    if (meta.manifest_map.doi) {
        def temp_doi_ref = ""
        String[] manifest_doi = meta.manifest_map.doi.tokenize(",")
        for (String doi_ref : manifest_doi) {
            def clean = doi_ref.replace('https://doi.org/', '').replace(' ', '')
            temp_doi_ref += "(doi: <a href='https://doi.org/${clean}'>${clean}</a>), "
        }
        meta["doi_text"] = temp_doi_ref[0..-3]
    } else {
        meta["doi_text"] = ""
    }
    meta["nodoi_text"]        = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"
    meta["tool_citations"]    = ""
    meta["tool_bibliography"] = ""

    def engine = new groovy.text.SimpleTemplateEngine()
    return engine.createTemplate(mqc_methods_yaml.text).make(meta).toString()
}
