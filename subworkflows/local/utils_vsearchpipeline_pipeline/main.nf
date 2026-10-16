//
// Subworkflow with functionality specific to the vsearchpipeline pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    primers           //  string: Path to primers CSV, or null
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

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
    def before_text = """
-\033[2m----------------------------------------------------\033[0m-
\033[0;34m  ${workflow.manifest.name} ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    def after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/', '')}" }.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/barbarahelena/vsearchpipeline/blob/main/CITATIONS.md
"""
    def command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        null
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
    def ch_reads = channel
        .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .map { meta, fastq_1, fastq_2 ->
            def new_meta = meta + [ single_end: !fastq_2 ]
            def reads    = fastq_2 ? [ fastq_1, fastq_2 ] : [ fastq_1 ]
            [ new_meta, reads ]
        }

    //
    // Create and validate primer channel (if primers are not skipped)
    //
    def ch_primers
    if (!params.skip_primers) {
        if (!primers) {
            error "[vsearchpipeline] params.primers is not set. Please provide a primers file with --primers or set --skip_primers true."
        }
        ch_primers = channel.fromPath(primers, checkIfExists: true)
            .map { csv -> parsePrimersheet(csv) }
            .first()
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

    main:

    def summary_params = paramsSummaryMap()

    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs)
        }
        completionSummary(monochrome_logs)
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
// Git commit of the pipeline code that is running.
// workflow.commitId is set when run from a repository (nextflow run barbarahelena/vsearchpipeline -r ...);
// otherwise ask git, and return 'unknown' rather than crash when git is missing (e.g. some HPC login nodes).
//
def gitCommitHash() {
    if (workflow.commitId) {
        return workflow.commitId
    }
    try {
        def proc = ['git', '-C', projectDir.toString(), 'rev-parse', 'HEAD'].execute()
        def out  = proc.text.trim()
        return proc.waitFor() == 0 && out ? out : 'unknown'
    } catch (Exception _e) {
        return 'unknown'
    }
}

//
// Generate methods description for MultiQC
//
def methodsDescriptionText(mqc_methods_yaml) {
    def meta          = [:]
    meta.workflow     = workflow.toMap()
    meta.manifest_map = workflow.manifest.toMap()

    if (meta.manifest_map.doi) {
        meta["doi_text"] = meta.manifest_map.doi.tokenize(",")
            .collect { doi_ref -> doi_ref.replace('https://doi.org/', '').replace(' ', '') }
            .collect { clean -> "(doi: <a href='https://doi.org/${clean}'>${clean}</a>)" }
            .join(", ")
    } else {
        meta["doi_text"] = ""
    }
    meta["nodoi_text"]        = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"
    meta["tool_citations"]    = ""
    meta["tool_bibliography"] = ""

    def engine = new groovy.text.SimpleTemplateEngine()
    return engine.createTemplate(mqc_methods_yaml.text).make(meta).toString()
}
