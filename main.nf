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
        params.help_full.toString().toBoolean(),   // CLI flags can arrive as the string "true"
        params.show_hidden.toString().toBoolean()
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
