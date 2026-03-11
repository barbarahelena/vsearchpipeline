/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowVsearchpipeline.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Local modules
//
include { SEQTK_TRIMFQ }                                            from '../modules/local/seqtk/trimfq'
include { VSEARCH_FASTQMERGEPAIRS }                                 from '../modules/local/vsearch/fastqmergepairs'
include { VSEARCH_FASTQFILTER }                                     from '../modules/local/vsearch/fastqfilter'
include { VSEARCH_DEREPFULLLENGTH }                                 from '../modules/local/vsearch/derepfulllength'
include { VSEARCH_DEREPFULLLENGTHALL }                              from '../modules/local/vsearch/derepfulllengthall'
include { VSEARCH_CLUSTERUNOISE }                                   from '../modules/local/vsearch/clusterunoise'
include { VSEARCH_SORT_REMOVE_SINGLETONS }                          from '../modules/local/vsearch/sort_remove_singletons'
include { VSEARCH_UCHIMEDENOVO }                                    from '../modules/local/vsearch/uchimedenovo'
include { VSEARCH_USEARCHGLOBAL }                                   from '../modules/local/vsearch/usearchglobal'
include { MAFFT }                                                   from '../modules/local/mafft'
include { FASTTREE }                                                from '../modules/local/fasttree'
include { SILVADATABASES }                                          from '../modules/local/silvadatabases'
include { DADA2_ASSIGNTAXONOMY }                                    from '../modules/local/dada2/assigntaxonomy'
include { PICRUST2 }                                                 from '../modules/local/picrust2'
include { PHYLOSEQ_MAKEOBJECT as PHYLOSEQ_COMPLETE_MAKEOBJECT }     from '../modules/local/phyloseq/makeobject'
include { PHYLOSEQ_FIXTAXONOMY as PHYLOSEQ_COMPLETE_FIXTAX }        from '../modules/local/phyloseq/fixtaxonomy'
include { PHYLOSEQ_METRICS as PHYLOSEQ_COMPLETE_METRICS }           from '../modules/local/phyloseq/metrics'
include { PHYLOSEQ_DECONTAM }                                        from '../modules/local/phyloseq/decontam'
include { PHYLOSEQ_RAREFACTION as PHYLOSEQ_RAREFIED }               from '../modules/local/phyloseq/rarefaction'
include { PHYLOSEQ_METRICS as PHYLOSEQ_RAREFIED_METRICS }           from '../modules/local/phyloseq/metrics'
include { PHYLOSEQ_FIXTAXONOMY as PHYLOSEQ_RAREFIED_FIXTAX }        from '../modules/local/phyloseq/fixtaxonomy'


//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { PRIMERS_CHECK } from '../subworkflows/local/primers_check'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow VSEARCHPIPELINE {
    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    if(!params.skip_primers){
        PRIMERS_CHECK (
            file(params.primers)
        )
        ch_primers = PRIMERS_CHECK.out.primers.first()
        ch_versions = ch_versions.mix(PRIMERS_CHECK.out.versions)
    }
    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    //
    // MODULE:  Seqtk trim primers in fastq files
    //
    if (!params.skip_primers) {
        SEQTK_TRIMFQ (
            INPUT_CHECK.out.reads, 
            ch_primers
        ).reads.set {ch_trimmed_reads}
    } else {
        ch_trimmed_reads = INPUT_CHECK.out.reads
    }    

    //
    // MODULE: VSEARCH merge fastq pairs
    //
    VSEARCH_FASTQMERGEPAIRS (
        ch_trimmed_reads,
        params.merge_allowmergestagger,
        params.merge_maxdiffs,
        params.merge_minlen,
        params.merge_maxlen,
        params.merge_maxdiffpct
    )

    //
    // MODULE: VSEARCH filter fastq files
    //
    VSEARCH_FASTQFILTER (
        VSEARCH_FASTQMERGEPAIRS.out.reads,
        params.filter_maxee,
        params.filter_minlen,
        params.filter_maxlen,
        params.filter_maxns 
    )
    //
    // MODULE: VSEARCH dereplicate per sample
    //
    VSEARCH_DEREPFULLLENGTH (
        VSEARCH_FASTQFILTER.out.reads,
        params.derep_strand
    )
    //
    // Combine all reads
    //
    fastq_files = VSEARCH_DEREPFULLLENGTH.out.reads
        .collect { it[1] }

    // 
    // MODULE: VSEARCH dereplicate for all reads
    //
    VSEARCH_DEREPFULLLENGTHALL (
        fastq_files,
        params.derep_all_strand,
        params.derep_all_fastawidth,
        params.derep_all_minunique
    )
    //
    // MODULE: VSEARCH cluster asvs
    //
    VSEARCH_CLUSTERUNOISE (
        VSEARCH_DEREPFULLLENGTHALL.out.reads,
        params.cluster_minsize,
        params.cluster_alpha
    )
    // 
    // MODULE: VSEARCH sort and remove singletons
    //
    VSEARCH_SORT_REMOVE_SINGLETONS (
        VSEARCH_CLUSTERUNOISE.out.asvs,
        params.sort_fastawidth,
        params.sort_minsize
    )
    //
    // MODULE: VSEARCH chimera detection
    //
    VSEARCH_UCHIMEDENOVO (
        VSEARCH_SORT_REMOVE_SINGLETONS.out.asvs,
        params.uchime_label
    )
    //
    // MODULE: VSEARCH make count table
    //
    VSEARCH_USEARCHGLOBAL (
        VSEARCH_DEREPFULLLENGTHALL.out.concatreads,
        VSEARCH_UCHIMEDENOVO.out.asvs,
        params.usearch_id
    )
    ch_versions = ch_versions.mix(VSEARCH_USEARCHGLOBAL.out.versions)
    
    if(params.skip_tree != true){
        //
        // MODULE: MAFFT for multiple sequence alignment
        //
        MAFFT (
            VSEARCH_UCHIMEDENOVO.out.asvs
        )
        ch_versions = ch_versions.mix(MAFFT.out.versions)

        //
        // MODULE: Build tree with FastTree
        //
            FASTTREE (
                MAFFT.out.msa
            )
            ch_versions = ch_versions.mix(FASTTREE.out.versions)
            ch_tree = FASTTREE.out.tree
    } else {
        ch_tree = Channel.fromPath("$projectDir/assets/NO_TREEFILE")
    }
    
    // 
    // MODULE: Download SILVA if not already present in db folder
    //
    SILVADATABASES()
    
    // 
    // MODULE: DADA2 Assign taxonomy with SILVA db
    //
    DADA2_ASSIGNTAXONOMY (
        VSEARCH_UCHIMEDENOVO.out.asvs,
        SILVADATABASES.out.asvdb,
        SILVADATABASES.out.speciesdb,
        params.dada2_minboot,
        params.dada2_allowmultiple,
        params.dada2_tryrevcompl
    )
    ch_versions = ch_versions.mix(DADA2_ASSIGNTAXONOMY.out.versions)

    //
    // MODULE: PICRUST
    //
    if(params.skip_picrust != true){
        PICRUST2 ( 
            VSEARCH_UCHIMEDENOVO.out.asvs,
            VSEARCH_USEARCHGLOBAL.out.counts
        )
        ch_versions = ch_versions.mix(PICRUST2.out.versions)
    }
    //
    // MODULE: Make phyloseq object
    //
    PHYLOSEQ_COMPLETE_MAKEOBJECT (
        VSEARCH_UCHIMEDENOVO.out.asvs,
        VSEARCH_USEARCHGLOBAL.out.counts,
        ch_tree,
        DADA2_ASSIGNTAXONOMY.out.taxtable
    )

    ch_phyloseq = PHYLOSEQ_COMPLETE_MAKEOBJECT.out.phyloseq
    ch_versions = ch_versions.mix(PHYLOSEQ_COMPLETE_MAKEOBJECT.out.versions)
    ch_taxtable = PHYLOSEQ_COMPLETE_MAKEOBJECT.out.taxtable
    ch_complete = true

    //
    // MODULE: Decontam (optional) — run before rarefaction
    //
    if (params.run_decontam) {
        PHYLOSEQ_DECONTAM (
            ch_phyloseq,
            file(params.input)
        )
        ch_versions = ch_versions.mix(PHYLOSEQ_DECONTAM.out.versions)
        // Use decontam output as basis for rarefaction and downstream steps
        ch_phyloseq_for_rarefaction = PHYLOSEQ_DECONTAM.out.phyloseq
    } else {
        ch_phyloseq_for_rarefaction = ch_phyloseq
    }

    //
    // MODULE: Fix taxonomy
    //
    if (!params.skip_fixtaxonomy) {
        PHYLOSEQ_COMPLETE_FIXTAX (
            ch_phyloseq,
            ch_complete
        )
        ch_taxtable = PHYLOSEQ_COMPLETE_FIXTAX.out.taxonomy
        // //
        // // MODULE: Overview metrics
        // //
        if (!params.skip_metrics) {
            PHYLOSEQ_COMPLETE_METRICS (
                ch_phyloseq,
                ch_taxtable,
                ch_complete
            )
        }
    }

    if (!params.skip_rarefaction) {
        ch_complete_new = false
        //
        // MODULE: Rarefaction — uses decontam output if run_decontam, else complete phyloseq
        //
        PHYLOSEQ_RAREFIED (
            ch_phyloseq_for_rarefaction,
            params.rarelevel,
        )
        
        ch_rarefied_phyloseq = PHYLOSEQ_RAREFIED.out.phyloseq

        if (!params.skip_fixtaxonomy) {
            PHYLOSEQ_RAREFIED_FIXTAX (
                ch_rarefied_phyloseq,
                ch_complete_new
            )
    
        ch_rarefied_taxtable = PHYLOSEQ_RAREFIED_FIXTAX.out.taxonomy
            if (!params.skip_metrics) {
                PHYLOSEQ_RAREFIED_METRICS (
                    ch_rarefied_phyloseq,
                    ch_rarefied_taxtable,
                    ch_complete_new
                )
            }
        }
    }
    
    //
    // MODULE: Collect software versions
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowVsearchpipeline.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowVsearchpipeline.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )
    multiqc_report = MULTIQC.out.report.toList()


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
