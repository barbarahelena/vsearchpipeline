/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap          } from 'plugin/nf-validation'
include { paramsSummaryMultiqc      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_vsearchpipeline_pipeline'

//
// MODULE: Local modules
//
include { SEQTK_TRIMFQ                                          } from '../modules/local/seqtk/trimfq'
include { VSEARCH_FASTQMERGEPAIRS                               } from '../modules/local/vsearch/fastqmergepairs'
include { VSEARCH_FASTQFILTER                                   } from '../modules/local/vsearch/fastqfilter'
include { VSEARCH_DEREPFULLLENGTH                               } from '../modules/local/vsearch/derepfulllength'
include { VSEARCH_DEREPFULLLENGTHALL                            } from '../modules/local/vsearch/derepfulllengthall'
include { VSEARCH_CLUSTERUNOISE                                 } from '../modules/local/vsearch/clusterunoise'
include { VSEARCH_SORT_REMOVE_SINGLETONS                        } from '../modules/local/vsearch/sort_remove_singletons'
include { VSEARCH_UCHIMEDENOVO                                  } from '../modules/local/vsearch/uchimedenovo'
include { VSEARCH_USEARCHGLOBAL                                 } from '../modules/local/vsearch/usearchglobal'
include { VSEARCH_MAPPINGRATE                                   } from '../modules/local/vsearch/mappingrate'
include { MAFFT                                                 } from '../modules/local/mafft'
include { FASTTREE                                              } from '../modules/local/fasttree'
include { SILVADATABASES                                        } from '../modules/local/silvadatabases'
include { DADA2_ASSIGNTAXONOMY                                  } from '../modules/local/dada2/assigntaxonomy'
include { PICRUST2                                              } from '../modules/local/picrust2'
include { PHYLOSEQ_MAKEOBJECT as PHYLOSEQ_COMPLETE_MAKEOBJECT   } from '../modules/local/phyloseq/makeobject'
include { PHYLOSEQ_FIXTAXONOMY as PHYLOSEQ_COMPLETE_FIXTAX      } from '../modules/local/phyloseq/fixtaxonomy'
include { PHYLOSEQ_METRICS as PHYLOSEQ_COMPLETE_METRICS         } from '../modules/local/phyloseq/metrics'
include { PHYLOSEQ_DECONTAM                                     } from '../modules/local/phyloseq/decontam'
include { PHYLOSEQ_RAREFACTION as PHYLOSEQ_RAREFIED             } from '../modules/local/phyloseq/rarefaction'
include { PHYLOSEQ_METRICS as PHYLOSEQ_RAREFIED_METRICS         } from '../modules/local/phyloseq/metrics'
include { PHYLOSEQ_FIXTAXONOMY as PHYLOSEQ_RAREFIED_FIXTAX      } from '../modules/local/phyloseq/fixtaxonomy'

//
// MODULE: nf-core modules
//
include { FASTQC  } from '../modules/nf-core/fastqc/main'
include { MULTIQC } from '../modules/nf-core/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VSEARCHPIPELINE {

    take:
    ch_reads   // channel: [ val(meta), [ reads ] ]
    ch_primers // channel: val(map with forward/reverse keys), or empty

    main:
    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_reads
    )
    //
    // MODULE:  Seqtk trim primers in fastq files
    //
    if (!params.skip_primers) {
        SEQTK_TRIMFQ (
            ch_reads, 
            ch_primers
        ).reads.set {ch_trimmed_reads}
    } else {
        ch_trimmed_reads = ch_reads
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
        .collect { f -> f[1] }

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

    //
    // MODULE: Calculate per-sample mapping rates from usearch_global outputs
    //
    VSEARCH_MAPPINGRATE (
        VSEARCH_USEARCHGLOBAL.out.counts,
        VSEARCH_USEARCHGLOBAL.out.mapping_stats,
        VSEARCH_FASTQFILTER.out.filter_stats.map { _meta, f -> f }.collect()
    )
    
    if(params.skip_tree != true){
        // MODULE: MAFFT for multiple sequence alignment
        //
        MAFFT (
            VSEARCH_UCHIMEDENOVO.out.asvs
        )

        //
        // MODULE: Build tree with FastTree
        //
            FASTTREE (
                MAFFT.out.msa
            )
            ch_tree = FASTTREE.out.tree
    } else {
        ch_tree = channel.fromPath("$projectDir/assets/NO_TREEFILE")
    }
    
    // 
    // MODULE: Download SILVA if not already present in db folder, or use supplied paths
    //
    if (params.silva_asv_db && params.silva_species_db) {
        ch_silva_asv_db     = channel.fromPath(params.silva_asv_db, checkIfExists: true)
        ch_silva_species_db = channel.fromPath(params.silva_species_db, checkIfExists: true)
    } else {
        SILVADATABASES()
        ch_silva_asv_db     = SILVADATABASES.out.asvdb
        ch_silva_species_db = SILVADATABASES.out.speciesdb
    }
    
    // 
    // MODULE: DADA2 Assign taxonomy with SILVA db
    //
    DADA2_ASSIGNTAXONOMY (
        VSEARCH_UCHIMEDENOVO.out.asvs,
        ch_silva_asv_db,
        ch_silva_species_db,
        params.dada2_minboot,
        params.dada2_allowmultiple,
        params.dada2_tryrevcompl
    )

    //
    // MODULE: PICRUST
    //
    if(params.skip_picrust != true){
        PICRUST2 ( 
            VSEARCH_UCHIMEDENOVO.out.asvs,
            VSEARCH_USEARCHGLOBAL.out.counts
        )
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
    // MODULE: MultiQC
    //
    def ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    def workflow_summary    = paramsSummaryMultiqc(paramsSummaryMap(workflow))
    def ch_workflow_summary = channel.value(workflow_summary).collectFile(name: 'workflow_summary_mqc.yaml')

    def methods_description = methodsDescriptionText(ch_multiqc_custom_methods_description)
    def ch_methods_description = channel.value(methods_description).collectFile(name: 'methods_description_mqc.yaml')

    def n_samples = ch_reads.count().map { n -> n * 2 }  // paired-end: 2 FastQC files per sample

    n_samples.subscribe { n ->
        if (n > 6000) {
            log.warn "[vsearchpipeline] Skipping MultiQC: number of FastQC files (${n}) exceeds 6000."
        }
    }

    ch_multiqc_files = channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary)
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description)
    ch_multiqc_files = ch_multiqc_files.mix(channel.topic('versions').collectFile(name: 'collated_versions.yml'))
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map { _meta, zips -> zips }.flatten())

    n_samples.filter { n -> n <= 6000 }
        .map { n -> [ n ] }
        .combine(ch_multiqc_files.collect().map { files -> [ files ] })
        .map { _n, files -> [ [id: 'multiqc'], files, [], [], [], [] ] }
        .set { ch_multiqc_input }

    MULTIQC ( ch_multiqc_input )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
