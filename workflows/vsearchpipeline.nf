/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { workflowVersionToYAML     } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_vsearchpipeline_pipeline'
include { gitCommitHash             } from '../subworkflows/local/utils_vsearchpipeline_pipeline'

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
include { FASTQC          } from '../modules/nf-core/fastqc/main'
include { MULTIQC         } from '../modules/nf-core/multiqc/main'

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
    // MODULE: Seqtk trim primers in fastq files
    //
    def ch_trimmed_reads = ch_reads
    def ch_out_seqtk     = channel.empty()
    if (!params.skip_primers) {
        SEQTK_TRIMFQ (
            ch_reads,
            ch_primers
        )
        ch_trimmed_reads = SEQTK_TRIMFQ.out.reads
        ch_out_seqtk     = SEQTK_TRIMFQ.out.reads
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
    // MODULE: VSEARCH dereplicate for all reads
    //
    VSEARCH_DEREPFULLLENGTHALL (
        // Sort by file name so the input order (and so the ASV numbering) does not depend on which task
        // finished first. Not sort: true, which compares full paths including the random work dir names.
        VSEARCH_DEREPFULLLENGTH.out.reads.map { _meta, fasta -> fasta }.collect(sort: { a, b -> a.name <=> b.name }),
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

    def ch_out_vsearch = channel.empty()
        .mix(
            VSEARCH_FASTQMERGEPAIRS.out.reads,
            VSEARCH_FASTQFILTER.out.reads,
            VSEARCH_FASTQFILTER.out.filter_stats,
            VSEARCH_DEREPFULLLENGTH.out.reads,
            VSEARCH_DEREPFULLLENGTHALL.out.reads,
            VSEARCH_DEREPFULLLENGTHALL.out.concatreads,
            VSEARCH_CLUSTERUNOISE.out.asvs,
            VSEARCH_SORT_REMOVE_SINGLETONS.out.asvs,
            VSEARCH_UCHIMEDENOVO.out.asvs,
            VSEARCH_UCHIMEDENOVO.out.chimeras,
            VSEARCH_USEARCHGLOBAL.out.counts,
            VSEARCH_USEARCHGLOBAL.out.mapping_stats,
            VSEARCH_MAPPINGRATE.out.summary,
            VSEARCH_MAPPINGRATE.out.overall
        )

    //
    // MODULE: MAFFT multiple sequence alignment + FastTree
    //
    def ch_tree          = channel.fromPath("${projectDir}/assets/NO_TREEFILE")
    def ch_out_mafft     = channel.empty()
    def ch_out_fasttree  = channel.empty()
    if (!params.skip_tree) {
        MAFFT (
            VSEARCH_UCHIMEDENOVO.out.asvs
        )
        FASTTREE (
            MAFFT.out.msa
        )
        ch_tree         = FASTTREE.out.tree
        ch_out_mafft    = MAFFT.out.msa
        ch_out_fasttree = FASTTREE.out.tree
    }

    //
    // MODULE: Download SILVA if not supplied
    //
    def ch_silva_asv_db     = channel.empty()
    def ch_silva_species_db = channel.empty()
    def ch_out_silva        = channel.empty()
    if (params.silva_asv_db && params.silva_species_db) {
        ch_silva_asv_db     = channel.fromPath(params.silva_asv_db, checkIfExists: true)
        ch_silva_species_db = channel.fromPath(params.silva_species_db, checkIfExists: true)
    } else {
        SILVADATABASES()
        ch_silva_asv_db     = SILVADATABASES.out.asvdb
        ch_silva_species_db = SILVADATABASES.out.speciesdb
        if (params.save_silva_db) {
            ch_out_silva = SILVADATABASES.out.asvdb.mix(SILVADATABASES.out.speciesdb)
        }
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
    def ch_out_picrust2 = channel.empty()
    if (!params.skip_picrust) {
        PICRUST2 (
            VSEARCH_UCHIMEDENOVO.out.asvs,
            VSEARCH_USEARCHGLOBAL.out.counts
        )
        ch_out_picrust2 = PICRUST2.out.outfolder
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

    def ch_phyloseq = PHYLOSEQ_COMPLETE_MAKEOBJECT.out.phyloseq
    def ch_out_phyloseq_complete = PHYLOSEQ_COMPLETE_MAKEOBJECT.out.phyloseq
        .mix(PHYLOSEQ_COMPLETE_MAKEOBJECT.out.taxtable)

    //
    // MODULE: Decontam (optional) — run before rarefaction
    //
    def ch_phyloseq_for_rarefaction = ch_phyloseq
    def ch_out_phyloseq_decontam    = channel.empty()
    if (params.run_decontam) {
        PHYLOSEQ_DECONTAM (
            ch_phyloseq,
            file(params.input)
        )
        // Use decontam output as basis for rarefaction and downstream steps
        ch_phyloseq_for_rarefaction = PHYLOSEQ_DECONTAM.out.phyloseq
        ch_out_phyloseq_decontam    = PHYLOSEQ_DECONTAM.out.phyloseq
            .mix(
                PHYLOSEQ_DECONTAM.out.report,
                PHYLOSEQ_DECONTAM.out.contaminants,
                PHYLOSEQ_DECONTAM.out.prev_plot
            )
    }

    //
    // MODULE: Fix taxonomy + overview metrics (complete dataset)
    //
    if (!params.skip_fixtaxonomy) {
        PHYLOSEQ_COMPLETE_FIXTAX (
            ch_phyloseq,
            'complete'
        )
        ch_out_phyloseq_complete = ch_out_phyloseq_complete.mix(
            PHYLOSEQ_COMPLETE_FIXTAX.out.taxonomy,
            PHYLOSEQ_COMPLETE_FIXTAX.out.phylevels,
            PHYLOSEQ_COMPLETE_FIXTAX.out.phylevelstop
        )

        if (!params.skip_metrics) {
            PHYLOSEQ_COMPLETE_METRICS (
                ch_phyloseq,
                PHYLOSEQ_COMPLETE_FIXTAX.out.taxonomy,
                'complete'
            )
            ch_out_phyloseq_complete = ch_out_phyloseq_complete.mix(
                PHYLOSEQ_COMPLETE_METRICS.out.speciescomp,
                PHYLOSEQ_COMPLETE_METRICS.out.genuscomp,
                PHYLOSEQ_COMPLETE_METRICS.out.famcomp,
                PHYLOSEQ_COMPLETE_METRICS.out.phylumcomp,
                PHYLOSEQ_COMPLETE_METRICS.out.shannon,
                PHYLOSEQ_COMPLETE_METRICS.out.richness,
                PHYLOSEQ_COMPLETE_METRICS.out.metrics
            )
        }
    }

    //
    // MODULE: Rarefaction — uses decontam output if run_decontam, else complete phyloseq
    //
    def ch_out_phyloseq_rarefied = channel.empty()
    if (!params.skip_rarefaction) {
        PHYLOSEQ_RAREFIED (
            ch_phyloseq_for_rarefaction,
            params.rarelevel
        )
        ch_out_phyloseq_rarefied = PHYLOSEQ_RAREFIED.out.phyloseq
            .mix(PHYLOSEQ_RAREFIED.out.rarecurve, PHYLOSEQ_RAREFIED.out.rarereport)

        if (!params.skip_fixtaxonomy) {
            PHYLOSEQ_RAREFIED_FIXTAX (
                PHYLOSEQ_RAREFIED.out.phyloseq,
                'rarefied'
            )
            ch_out_phyloseq_rarefied = ch_out_phyloseq_rarefied.mix(
                PHYLOSEQ_RAREFIED_FIXTAX.out.taxonomy,
                PHYLOSEQ_RAREFIED_FIXTAX.out.phylevels,
                PHYLOSEQ_RAREFIED_FIXTAX.out.phylevelstop
            )

            if (!params.skip_metrics) {
                PHYLOSEQ_RAREFIED_METRICS (
                    PHYLOSEQ_RAREFIED.out.phyloseq,
                    PHYLOSEQ_RAREFIED_FIXTAX.out.taxonomy,
                    'rarefied'
                )
                ch_out_phyloseq_rarefied = ch_out_phyloseq_rarefied.mix(
                    PHYLOSEQ_RAREFIED_METRICS.out.speciescomp,
                    PHYLOSEQ_RAREFIED_METRICS.out.genuscomp,
                    PHYLOSEQ_RAREFIED_METRICS.out.famcomp,
                    PHYLOSEQ_RAREFIED_METRICS.out.phylumcomp,
                    PHYLOSEQ_RAREFIED_METRICS.out.shannon,
                    PHYLOSEQ_RAREFIED_METRICS.out.richness,
                    PHYLOSEQ_RAREFIED_METRICS.out.metrics
                )
            }
        }
    }

    //
    // Record the pipeline git commit. Written by Nextflow itself: no task, container or git in a container needed.
    //
    def ch_commit_hash = channel.value("""
        You can obtain the exact version of the vsearchpipeline used for this run
        by cloning the repo and running the Git command:

        git checkout ${gitCommitHash()}
        """.stripIndent())
        .collectFile(name: 'git_commit_hash.txt')

    //
    // Collect software versions from the 'versions' topic channel.
    // Each entry is [ process, tool, version ]: the version is a string (eval) or a file (R templates).
    //
    def ch_versions = channel.topic('versions')
        .map { process, tool, version ->
            def v = (version instanceof Path) ? version.text.trim() : version.toString().trim()
            [ process.tokenize(':')[-1], tool, v ]
        }
        .unique()

    // Plain table (always written, independent of MultiQC)
    ch_versions
        .map { process, tool, version -> "${process}\t${tool}\t${version}" }
        .collectFile(
            name: 'software_versions.tsv',
            storeDir: "${params.outdir}/pipeline_info",
            seed: "process\tsoftware\tversion",
            sort: true,
            newLine: true
        )

    // YAML grouped per process, as MultiQC's software versions section expects.
    // Map to strings *before* collectFile: given tuples, collectFile uses the first element as a file name.
    def ch_versions_yaml = ch_versions
        .map { process, tool, version -> [ process, "    ${tool}: ${version}" ] }
        .groupTuple()
        .map { process, lines -> "${process}:\n${lines.sort().join('\n')}" }
        .mix(channel.of(workflowVersionToYAML()))
        .collectFile(
            name: 'vsearchpipeline_software_mqc_versions.yml',
            storeDir: "${params.outdir}/pipeline_info",
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    def ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    def workflow_summary    = paramsSummaryMultiqc(paramsSummaryMap())
    def ch_workflow_summary = channel.value(workflow_summary).collectFile(name: 'workflow_summary_mqc.yaml')

    def methods_description = methodsDescriptionText(ch_multiqc_custom_methods_description)
    def ch_methods_description = channel.value(methods_description).collectFile(name: 'methods_description_mqc.yaml')

    def n_samples = ch_reads.count().map { n -> n * 2 }  // paired-end: 2 FastQC files per sample

    n_samples.subscribe { n ->
        if (n > 6000) {
            log.warn "[vsearchpipeline] Skipping MultiQC: number of FastQC files (${n}) exceeds 6000."
        }
    }

    def ch_multiqc_files = channel.empty()
        .mix(
            ch_workflow_summary,
            ch_methods_description,
            ch_versions_yaml,
            FASTQC.out.zip.map { _meta, zips -> zips }.flatten()
        )

    def ch_multiqc_input = n_samples.filter { n -> n <= 6000 }
        .map { n -> [ n ] }
        .combine(ch_multiqc_files.collect().map { files -> [ files ] })
        .map { _n, files -> [ [id: 'multiqc'], files, [], [], [], [] ] }

    MULTIQC ( ch_multiqc_input )

    emit:
    // One channel per results folder; the output block in main.nf decides where each goes
    fastqc            = FASTQC.out.html.mix(FASTQC.out.zip)
    seqtk             = ch_out_seqtk
    vsearch           = ch_out_vsearch
    mafft             = ch_out_mafft
    fasttree          = ch_out_fasttree
    silva_db          = ch_out_silva
    dada2             = DADA2_ASSIGNTAXONOMY.out.taxtable
    picrust2          = ch_out_picrust2
    phyloseq_complete = ch_out_phyloseq_complete
    phyloseq_decontam = ch_out_phyloseq_decontam
    phyloseq_rarefied = ch_out_phyloseq_rarefied
    multiqc           = MULTIQC.out.report.mix(MULTIQC.out.data, MULTIQC.out.plots)
    pipeline_info     = ch_commit_hash
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
