nextflow.enable.types = true

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
    SAMPLES
    One record per samplesheet row. Per-sample steps take only the fields they need
    and return a record with the sample id and their new file(s).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record Sample {
    id: String
    single_end: Boolean
    fastq: List<Path>   // [ fastq_1, fastq_2 ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OPTIONS
    The parameters this workflow uses. main.nf passes `params`; only these fields are
    read, and each must have the declared type.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

record VsearchOptions {
    input: Path?
    outdir: String
    skip_primers: Boolean
    merge_allowmergestagger: Boolean
    merge_maxdiffs: Integer
    merge_minlen: Integer
    merge_maxlen: Integer
    merge_maxdiffpct: Integer
    filter_maxee: Float
    filter_minlen: Integer
    filter_maxlen: Integer
    filter_maxns: Integer
    derep_strand: String
    derep_all_strand: String
    derep_all_fastawidth: Integer
    derep_all_minunique: Integer
    cluster_minsize: Integer
    cluster_alpha: Float
    sort_fastawidth: Integer
    sort_minsize: Integer
    uchime_label: String
    usearch_id: Float
    skip_tree: Boolean
    silva_asv_db: Path?
    silva_species_db: Path?
    save_silva_db: Boolean
    dada2_minboot: Integer
    dada2_allowmultiple: Integer
    dada2_tryrevcompl: Boolean
    skip_picrust: Boolean
    run_decontam: Boolean
    skip_fixtaxonomy: Boolean
    skip_metrics: Boolean
    skip_rarefaction: Boolean
    rarelevel: Integer
    multiqc_methods_description: String?
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VSEARCHPIPELINE {

    take:
    ch_samples: Channel<Sample>
    ch_primers: Value<Map>     // [ forward: ..., reverse: ... ]; [:] when --skip_primers
    opts: VsearchOptions

    main:
    //
    // MODULE: Run FastQC
    //
    // nf-core FastQC is untyped and expects tuple(meta, reads)
    def fastqc_results = FASTQC(ch_samples.map { s -> tuple([id: s.id, single_end: s.single_end], s.fastq) })   // not 'fastqc': an emit name below

    //
    // MODULE: Seqtk trim primers in fastq files
    //
    // record(...) keeps the field types; subMap() would return an untyped Record
    def ch_trimmed_reads = ch_samples.map { s -> record(id: s.id, fastq: s.fastq) }
    def ch_out_seqtk     = channel.empty()
    if (!opts.skip_primers) {
        def trimmed = SEQTK_TRIMFQ(ch_samples, ch_primers)
        ch_trimmed_reads = trimmed
        ch_out_seqtk     = trimmed.flatMap { r -> r.fastq }
    }

    //
    // MODULE: VSEARCH merge read pairs, quality filter, dereplicate per sample
    //
    def merged = VSEARCH_FASTQMERGEPAIRS(
        ch_trimmed_reads,
        opts.merge_allowmergestagger,
        opts.merge_maxdiffs,
        opts.merge_minlen,
        opts.merge_maxlen,
        opts.merge_maxdiffpct
    )
    def filtered = VSEARCH_FASTQFILTER(
        merged,
        opts.filter_maxee,
        opts.filter_minlen,
        opts.filter_maxlen,
        opts.filter_maxns
    )
    def derep = VSEARCH_DEREPFULLLENGTH(
        filtered,
        opts.derep_strand
    )

    //
    // MODULE: VSEARCH dereplicate all samples together
    //
    // Sort by file name so the input order (and so the ASV numbering) does not depend on which
    // task finished first. collect() returns an unordered Bag; toSorted gives an ordered List.
    def ch_derep_sorted = derep
        .map { r -> r.derep }
        .collect()
        .map { fastas -> fastas.toSorted { a, b -> a.name <=> b.name } }
    def derep_all = VSEARCH_DEREPFULLLENGTHALL(
        ch_derep_sorted,
        opts.derep_all_strand,
        opts.derep_all_fastawidth,
        opts.derep_all_minunique
    )

    //
    // MODULE: VSEARCH denoise (UNOISE), remove singletons and chimeras, count table, mapping rate
    //
    def clustered = VSEARCH_CLUSTERUNOISE(
        derep_all.reads,
        opts.cluster_minsize,
        opts.cluster_alpha
    )
    def nonsingle = VSEARCH_SORT_REMOVE_SINGLETONS(
        clustered,
        opts.sort_fastawidth,
        opts.sort_minsize
    )
    def nonchimeric = VSEARCH_UCHIMEDENOVO(
        nonsingle,
        opts.uchime_label
    )
    def ch_asvs = nonchimeric.asvs
    def counts = VSEARCH_USEARCHGLOBAL(
        derep_all.concatreads,
        ch_asvs,
        opts.usearch_id
    )
    def mappingrate = VSEARCH_MAPPINGRATE(
        counts.counts,
        counts.mapping_stats,
        filtered.map { r -> r.filter_stats }.collect().map { files -> files.toList() }
    )

    // mix() takes one channel or value at a time
    def ch_out_vsearch = channel.empty()
        .mix(merged.map { r -> r.merged })
        .mix(filtered.map { r -> r.filtered })
        .mix(filtered.map { r -> r.filter_stats })
        .mix(derep.map { r -> r.derep })
        .mix(derep_all.reads)
        .mix(derep_all.concatreads)
        .mix(clustered)
        .mix(nonsingle)
        .mix(nonchimeric.asvs)
        .mix(nonchimeric.chimeras)
        .mix(counts.counts)
        .mix(counts.mapping_stats)
        .mix(mappingrate.summary)
        .mix(mappingrate.overall)

    //
    // MODULE: MAFFT multiple sequence alignment + FastTree
    //
    def ch_tree         = channel.value(file("${projectDir}/assets/NO_TREEFILE"))
    def ch_out_mafft    = channel.empty()
    def ch_out_fasttree = channel.empty()
    if (!opts.skip_tree) {
        def msa  = MAFFT(ch_asvs)
        def tree = FASTTREE(msa)
        ch_tree         = tree
        ch_out_mafft    = channel.empty().mix(msa)    // value -> channel, for publishing
        ch_out_fasttree = channel.empty().mix(tree)
    }

    //
    // MODULE: Download SILVA if not supplied
    //
    def ch_silva_asv_db     = channel.value(opts.silva_asv_db)
    def ch_silva_species_db = channel.value(opts.silva_species_db)
    def ch_out_silva        = channel.empty()
    if (!opts.silva_asv_db || !opts.silva_species_db) {
        def silva = SILVADATABASES()
        ch_silva_asv_db     = silva.asvdb
        ch_silva_species_db = silva.speciesdb
        if (opts.save_silva_db) {
            ch_out_silva = channel.empty().mix(silva.asvdb).mix(silva.speciesdb)
        }
    }

    //
    // MODULE: DADA2 Assign taxonomy with SILVA db
    //
    def taxonomy = DADA2_ASSIGNTAXONOMY(
        ch_asvs,
        ch_silva_asv_db,
        ch_silva_species_db,
        opts.dada2_minboot,
        opts.dada2_allowmultiple,
        opts.dada2_tryrevcompl
    )

    //
    // MODULE: PICRUST
    //
    def ch_out_picrust2 = channel.empty()
    if (!opts.skip_picrust) {
        def picrust = PICRUST2(ch_asvs, counts.counts)
        ch_out_picrust2 = picrust.flatMap { files -> files }
    }

    //
    // MODULE: Make phyloseq object
    //
    def phyloseq = PHYLOSEQ_COMPLETE_MAKEOBJECT(
        ch_asvs,
        counts.counts,
        ch_tree,
        taxonomy
    )
    def ch_out_phyloseq_complete = channel.empty().mix(phyloseq.phyloseq).mix(phyloseq.taxtable)

    //
    // MODULE: Decontam (optional) — run before rarefaction
    //
    def ch_phyloseq_for_rarefaction = phyloseq.phyloseq
    def ch_out_phyloseq_decontam    = channel.empty()
    if (opts.run_decontam) {
        def decontam = PHYLOSEQ_DECONTAM(phyloseq.phyloseq, opts.input)
        // Use decontam output as basis for rarefaction and downstream steps
        ch_phyloseq_for_rarefaction = decontam.phyloseq
        ch_out_phyloseq_decontam    = channel.empty()
            .mix(decontam.phyloseq).mix(decontam.report).mix(decontam.contaminants).mix(decontam.prev_plot)
    }

    //
    // MODULE: Fix taxonomy + overview metrics (complete dataset)
    //
    if (!opts.skip_fixtaxonomy) {
        def fixtax = PHYLOSEQ_COMPLETE_FIXTAX(phyloseq.phyloseq, 'complete')
        ch_out_phyloseq_complete = ch_out_phyloseq_complete
            .mix(fixtax.taxonomy).mix(fixtax.phylevels).mix(fixtax.phylevelstop)

        if (!opts.skip_metrics) {
            def metrics = PHYLOSEQ_COMPLETE_METRICS(phyloseq.phyloseq, fixtax.taxonomy, 'complete')
            ch_out_phyloseq_complete = ch_out_phyloseq_complete
                .mix(metrics.speciescomp).mix(metrics.genuscomp).mix(metrics.famcomp).mix(metrics.phylumcomp)
                .mix(metrics.shannon).mix(metrics.richness).mix(metrics.metrics)
        }
    }

    //
    // MODULE: Rarefaction — uses decontam output if run_decontam, else complete phyloseq
    //
    def ch_out_phyloseq_rarefied = channel.empty()
    if (!opts.skip_rarefaction) {
        def rarefied = PHYLOSEQ_RAREFIED(ch_phyloseq_for_rarefaction, opts.rarelevel)
        ch_out_phyloseq_rarefied = channel.empty().mix(rarefied.phyloseq).mix(rarefied.rarecurve).mix(rarefied.rarereport)

        if (!opts.skip_fixtaxonomy) {
            def fixtax_r = PHYLOSEQ_RAREFIED_FIXTAX(rarefied.phyloseq, 'rarefied')
            ch_out_phyloseq_rarefied = ch_out_phyloseq_rarefied
                .mix(fixtax_r.taxonomy).mix(fixtax_r.phylevels).mix(fixtax_r.phylevelstop)

            if (!opts.skip_metrics) {
                def metrics_r = PHYLOSEQ_RAREFIED_METRICS(rarefied.phyloseq, fixtax_r.taxonomy, 'rarefied')
                ch_out_phyloseq_rarefied = ch_out_phyloseq_rarefied
                    .mix(metrics_r.speciescomp).mix(metrics_r.genuscomp).mix(metrics_r.famcomp).mix(metrics_r.phylumcomp)
                    .mix(metrics_r.shannon).mix(metrics_r.richness).mix(metrics_r.metrics)
            }
        }
    }

    //
    // Record the pipeline git commit. Written by Nextflow itself: no task, container or git in a container needed.
    //
    def ch_commit_hash = channel.of("""
        You can obtain the exact version of the vsearchpipeline used for this run
        by cloning the repo and running the Git command:

        git checkout ${gitCommitHash()}
        """.stripIndent())
        .collectFile(name: 'git_commit_hash.txt')

    //
    // Collect software versions from the 'versions' topic channel.
    // Each entry is [ process, tool, version ]: the version is a string (eval) or a file (R templates).
    //
    def ch_versions_topic: Channel<Tuple> = channel.topic('versions')
    def ch_versions = ch_versions_topic
        .map { process, tool, version ->
            def v = (version instanceof Path) ? version.text.trim() : version.toString().trim()
            tuple(process.tokenize(':')[-1], tool, v)   // tuple(), not a list: the type checker then knows its 3 parts
        }
        .unique()

    // Plain table (always written, independent of MultiQC)
    ch_versions
        .map { process, tool, version -> "${process}\t${tool}\t${version}" }
        .collectFile(
            name: 'software_versions.tsv',
            storeDir: "${opts.outdir}/pipeline_info",
            seed: "process\tsoftware\tversion",
            sort: true,
            newLine: true
        )

    // YAML grouped per process, as MultiQC's software versions section expects.
    // Map to strings *before* collectFile: given tuples, collectFile uses the first element as a file name.
    def ch_versions_yaml = ch_versions
        .map { process, tool, version -> tuple(process, "    ${tool}: ${version}") }
        .groupBy()
        .map { process, lines -> "${process}:\n${lines.toSorted().join('\n')}" }
        .mix(channel.of(workflowVersionToYAML()))
        .collectFile(
            name: 'vsearchpipeline_software_mqc_versions.yml',
            storeDir: "${opts.outdir}/pipeline_info",
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    def methods_description_template = opts.multiqc_methods_description
        ? file(opts.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)

    def ch_workflow_summary = channel.of(paramsSummaryMultiqc(paramsSummaryMap()))
        .collectFile(name: 'workflow_summary_mqc.yaml')
    def ch_methods_description = channel.of(methodsDescriptionText(methods_description_template))
        .collectFile(name: 'methods_description_mqc.yaml')

    // Skip MultiQC for very large runs (> 6000 FastQC reports): it becomes too slow to be useful
    def ch_multiqc_input = channel.empty()
        .mix(ch_workflow_summary)
        .mix(ch_methods_description)
        .mix(ch_versions_yaml)
        .mix(fastqc_results.zip.flatMap { _meta, zips -> [ zips ].flatten() })   // one file (single-end) or a list
        .collect()
        .flatMap { files ->
            def n_fastqc = files.findAll { f -> f.name.endsWith('_fastqc.zip') }.size()
            if (n_fastqc > 6000) {
                log.warn "[vsearchpipeline] Skipping MultiQC: number of FastQC files (${n_fastqc}) exceeds 6000."
                return []
            }
            return [ tuple([id: 'multiqc'], files.toList(), [], [], [], []) ]
        }

    def multiqc_results = MULTIQC(ch_multiqc_input)

    emit:
    // One channel of files per results folder; the output block in main.nf decides where each goes
    fastqc: Channel<Path>            = fastqc_results.html.mix(fastqc_results.zip).flatMap { _meta, files -> [ files ].flatten() }
    seqtk: Channel<Path>             = ch_out_seqtk
    vsearch: Channel<Path>           = ch_out_vsearch
    mafft: Channel<Path>             = ch_out_mafft
    fasttree: Channel<Path>          = ch_out_fasttree
    silva_db: Channel<Path>          = ch_out_silva
    dada2: Channel<Path>             = channel.empty().mix(taxonomy)
    picrust2: Channel<Path>          = ch_out_picrust2
    phyloseq_complete: Channel<Path> = ch_out_phyloseq_complete
    phyloseq_decontam: Channel<Path> = ch_out_phyloseq_decontam
    phyloseq_rarefied: Channel<Path> = ch_out_phyloseq_rarefied
    multiqc: Channel<Path>           = multiqc_results.report.mix(multiqc_results.data).mix(multiqc_results.plots).map { _meta, f -> f }
    pipeline_info: Channel<Path>     = ch_commit_hash
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
