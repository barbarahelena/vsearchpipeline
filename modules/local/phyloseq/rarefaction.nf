process PHYLOSEQ_RAREFACTION {
    label 'process_multi_low'
    label 'phyloseq'

    input:
    path    phyloseq
    val     rarelevel

    output:
    path "phyloseq_rarefied.RDS"    , emit: phyloseq
    path "rarefaction_plot.pdf"     , emit: rarecurve
    path "rarefaction_report.txt"   , emit: rarereport
    tuple val("${task.process}"), val('r-base'),   path("version_r.txt"),        emit: versions_r,        topic: versions
    tuple val("${task.process}"), val('phyloseq'), path("version_phyloseq.txt"), emit: versions_phyloseq, topic: versions
    tuple val("${task.process}"), val('ggplot2'),  path("version_ggplot2.txt"),  emit: versions_ggplot2,  topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    nf_seed      = task.ext.seed ?: '1234'
    nf_rarelevel = rarelevel ?: 0
    template 'rarefaction.R'

    stub:
    nf_seed      = task.ext.seed ?: '1234'
    nf_rarelevel = rarelevel ?: 0
    """
    touch phyloseq_rarefied.RDS
    touch rarefaction_plot.pdf
    touch rarefaction_report.txt
    touch version_r.txt
    touch version_phyloseq.txt
    touch version_ggplot2.txt
    """
}
