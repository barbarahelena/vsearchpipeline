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
    path "versions.yml"             , emit: versions, topic: versions

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
    touch versions.yml
    """
}
