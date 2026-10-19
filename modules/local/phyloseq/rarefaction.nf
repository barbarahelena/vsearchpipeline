nextflow.enable.types = true

process PHYLOSEQ_RAREFACTION {
    label 'process_multi_low'
    label 'phyloseq'

    input:
    phyloseq: Path
    rarelevel: Integer   // 0 = choose a level automatically

    output:
    phyloseq   = file("phyloseq_rarefied.RDS")
    rarecurve  = file("rarefaction_plot.pdf")
    rarereport = file("rarefaction_report.txt")

    topic:
    tuple(task.process, 'r-base',   file("version_r.txt"))        >> 'versions'
    tuple(task.process, 'phyloseq', file("version_phyloseq.txt")) >> 'versions'
    tuple(task.process, 'ggplot2',  file("version_ggplot2.txt"))  >> 'versions'

    script:
    template 'rarefaction.R'

    stub:
    """
    touch phyloseq_rarefied.RDS
    touch rarefaction_plot.pdf
    touch rarefaction_report.txt
    touch version_r.txt
    touch version_phyloseq.txt
    touch version_ggplot2.txt
    """
}
