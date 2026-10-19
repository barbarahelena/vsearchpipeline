nextflow.enable.types = true

process PHYLOSEQ_FIXTAXONOMY {
    label 'process_multi_low'
    label 'phyloseq'

    input:
    phyloseq: Path
    postfix: String   // 'complete' or 'rarefied', used in output file names

    output:
    taxonomy     = file("taxtable_${postfix}.RDS")
    phylevels    = file("phylogen_levels_${postfix}.csv")
    phylevelstop = file("phylogen_levels_top300_${postfix}.csv")

    topic:
    tuple(task.process, 'r-base',   file("version_r.txt"))        >> 'versions'
    tuple(task.process, 'phyloseq', file("version_phyloseq.txt")) >> 'versions'

    script:
    template 'fixtaxonomy.R'

    stub:
    """
    touch taxtable_${postfix}.RDS
    touch phylogen_levels_${postfix}.csv
    touch phylogen_levels_top300_${postfix}.csv
    touch version_r.txt
    touch version_phyloseq.txt
    """
}
