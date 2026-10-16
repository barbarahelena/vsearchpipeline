nextflow.enable.types = true

process VSEARCH_USEARCHGLOBAL {
    label 'process_highcpu'
    label 'error_retry'
    label 'vsearch'

    input:
    allreads: Path   // all reads concatenated
    asvs: Path
    id: Float        // minimum identity, e.g. 0.97

    output:
    counts        = file("count_table.txt")
    mapping_stats = file("mapping_stats.txt")

    topic:
    tuple(task.process, 'vsearch', eval('vsearch --version 2>&1 | sed -nE \'1s/^vsearch v([0-9.]+).*/\\1/p\'')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    """
    vsearch \\
        --usearch_global ${allreads} \\
        --db ${asvs} \\
        --id ${id} \\
        --threads ${task.cpus} \\
        ${args} \\
        --otutabout count_table.txt \\
        2> mapping_stats.txt
    """

    stub:
    """
    touch count_table.txt
    touch mapping_stats.txt
    """
}
