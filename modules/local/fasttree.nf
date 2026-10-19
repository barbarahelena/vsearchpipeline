nextflow.enable.types = true

process FASTTREE {
    label 'process_single_medium'
    label 'error_retry'
    label 'fasttree'

    input:
    msa: Path

    output:
    tree = file("asvs.msa.tree")

    topic:
    tuple(task.process, 'fasttree', eval('fasttree -help 2>&1 | sed -nE \'1s/^FastTree ([0-9.]+).*/\\1/p\'')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    """
    fasttree \\
        -nt \\
        -gtr \\
        -gamma \\
        ${args} \\
        ${msa} \\
        > asvs.msa.tree
    """

    stub:
    """
    touch asvs.msa.tree
    """
}
