process VSEARCH_SORT_REMOVE_SINGLETONS {
    label 'process_single_low'
    label 'vsearch'
    label 'error_retry'

    input:
    path(reads)
    val fasta_width
    val minsize

    output:
    path "asvs_nonsingle.fasta"                                                                                           , emit: asvs
    tuple val("${task.process}"), val('vsearch'), eval('vsearch --version 2>&1 | sed -nE \'1s/^vsearch v([0-9.]+).*/\\1/p\''), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    vsearch --sortbysize $reads \
        --threads $task.cpus \
        --sizein \
        --sizeout \
        --fasta_width $fasta_width \
        --minsize $minsize \
        --output asvs_nonsingle.fasta
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch asvs_nonsingle.fasta
    """
}
