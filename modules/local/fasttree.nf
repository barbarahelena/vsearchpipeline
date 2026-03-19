process FASTTREE {
    label 'process_single_medium'
    label 'error_retry'
    label 'fasttree'

    input:
    path msa

    output:
    path "asvs.msa.tree"                                                                                                                            , emit: tree
    tuple val("${task.process}"), val('fasttree'), eval("fasttree -help 2>&1 | head -n 1 | sed -n 's/.*\\([0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\).*/\\1/p'"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    fasttree \\
        -nt \\
        -gtr \\
        -gamma \\
        $msa \\
        > asvs.msa.tree
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch asvs.msa.tree
    """
}
