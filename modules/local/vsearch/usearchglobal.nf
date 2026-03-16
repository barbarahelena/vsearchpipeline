process VSEARCH_USEARCHGLOBAL {
    label 'process_highcpu'
    label 'error_retry'
    label 'vsearch'

    input:
    path(allreads)
    path(asvs)
    val id
    
    output:
    path "all.concat.fasta"                                                                                               , emit: concatfasta
    path "count_table.txt"                                                                                                , emit: counts
    tuple val("${task.process}"), val('vsearch'), eval('vsearch --version 2>&1 | head -n 1 | sed \'s/vsearch //; s/,.*//\''), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """  
    vsearch \\
        --usearch_global $allreads \\
        --db $asvs \\
        --id $id \\
        --threads $task.cpus \\
        --otutabout count_table.txt
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch all.concat.fasta
    touch count_table.txt
    """
}
