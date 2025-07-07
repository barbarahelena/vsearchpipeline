process FASTTREE {
    label 'process_single_medium'
    label 'error_retry'

    conda "bioconda::fasttree=2.1.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fasttree%3A2.1.11--hec16e2b_1':
        'biocontainers/fasttree:2.1.11--hec16e2b_1' }"

    input:
    path msa

    output:
    path "asvs.msa.tree"        , emit: tree
    path "versions.yml"         , emit: versions

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
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fasttree: \$(fasttree -help 2>&1 | head -n 1 | sed -n 's/.*\\([0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\).*/\\1/p')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    
    """
}
