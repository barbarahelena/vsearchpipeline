nextflow.enable.types = true

process VSEARCH_CLUSTERUNOISE {
    label 'process_multi_verylow'
    label 'error_retry'
    label 'vsearch'

    input:
    reads: Path
    minsize: Integer
    alpha: Number

    output:
    asvs = file("asvs.clustered.fasta")

    topic:
    tuple(task.process, 'vsearch', eval('vsearch --version 2>&1 | sed -nE \'1s/^vsearch v([0-9.]+).*/\\1/p\'')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    """
    vsearch \\
        --cluster_unoise ${reads} \\
        --centroids asvs.clustered.fasta \\
        --minsize ${minsize} \\
        --threads ${task.cpus} \\
        ${args} \\
        --unoise_alpha ${alpha}
    """

    stub:
    """
    touch asvs.clustered.fasta
    """
}
