nextflow.enable.types = true

process VSEARCH_SORT_REMOVE_SINGLETONS {
    label 'process_single_low'
    label 'vsearch'
    label 'error_retry'

    input:
    reads: Path
    fasta_width: Integer   // 0 = no line wrapping
    minsize: Integer

    output:
    asvs = file("asvs_nonsingle.fasta")

    topic:
    tuple(task.process, 'vsearch', eval('vsearch --version 2>&1 | sed -nE \'1s/^vsearch v([0-9.]+).*/\\1/p\'')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    """
    vsearch \\
        --sortbysize ${reads} \\
        --threads ${task.cpus} \\
        --sizein \\
        --sizeout \\
        --fasta_width ${fasta_width} \\
        --minsize ${minsize} \\
        ${args} \\
        --output asvs_nonsingle.fasta
    """

    stub:
    """
    touch asvs_nonsingle.fasta
    """
}
