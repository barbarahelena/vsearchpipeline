nextflow.enable.types = true

process VSEARCH_DEREPFULLLENGTH {
    tag id
    label 'process_single_low'
    label 'vsearch'
    label 'error_retry'

    input:
    record(id: String, filtered: Path)
    strand: String   // 'plus' or 'both'

    output:
    sample = record(id: id, derep: file("*.derep.fasta"))

    topic:
    tuple(task.process, 'vsearch', eval('vsearch --version 2>&1 | sed -nE \'1s/^vsearch v([0-9.]+).*/\\1/p\'')) >> 'versions'

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: id
    """
    vsearch \\
        --fastx_uniques ${filtered} \\
        --fastaout ${prefix}.derep.fasta \\
        --strand ${strand} \\
        --sizeout \\
        ${args} \\
        --relabel ${prefix}.
    """

    stub:
    def prefix = task.ext.prefix ?: id
    """
    touch ${prefix}.derep.fasta
    """
}
