nextflow.enable.types = true

process VSEARCH_DEREPFULLLENGTHALL {
    label 'process_single_med'
    label 'vsearch'
    label 'error_retry'

    input:
    reads: List<Path>     // per-sample dereplicated fasta files
    strand: String        // 'plus' or 'both'
    fastawidth: Integer   // 0 = no line wrapping
    minunique: Integer

    output:
    concatreads = file("all.concat.fasta")
    reads       = file("all.derep.fasta")

    topic:
    tuple(task.process, 'vsearch', eval('vsearch --version 2>&1 | sed -nE \'1s/^vsearch v([0-9.]+).*/\\1/p\'')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    """
    cat ${reads.join(' ')} > all.concat.fasta

    vsearch \\
        --fastx_uniques all.concat.fasta \\
        --fastaout all.derep.fasta \\
        --strand ${strand} \\
        --sizein \\
        --sizeout \\
        --fasta_width ${fastawidth} \\
        ${args} \\
        --minuniquesize ${minunique}
    """

    stub:
    """
    touch all.concat.fasta
    touch all.derep.fasta
    """
}
