process VSEARCH_DEREPFULLLENGTHALL {
    //tag "$meta.id"
    label 'process_single_med'
    label 'vsearch'
    label 'error_retry'

    input:
    path(reads)
    val strand
    val fastawidth
    val minunique

    output:
    path("all.concat.fasta") , emit: concatreads
    path("all.derep.fasta")  , emit: reads
    path "versions.yml"      , emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    cat $reads > all.concat.fasta
    
    vsearch \\
        --fastx_uniques all.concat.fasta \\
        --fastaout all.derep.fasta \\
        --strand $strand \\
        --sizein \\
        --sizeout \\
        --fasta_width $fastawidth \\
        --minuniquesize $minunique

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n 1 | sed 's/vsearch //; s/,.*//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch all.concat.fasta
    touch all.derep.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n 1 | sed 's/vsearch //; s/,.*//')
    END_VERSIONS
    """
}
