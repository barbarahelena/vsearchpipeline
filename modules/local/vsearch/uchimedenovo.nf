nextflow.enable.types = true

process VSEARCH_UCHIMEDENOVO {
    label 'process_single_low'
    label 'vsearch'

    input:
    reads: Path
    label: String   // prefix for ASV names, e.g. 'ASV_'

    output:
    chimeras = file("chimeras.fasta")
    asvs     = file("asvs_nonchimeras.fasta")

    topic:
    tuple(task.process, 'vsearch', eval('vsearch --version 2>&1 | sed -nE \'1s/^vsearch v([0-9.]+).*/\\1/p\'')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    """
    vsearch \\
        --uchime3_denovo ${reads} \\
        --chimeras chimeras.fasta \\
        --nonchimeras asvs_nonchimeras.fasta \\
        --threads ${task.cpus} \\
        ${args} \\
        --relabel ${label}
    """

    stub:
    """
    touch chimeras.fasta
    touch asvs_nonchimeras.fasta
    """
}
