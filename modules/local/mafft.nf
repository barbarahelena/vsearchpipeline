nextflow.enable.types = true

process MAFFT {
    label 'process_single_low'
    label 'mafft'

    input:
    asvs: Path

    output:
    msa = file("asvs.msa")

    topic:
    tuple(task.process, 'mafft', eval('mafft --version 2>&1 | sed -nE \'s/^v([0-9.]+).*/\\1/p\'')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    """
    mafft --auto ${args} ${asvs} > asvs.msa
    """

    stub:
    """
    touch asvs.msa
    """
}
