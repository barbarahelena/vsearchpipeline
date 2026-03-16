process MAFFT {
    label 'process_single_low'
    label 'mafft'

    input:
    path asvs

    output:
    path "asvs.msa"                                                                                    , emit: msa
    tuple val("${task.process}"), val('mafft'), eval("mafft --version 2>&1 | sed 's/^v//; s/ (.*//'" ), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mafft --auto $asvs > asvs.msa
    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    touch asvs.msa
    """
}
