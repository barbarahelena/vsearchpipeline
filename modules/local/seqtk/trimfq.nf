process SEQTK_TRIMFQ {
    tag "$meta.id"
    label 'process_single_low'

    conda "bioconda::seqtk=1.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.4--he4a0461_1':
        'biocontainers/seqtk:1.4--h7132678_0' }"

    input:
    tuple val(meta), path(reads)
    val primers

    output:
    tuple val(meta), path("*.trim.fastq.gz")                                                       , emit: reads
    tuple val("${task.process}"), val('seqtk'), eval('echo $(seqtk --version 2>&1)'), emit: versions, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fwd_chars = primers.forward.size()
    def rev_chars = primers.reverse.size()
    def fwd_reads = reads[0]
    def rev_reads = reads[1]
    def fwd_trimmed = "${prefix}_1.trim.fastq.gz"
    def rev_trimmed = "${prefix}_2.trim.fastq.gz"

    """
    primer_fwd_rep=\$(echo "${primers.forward}" | sed 's/[^ACGT]/./g')
    primer_rev_rep=\$(echo "${primers.reverse}" | sed 's/[^ACGT]/./g')

    if [ "\$(zcat "$fwd_reads" | grep -c "\$primer_fwd_rep")" -ge 10 ] && [ "\$(zcat "$rev_reads" | grep -c "\$primer_rev_rep")" -ge 10 ] 
    then
        seqtk trimfq -b $fwd_chars "$fwd_reads" | gzip > $fwd_trimmed
        seqtk trimfq -b $rev_chars "$rev_reads" | gzip > $rev_trimmed
    else
        echo 'Warning: Few primer matches found, no trimming performed' >&2
        cp "$fwd_reads" "$fwd_trimmed"
        cp "$rev_reads" "$rev_trimmed"
    fi
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fwd_chars = primers.forward.size()
    def rev_chars = primers.reverse.size()
    def fwd_reads = reads[0]
    def rev_reads = reads[1]
    def fwd_trimmed = "${prefix}_1.trim.fastq.gz"
    def rev_trimmed = "${prefix}_2.trim.fastq.gz"

    """
    touch $fwd_trimmed
    touch $rev_trimmed
    """
}
