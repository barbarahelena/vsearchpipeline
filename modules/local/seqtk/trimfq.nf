nextflow.enable.types = true

process SEQTK_TRIMFQ {
    tag id
    label 'process_single_low'
    label 'seqtk'

    input:
    record(id: String, fastq: List<Path>)   // fastq: [ forward, reverse ]
    primers: Map                          // [ forward: 'ACGT...', reverse: 'ACGT...' ]

    output:
    // same shape as the input, so the merge step does not depend on whether trimming ran.
    // An explicit list keeps forward/reverse order; files() would return an unordered set.
    sample = record(id: id, fastq: [ file("*_1.trim.fastq.gz"), file("*_2.trim.fastq.gz") ])

    topic:
    tuple(task.process, 'seqtk', eval('seqtk 2>&1 | sed -nE \'s/^Version: ([^ ]+).*/\\1/p\'')) >> 'versions'

    script:
    def prefix      = task.ext.prefix ?: id
    def fwd_chars   = primers.forward.size()
    def rev_chars   = primers.reverse.size()
    def fwd_reads   = fastq[0]
    def rev_reads   = fastq[1]
    def fwd_trimmed = "${prefix}_1.trim.fastq.gz"
    def rev_trimmed = "${prefix}_2.trim.fastq.gz"
    """
    primer_fwd_rep=\$(echo "${primers.forward}" | sed 's/[^ACGT]/./g')
    primer_rev_rep=\$(echo "${primers.reverse}" | sed 's/[^ACGT]/./g')

    # gzip -dc, not zcat: macOS zcat only reads .Z files, which silently skipped trimming on Macs
    if [ "\$(gzip -dc "${fwd_reads}" | grep -c "\$primer_fwd_rep")" -ge 10 ] && [ "\$(gzip -dc "${rev_reads}" | grep -c "\$primer_rev_rep")" -ge 10 ]
    then
        seqtk trimfq -b ${fwd_chars} "${fwd_reads}" | gzip > ${fwd_trimmed}
        seqtk trimfq -b ${rev_chars} "${rev_reads}" | gzip > ${rev_trimmed}
    else
        echo 'Warning: Few primer matches found, no trimming performed' >&2
        cp "${fwd_reads}" "${fwd_trimmed}"
        cp "${rev_reads}" "${rev_trimmed}"
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: id
    """
    touch ${prefix}_1.trim.fastq.gz
    touch ${prefix}_2.trim.fastq.gz
    """
}
