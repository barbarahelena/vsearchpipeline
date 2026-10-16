nextflow.enable.types = true

process VSEARCH_FASTQFILTER {
    tag meta.id
    label 'process_single_low'
    label 'vsearch'
    label 'error_retry'

    input:
    tuple(meta: Map, reads: Path)
    maxee: Float        // expected errors can be fractional, e.g. 0.5
    minlength: Integer  // 0 = no minimum
    maxlength: Integer  // 0 = no maximum
    maxns: Integer      // 0 = vsearch default

    output:
    reads        = tuple(meta, file("*.filtered.fastq"))
    filter_stats = tuple(meta, file("*.filter_stats.txt"))

    topic:
    tuple(task.process, 'vsearch', eval('vsearch --version 2>&1 | sed -nE \'1s/^vsearch v([0-9.]+).*/\\1/p\'')) >> 'versions'

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def min       = minlength != 0 ? "--fastq_minlen ${minlength}" : ""
    def max       = maxlength != 0 ? "--fastq_maxlen ${maxlength}" : ""
    def maxns_arg = maxns != 0 ? "--fastq_maxns ${maxns}" : ""
    """
    vsearch \\
        --fastq_filter ${reads} \\
        --fastq_maxee ${maxee} \\
        ${min} \\
        ${max} \\
        ${maxns_arg} \\
        ${args} \\
        --fastqout ${prefix}.filtered.fastq \\
        2> ${prefix}.filter_stats.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.filtered.fastq
    touch ${prefix}.filter_stats.txt
    """
}
