process VSEARCH_FASTQFILTER {
    tag "$meta.id"
    label 'process_single_low'
    label 'vsearch'
    label 'error_retry'

    input:
    tuple val(meta), path(reads)
    val maxee
    val minlength
    val maxlength
    val maxns

    output:
    tuple val(meta), path("*.filtered.fastq")                                                                         , emit: reads
    tuple val(meta), path("*.filter_stats.txt")                                                                       , emit: filter_stats
    tuple val("${task.process}"), val('vsearch'), eval('vsearch --version 2>&1 | head -n 1 | sed \'s/vsearch //; s/,.*//\''), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def filtered = "${prefix}.filtered.fastq"
    def min = minlength != 0 ? "--fastq_minlen ${minlength}" : ""
    def max = maxlength != 0 ? "--fastq_maxlen ${maxlength}" : ""
    def maxns = maxns != 0 ? "--fastq_maxns ${maxns}" : ""

    """
    vsearch \\
        --fastq_filter $reads \\
        --fastq_maxee $maxee \\
        $min \\
        $max \\
        $maxns \\
        --fastqout $filtered \\
        2> ${prefix}.filter_stats.txt
    """

    stub:
    def stub_prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${stub_prefix}.filtered.fastq"
    touch "${stub_prefix}.filter_stats.txt"
    """
}
