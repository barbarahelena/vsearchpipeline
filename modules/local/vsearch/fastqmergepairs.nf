process VSEARCH_FASTQMERGEPAIRS {
    tag "$meta.id"
    label 'process_multi_low'
    label 'vsearch'
    label 'error_retry'

    input:
    tuple val(meta), path(reads)
    val allowmergestagger
    val maxdiffs
    val minlength
    val maxlength
    val maxdiffpct

    output:
    tuple val(meta), path("*.merged.fastq.gz")                                                                            , emit: reads
    tuple val("${task.process}"), val('vsearch'), eval('vsearch --version 2>&1 | head -n 1 | sed \'s/vsearch //; s/,.*//\''), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fwd_reads = reads[0]
    def rev_reads = reads[1]
    def merged = "${prefix}.merged.fastq.gz"
    def allowmergestagger = allowmergestagger ? "--fastq_allowmergestagger" : ''
    def maxdiffs = maxdiffs ? "--fastq_maxdiffs ${maxdiffs}" : ''
    def minmerge = minlength != 0 ? "--fastq_minmergelen ${minlength}" : ''
    def maxmerge = maxlength != 0 ? "--fastq_maxmergelen ${maxlength}" : ''
    def maxdiffpct = maxdiffpct != 100 ? "--fastq_maxdiffpct ${maxdiffpct}" : ''

    """
    vsearch \\
        --fastq_mergepairs ${fwd_reads} \\
        --reverse ${rev_reads} \\
        ${maxdiffs} \\
        ${allowmergestagger} \\
        ${minmerge} \\
        ${maxmerge} \\
        ${maxdiffpct} \\
        --fastqout ${merged}
    """

    stub:
    def stub_prefix = task.ext.prefix ?: "${meta.id}"
    def stub_merged = "${stub_prefix}.merged.fastq.gz"
    """
    touch ${stub_merged}
    """
}