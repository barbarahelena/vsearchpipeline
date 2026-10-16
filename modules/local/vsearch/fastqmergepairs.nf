nextflow.enable.types = true

process VSEARCH_FASTQMERGEPAIRS {
    tag meta.id
    label 'process_multi_low'
    label 'vsearch'
    label 'error_retry'

    input:
    tuple(meta: Map, reads: List<Path>)   // [ forward, reverse ]
    allowmergestagger: Boolean
    maxdiffs: Integer                     // 0 = vsearch default
    minlength: Integer                    // 0 = no minimum
    maxlength: Integer                    // 0 = no maximum
    maxdiffpct: Integer                   // 100 = vsearch default

    output:
    reads = tuple(meta, file("*.merged.fastq.gz"))

    topic:
    tuple(task.process, 'vsearch', eval('vsearch --version 2>&1 | sed -nE \'1s/^vsearch v([0-9.]+).*/\\1/p\'')) >> 'versions'

    script:
    def args                  = task.ext.args ?: ''
    def prefix                = task.ext.prefix ?: "${meta.id}"
    def allowmergestagger_arg = allowmergestagger ? "--fastq_allowmergestagger" : ''
    def maxdiffs_arg          = maxdiffs ? "--fastq_maxdiffs ${maxdiffs}" : ''
    def minmerge              = minlength != 0 ? "--fastq_minmergelen ${minlength}" : ''
    def maxmerge              = maxlength != 0 ? "--fastq_maxmergelen ${maxlength}" : ''
    def maxdiffpct_arg        = maxdiffpct != 100 ? "--fastq_maxdiffpct ${maxdiffpct}" : ''
    """
    vsearch \\
        --fastq_mergepairs ${reads[0]} \\
        --reverse ${reads[1]} \\
        ${maxdiffs_arg} \\
        ${allowmergestagger_arg} \\
        ${minmerge} \\
        ${maxmerge} \\
        ${maxdiffpct_arg} \\
        ${args} \\
        --fastqout ${prefix}.merged.fastq

    # --fastqout writes plain text whatever the extension, so compress explicitly
    gzip ${prefix}.merged.fastq
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged.fastq.gz
    """
}
