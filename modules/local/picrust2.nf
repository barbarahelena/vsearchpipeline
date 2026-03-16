process PICRUST2 {
    label 'process_multi_high'
    label 'picrust2'
    label 'error_retry'

    input:
    path(asvfasta)
    path(asvtab)

    output:
    path("picrust_output/*")                                                                                                   , emit: outfolder
    tuple val("${task.process}"), val('picrust2'), eval('picrust2_pipeline.py -v | sed "s/picrust2_pipeline.py //"'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    picrust2_pipeline.py $args -s $asvfasta -i $asvtab -o picrust_output -p ${task.cpus} --in_traits EC,KO --verbose
    """

    stub:
    def args = task.ext.args ?: ''
    """
    """
}
