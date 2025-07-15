process PICRUST2 {
    label 'process_multi_high'
    label 'picrust2'
    label 'error_retry'

    input:
    path(asvfasta)
    path(asvtab)

    output:
    path("picrust_output/*") , emit: outfolder
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    picrust2_pipeline.py $args -s $asvfasta -i $asvtab -o picrust_output -p ${task.cpus} --in_traits EC,KO --verbose

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        picrust2: \$( picrust2_pipeline.py -v | sed -e "s/picrust2_pipeline.py //g" )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    """
}
