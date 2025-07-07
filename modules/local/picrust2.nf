process PICRUST2 {
    label 'process_high'
    label 'picrust2'

    input:
    path(asvfasta)
    path(asvtab)

    output:
    path("picrust_output/*") , emit: outfolder
    path("*_descrip.tsv"), emit: pathways
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    picrust2_pipeline.py $args -s $asvfasta -i $asvtab -o picrust_output -p ${task.cpus} --in_traits EC,KO --verbose
    add_descriptions.py -i picrust_output/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC -o EC_pred_metagenome_unstrat_descrip.tsv
    add_descriptions.py -i picrust_output/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO -o KO_pred_metagenome_unstrat_descrip.tsv
    add_descriptions.py -i picrust_output/pathways_out/path_abun_unstrat.tsv.gz -m METACYC -o METACYC_path_abun_unstrat_descrip.tsv

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
