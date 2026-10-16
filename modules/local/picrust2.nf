nextflow.enable.types = true

process PICRUST2 {
    label 'process_multi_high'
    label 'picrust2'
    label 'error_retry'

    input:
    asvfasta: Path
    asvtab: Path

    output:
    outfolder = files("picrust_output/*")

    topic:
    tuple(task.process, 'picrust2', eval('picrust2_pipeline.py -v 2>&1 | sed -nE \'s/^picrust2_pipeline.py ([0-9.]+).*/\\1/p\'')) >> 'versions'

    script:
    def args = task.ext.args ?: ''
    """
    picrust2_pipeline.py ${args} -s ${asvfasta} -i ${asvtab} -o picrust_output -p ${task.cpus} --in_traits EC,KO --verbose
    """

    stub:
    """
    mkdir -p picrust_output/EC_metagenome_out picrust_output/KO_metagenome_out
    touch picrust_output/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz
    touch picrust_output/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz
    """
}
