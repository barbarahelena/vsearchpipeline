process PRIMERSHEET_CHECK {
    tag "$primersheet"
    label 'process_multi_low'

    conda "conda-forge::pandas=2.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    path primersheet

    output:
    path '*.csv'       , emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in barbarahelena/vsearchpipeline/modules/local/templates/
    template 'check_primersheet.py'
}
