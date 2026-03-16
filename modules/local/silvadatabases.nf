process SILVADATABASES {
    label 'process_single_med'

    output:
    path "SILVA_asv_db.fa.gz"                                        , emit: asvdb
    path "SILVA_species_db.fa.gz"                                    , emit: speciesdb
    tuple val("${task.process}"), val('SILVA'), eval('echo "138.2"') , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    wget -c https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toGenus_trainset.fa.gz?download=1 \\
        -O SILVA_asv_db.fa.gz
    wget -c https://zenodo.org/records/14169026/files/silva_v138.2_assignSpecies.fa.gz?download=1 \\
        -O SILVA_species_db.fa.gz
    """

    stub:
  
    """
    touch SILVA_asv_db.fa.gz
    touch SILVA_species_db.fa.gz
    """
}
