process SILVADATABASES {
    label 'process_single_med'
    storeDir 'db'

    output:
    path "SILVA_asv_db.fa.gz"             , emit: asvdb
    path "SILVA_species_db.fa.gz"         , emit: speciesdb
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    wget -c https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toGenus_trainset.fa.gz?download=1 \\
        -O SILVA_asv_db.fa.gz
    wget -c https://zenodo.org/records/14169026/files/silva_v138.2_assignSpecies.fa.gz?download=1 \\
        -O SILVA_species_db.fa.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SILVA: 138.2
    END_VERSIONS
    """

    stub:
  
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SILVA: 138.2
    END_VERSIONS
    """
}
