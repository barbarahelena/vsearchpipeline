nextflow.enable.types = true

process SILVADATABASES {
    label 'process_single_med'

    output:
    asvdb     = file("SILVA_asv_db.fa.gz")
    speciesdb = file("SILVA_species_db.fa.gz")

    topic:
    tuple(task.process, 'SILVA', '138.2') >> 'versions'

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
