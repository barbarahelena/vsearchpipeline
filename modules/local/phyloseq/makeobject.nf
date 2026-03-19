process PHYLOSEQ_MAKEOBJECT {
    label 'process_multi_low'
    label 'phyloseq'

    input:
    path asvs
    path counttable
    path tree
    path taxtable

    output:
    path "phyloseq.RDS"             , emit: phyloseq
    path "phylo_raw_taxtable.csv"   , emit: taxtable
    tuple val("${task.process}"), val('r-base'),    path("version_r.txt"),          emit: versions_r,          topic: versions
    tuple val("${task.process}"), val('phyloseq'),  path("version_phyloseq.txt"),   emit: versions_phyloseq,   topic: versions
    tuple val("${task.process}"), val('Biostrings'),path("version_biostrings.txt"), emit: versions_biostrings, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    treepresent = tree.name != 'NO_TREEFILE' ? "TRUE" : "FALSE"
    template 'makeobject.R'

    stub:
    """
    touch phyloseq.RDS
    touch phylo_raw_taxtable.csv
    touch version_r.txt
    touch version_phyloseq.txt
    touch version_biostrings.txt
    """
}
