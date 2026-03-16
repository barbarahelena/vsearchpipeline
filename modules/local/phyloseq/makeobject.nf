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
    path "versions.yml"             , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    treepresent = tree.name != 'NO_TREEFILE' ? "TRUE" : "FALSE"
    template 'makeobject.R'

    stub:
    """
    touch phyloseq.RDS
    touch phylo_raw_taxtable.csv
    touch versions.yml
    """
}
