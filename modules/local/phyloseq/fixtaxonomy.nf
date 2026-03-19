process PHYLOSEQ_FIXTAXONOMY {
    label 'process_multi_low'
    label 'phyloseq'

    input:
    path    phyloseq
    val     complete
    
    output:
    path "taxtable_*.RDS"           , emit: taxonomy
    path "phylogen_levels_*.csv"    , emit: phylevels
    path "phylogen_levels_top300_*.csv" , emit: phylevelstop
    tuple val("${task.process}"), val('r-base'),   path("version_r.txt"),        emit: versions_r,        topic: versions
    tuple val("${task.process}"), val('phyloseq'), path("version_phyloseq.txt"), emit: versions_phyloseq, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    postfix = complete ? "complete" : "rarefied"
    template 'fixtaxonomy.R'
    
    stub:
    postfix = complete ? "complete" : "rarefied"
    """
    touch taxtable_${postfix}.RDS
    touch phylogen_levels_${postfix}.csv
    touch phylogen_levels_top300_${postfix}.csv
    touch version_r.txt
    touch version_phyloseq.txt
    """
}
