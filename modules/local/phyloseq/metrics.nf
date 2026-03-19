process PHYLOSEQ_METRICS {
    label 'process_multi_med'
    label 'phyloseq'

    input:
    path    phyloseq
    path    taxtable
    val     complete

    output:
    path "composition_species_*.pdf"    , emit: speciescomp
    path "composition_genus_*.pdf"      , emit: genuscomp
    path "composition_family_*.pdf"     , emit: famcomp
    path "composition_phylum_*.pdf"     , emit: phylumcomp
    path "shannon_index_*.pdf"          , emit: shannon
    path "species_richness_*.pdf"       , emit: richness
    path "metrics_overview_*.txt"       , emit: metrics
    tuple val("${task.process}"), val('r-base'),   path("version_r.txt"),        emit: versions_r,        topic: versions
    tuple val("${task.process}"), val('phyloseq'), path("version_phyloseq.txt"), emit: versions_phyloseq, topic: versions
    tuple val("${task.process}"), val('vegan'),    path("version_vegan.txt"),    emit: versions_vegan,    topic: versions
    tuple val("${task.process}"), val('ggplot2'),  path("version_ggplot2.txt"),  emit: versions_ggplot2,  topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    postfix = complete ? "complete" : "rarefied"
    template 'metrics.R'

    stub:
    postfix = complete ? "complete" : "rarefied"
    """
    touch composition_species_${postfix}.pdf
    touch composition_genus_${postfix}.pdf
    touch composition_family_${postfix}.pdf
    touch composition_phylum_${postfix}.pdf
    touch shannon_index_${postfix}.pdf
    touch species_richness_${postfix}.pdf
    touch metrics_overview_${postfix}.txt
    touch version_r.txt
    touch version_phyloseq.txt
    touch version_vegan.txt
    touch version_ggplot2.txt
    """
}
