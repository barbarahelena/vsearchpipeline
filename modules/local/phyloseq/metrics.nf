nextflow.enable.types = true

process PHYLOSEQ_METRICS {
    label 'process_multi_med'
    label 'phyloseq'

    input:
    phyloseq: Path
    taxtable: Path
    postfix: String   // 'complete' or 'rarefied', used in output file names

    output:
    speciescomp = file("composition_species_${postfix}.pdf")
    genuscomp   = file("composition_genus_${postfix}.pdf")
    famcomp     = file("composition_family_${postfix}.pdf")
    phylumcomp  = file("composition_phylum_${postfix}.pdf")
    shannon     = file("shannon_index_${postfix}.pdf")
    richness    = file("species_richness_${postfix}.pdf")
    metrics     = file("metrics_overview_${postfix}.txt")

    topic:
    tuple(task.process, 'r-base',   file("version_r.txt"))        >> 'versions'
    tuple(task.process, 'phyloseq', file("version_phyloseq.txt")) >> 'versions'
    tuple(task.process, 'vegan',    file("version_vegan.txt"))    >> 'versions'
    tuple(task.process, 'ggplot2',  file("version_ggplot2.txt"))  >> 'versions'

    script:
    template 'metrics.R'

    stub:
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
