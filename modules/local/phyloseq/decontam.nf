process PHYLOSEQ_DECONTAM {
    label 'process_multi_low'
    label 'phyloseq'

    input:
    path phyloseq
    path samplesheet

    output:
    path "phyloseq_decontam.RDS"        , emit: phyloseq
    path "decontam_report.txt"          , emit: report
    path "decontam_contaminants.csv"    , emit: contaminants
    path "decontam_prev_plot.pdf"       , emit: prev_plot
    tuple val("${task.process}"), val('r-base'),   path("version_r.txt"),        emit: versions_r,        topic: versions
    tuple val("${task.process}"), val('phyloseq'), path("version_phyloseq.txt"), emit: versions_phyloseq, topic: versions
    tuple val("${task.process}"), val('decontam'), path("version_decontam.txt"), emit: versions_decontam, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'decontam.R'

    stub:
    """
    touch phyloseq_decontam.RDS
    touch decontam_report.txt
    touch decontam_contaminants.csv
    touch decontam_prev_plot.pdf
    touch version_r.txt
    touch version_phyloseq.txt
    touch version_decontam.txt
    """
}
