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
    path "versions.yml"                 , emit: versions, topic: versions

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
    touch versions.yml
    """
}
