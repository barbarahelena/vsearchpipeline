nextflow.enable.types = true

process PHYLOSEQ_DECONTAM {
    label 'process_multi_low'
    label 'phyloseq'

    input:
    phyloseq: Path
    samplesheet: Path   // needs NEGCON-labelled samples and/or a nucl_acid_conc column

    output:
    phyloseq     = file("phyloseq_decontam.RDS")
    report       = file("decontam_report.txt")
    contaminants = file("decontam_contaminants.csv")
    prev_plot    = file("decontam_prev_plot.pdf")

    topic:
    tuple(task.process, 'r-base',   file("version_r.txt"))        >> 'versions'
    tuple(task.process, 'phyloseq', file("version_phyloseq.txt")) >> 'versions'
    tuple(task.process, 'decontam', file("version_decontam.txt")) >> 'versions'

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
