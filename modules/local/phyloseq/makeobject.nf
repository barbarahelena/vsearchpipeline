nextflow.enable.types = true

process PHYLOSEQ_MAKEOBJECT {
    label 'process_multi_low'
    label 'phyloseq'

    input:
    asvs: Path
    counttable: Path
    tree: Path       // assets/NO_TREEFILE when --skip_tree
    taxtable: Path

    output:
    phyloseq = file("phyloseq.RDS")
    taxtable = file("phylo_raw_taxtable.csv")

    topic:
    // template scripts are R, not Bash, so they write version files instead of using eval()
    tuple(task.process, 'r-base',     file("version_r.txt"))          >> 'versions'
    tuple(task.process, 'phyloseq',   file("version_phyloseq.txt"))   >> 'versions'
    tuple(task.process, 'Biostrings', file("version_biostrings.txt")) >> 'versions'

    script:
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
