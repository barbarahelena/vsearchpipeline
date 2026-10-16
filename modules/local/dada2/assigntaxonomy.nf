nextflow.enable.types = true

process DADA2_ASSIGNTAXONOMY {
    label 'process_multi_high'
    label 'dada2'

    input:
    asvs: Path
    silva_asv_db: Path
    silva_species_db: Path
    minboot: Integer        // minimum bootstrap confidence, 0-100
    allowmultiple: Integer  // max number of species assignments per ASV
    tryrevcompl: Boolean

    output:
    taxtable = file("taxtable.csv")

    topic:
    tuple(task.process, 'R',     eval('Rscript -e "cat(paste(R.version[c(\'major\',\'minor\')], collapse=\'.\'))"')) >> 'versions'
    tuple(task.process, 'dada2', eval('Rscript -e "cat(as.character(packageVersion(\'dada2\')))"')) >> 'versions'

    script:
    def tryrc = tryrevcompl ? 'TRUE' : 'FALSE'
    def seed  = task.ext.seed ?: '1234'
    """
    Rscript - <<EOF
    suppressPackageStartupMessages(library(dada2))
    set.seed(${seed})
    taxtable <- assignTaxonomy( "${asvs}",
                                "${silva_asv_db}",
                                multithread = TRUE,
                                minBoot = ${minboot},
                                verbose = TRUE)
    taxa <- addSpecies( taxtable,
                        "${silva_species_db}",
                        verbose = TRUE,
                        allowMultiple = ${allowmultiple},
                        tryRC = ${tryrc})
    write.csv(taxa, file = "taxtable.csv", quote=FALSE)
    EOF
    """

    stub:
    """
    touch taxtable.csv
    """
}
