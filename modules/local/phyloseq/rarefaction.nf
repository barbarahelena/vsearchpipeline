process PHYLOSEQ_RAREFACTION {
    label 'process_multi_low'
    label 'phyloseq'

    input:
    path    phyloseq
    val     rarelevel

    output:
    path "phyloseq_rarefied.RDS"    , emit: phyloseq
    path "rarefaction_plot.pdf"     , emit: rarecurve
    path "rarefaction_report.txt"   , emit: rarereport

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def seed = task.ext.seed ?: '1234'
    def rarelevel = rarelevel ? "${rarelevel}" : 0

    """
    #!/usr/bin/env Rscript
    library(phyloseq)
    library(ggplot2)

    ## Open data
    phylo <- readRDS('$phyloseq')
    numbers_before_rarefaction <- paste0('The phyloseq object has ', nsamples(phylo), ' samples and ',
                 ntaxa(phylo), ' taxa. \n')
    print(numbers_before_rarefaction)

    ## Remove any remaining negative/blank control samples before rarefaction
    ctrl_samples <- sample_names(phylo)[grepl("NEGCON|^SB", sample_names(phylo))]
    if (length(ctrl_samples) > 0) {
        cat(paste0("Removing ", length(ctrl_samples), " control sample(s) before rarefaction: ",
                   paste(ctrl_samples, collapse = ", "), "\n"))
        phylo <- prune_samples(!sample_names(phylo) %in% ctrl_samples, phylo)
        cat(paste0("After control removal: ", nsamples(phylo), " samples.\n"))
    }
    
    ## Calculated rarefaction level
    rarelevel <- mean(colSums(phylo@otu_table)) - 3*sd(rowSums(phylo@otu_table))
    if(rarelevel <= 15000){ rarelevel <- median(colSums(phylo@otu_table)) - IQR(colSums(phylo@otu_table)) }
    if(rarelevel <= 15000){ rarelevel <- 15000}
    if(all(colSums(phylo@otu_table) < 15000)){ rarelevel <- min(colSums(phylo@otu_table))}

    print(paste0('Max counts: ', max(colSums(phylo@otu_table))))
    print(paste0('Min counts: ', min(colSums(phylo@otu_table))))
    print(paste0('Mean counts: ', mean(colSums(phylo@otu_table))))
    print(paste0('SD counts: ', sd(colSums(phylo@otu_table))))

    ## Rarefaction
    if(sum(phylo@otu_table[1,]) != sum(phylo@otu_table[2,])){
        rarefaction_yesno <- 'Rowsums are unequal, the data has not been rarefied yet.\n'
        print(rarefaction_yesno)
        if($rarelevel == 0){
            rarefaction_outcome <- paste0('Rarefaction level: ', rarelevel, '\n')
        } else{
            rarelevel <- $rarelevel
            rarefaction_outcome <- paste0('Rarefaction level was user-defined: ', rarelevel, '\n')
        }
        print(rarefaction_outcome)
    
        phylo_rare <- rarefy_even_depth(phylo, sample.size = rarelevel, 
                                        rngseed = $seed, replace = FALSE, 
                                        trimOTUs = TRUE, verbose = TRUE)
        
        numbers_after_rarefaction <- paste0('The phyloseq object has ', nsamples(phylo_rare), ' samples and ',
                    ntaxa(phylo_rare), ' taxa. \n')
        print(numbers_after_rarefaction)
        report <- paste0(numbers_before_rarefaction, 
                            rarefaction_yesno, 
                            rarefaction_outcome,
                            numbers_after_rarefaction)
    } else{
        rarefaction_yesno <- 'Phyloseq object seems already rarefied.\n'
        rarelevel <- ifelse($rarelevel !=0, $rarelevel, 0)
        print(rarefaction_yesno)
        report <- paste0(numbers_before_rarefaction, 
                            rarefaction_yesno)
        phylo_rare <- phylo
    }

    ## Plot histogram
    raredf <- data.frame(rarefaction = colSums(phylo@otu_table))
    rarepl <- ggplot(raredf, aes(x = rarefaction)) +
                geom_histogram(color = "black", 
                            fill = "royalblue", alpha = 0.8) +
                geom_vline(aes(xintercept = rarelevel), color = "firebrick") + 
                theme_minimal() +
                xlab("Number of counts") +
                ggtitle("Total counts distribution")
    ggsave(plot = rarepl, 'rarefaction_plot.pdf')
    
    writeLines(report, 'rarefaction_report.txt')
    saveRDS(phylo_rare, 'phyloseq_rarefied.RDS')
    """

    stub:
    def args = task.ext.args ?: ''
    def seed = task.ext.seed ?: '1234'
    def rarelevel = rarelevel ? "${rarelevel}" : 0

    """
    touch phyloseq_rarefied.RDS
    touch rarehist.pdf
    touch rarereport.txt
    """
}
