#!/usr/bin/env Rscript

library(phyloseq)
library(ggplot2)

## Open data
phylo <- readRDS("${phyloseq}")
numbers_before_rarefaction <- paste0(
    'The phyloseq object has ', nsamples(phylo), ' samples and ',
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
rarelevel <- mean(colSums(phylo@otu_table)) - 3 * sd(rowSums(phylo@otu_table))
if (rarelevel <= 15000) { rarelevel <- median(colSums(phylo@otu_table)) - IQR(colSums(phylo@otu_table)) }
if (rarelevel <= 15000) { rarelevel <- 15000 }
if (all(colSums(phylo@otu_table) < 15000)) { rarelevel <- min(colSums(phylo@otu_table)) }

print(paste0('Max counts: ',  max(colSums(phylo@otu_table))))
print(paste0('Min counts: ',  min(colSums(phylo@otu_table))))
print(paste0('Mean counts: ', mean(colSums(phylo@otu_table))))
print(paste0('SD counts: ',   sd(colSums(phylo@otu_table))))

## Rarefaction
if (sum(phylo@otu_table[1,]) != sum(phylo@otu_table[2,])) {
    rarefaction_yesno <- 'Rowsums are unequal, the data has not been rarefied yet.\n'
    print(rarefaction_yesno)
    user_rarelevel <- ${nf_rarelevel}
    min_counts     <- min(colSums(phylo@otu_table))
    if (user_rarelevel == 0) {
        rarefaction_outcome <- paste0('Rarefaction level: ', rarelevel, '\n')
    } else if (user_rarelevel > min_counts) {
        ## User-defined level exceeds the minimum sample depth — fall back to
        ## the auto-calculated level so we don't lose all samples (e.g. in tests).
        warning(paste0(
            'User-defined rarefaction level (', user_rarelevel, ') exceeds the ',
            'minimum sample count (', min_counts, '). ',
            'Falling back to the automatically calculated rarefaction level.'
        ))
        rarefaction_outcome <- paste0(
            'User-defined rarefaction level (', user_rarelevel,
            ') exceeded minimum sample count (', min_counts,
            '). Auto-calculated rarefaction level used: ', rarelevel, '\n'
        )
    } else {
        rarelevel <- user_rarelevel
        rarefaction_outcome <- paste0('Rarefaction level was user-defined: ', rarelevel, '\n')
    }
    print(rarefaction_outcome)

    phylo_rare <- tryCatch({
        rare <- rarefy_even_depth(phylo, sample.size = rarelevel,
                                  rngseed = ${nf_seed}, replace = FALSE,
                                  trimOTUs = TRUE, verbose = TRUE)
        if (nsamples(rare) == 0) {
            stop(paste0(
                'No samples remained after rarefaction at depth ', rarelevel,
                '. Minimum sample count is ', min_counts,
                '. Consider lowering the rarefaction level or setting rarelevel = 0 ',
                'to use the auto-calculated level.'
            ))
        }
        rare
    }, error = function(e) {
        stop(paste0(
            'Rarefaction failed: ', conditionMessage(e), '\n',
            'Minimum sample count: ', min_counts, '. ',
            'Rarefaction level attempted: ', rarelevel, '.\n',
            'Consider lowering the rarefaction level or setting it to 0 ',
            'to use the auto-calculated level.'
        ))
    })

    numbers_after_rarefaction <- paste0(
        'The phyloseq object has ', nsamples(phylo_rare), ' samples and ',
        ntaxa(phylo_rare), ' taxa. \n')
    print(numbers_after_rarefaction)
    report <- paste0(numbers_before_rarefaction,
                     rarefaction_yesno,
                     rarefaction_outcome,
                     numbers_after_rarefaction)
} else {
    rarefaction_yesno <- 'Phyloseq object seems already rarefied.\n'
    user_rarelevel <- ${nf_rarelevel}
    rarelevel <- ifelse(user_rarelevel != 0, user_rarelevel, 0)
    print(rarefaction_yesno)
    report <- paste0(numbers_before_rarefaction,
                     rarefaction_yesno)
    phylo_rare <- phylo
}

## Plot histogram
raredf <- data.frame(rarefaction = colSums(phylo@otu_table))
rarepl <- ggplot(raredf, aes(x = rarefaction)) +
    geom_histogram(color = "black", fill = "royalblue", alpha = 0.8) +
    geom_vline(aes(xintercept = rarelevel), color = "firebrick") +
    theme_minimal() +
    xlab("Number of counts") +
    ggtitle("Total counts distribution")
ggsave(plot = rarepl, 'rarefaction_plot.pdf')

writeLines(report, 'rarefaction_report.txt')
saveRDS(phylo_rare, 'phyloseq_rarefied.RDS')

################################################
## VERSIONS                                   ##
################################################

r.version       <- paste(R.version[['major']], R.version[['minor']], sep = '.')
phyloseq.version <- as.character(packageVersion('phyloseq'))
ggplot2.version  <- as.character(packageVersion('ggplot2'))

writeLines(
    c('"${task.process}":',
      paste0('    r-base: ',   r.version),
      paste0('    phyloseq: ', phyloseq.version),
      paste0('    ggplot2: ',  ggplot2.version)
    ),
    'versions.yml'
)
