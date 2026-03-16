#!/usr/bin/env Rscript

library(Biostrings)
library(phyloseq)

dna        <- readDNAStringSet("${asvs}")
counttable <- read.delim("${counttable}")
rownames(counttable) <- counttable[, 1]
counttable[, 1] <- NULL
taxtable <- read.csv("${taxtable}")
rownames(taxtable) <- taxtable[, 1]
taxtable[, 1] <- NULL

asv_seqs   <- as.character(dna)
asvs_names <- names(asv_seqs)
rownames(taxtable) <- asvs_names[match(rownames(taxtable), asv_seqs)]
taxtable   <- as.matrix(taxtable)
counttable <- as.matrix(counttable)

if (as.logical("${treepresent}")) {
    tree        <- ape::read.tree("${tree}")
    tree_rooted <- phytools::midpoint.root(tree)
    ps <- phyloseq(otu_table(counttable, taxa_are_rows = TRUE),
                   tax_table(taxtable), asv_seqs, tree_rooted)
} else {
    ps <- phyloseq(otu_table(counttable, taxa_are_rows = TRUE),
                   tax_table(taxtable), asv_seqs)
}
ps@refseq <- dna
ps
nsamples(ps)
ntaxa(ps)
sample_names(ps)

saveRDS(ps, "phyloseq.RDS")
write.csv(ps@tax_table, "phylo_raw_taxtable.csv")

################################################
## VERSIONS                                   ##
################################################

r.version          <- paste(R.version[['major']], R.version[['minor']], sep = '.')
phyloseq.version   <- as.character(packageVersion('phyloseq'))
biostrings.version <- as.character(packageVersion('Biostrings'))

writeLines(
    c('"${task.process}":',
      paste0('    r-base: ',     r.version),
      paste0('    phyloseq: ',   phyloseq.version),
      paste0('    Biostrings: ', biostrings.version)
    ),
    'versions.yml'
)
