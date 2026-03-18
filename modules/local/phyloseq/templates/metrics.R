#!/usr/bin/env Rscript

library(phyloseq)
library(dplyr)
library(ggplot2)
library(tidyr)
library(vegan)
library(forcats)

cols <- c("darkgreen", 'firebrick', "navy", "dodgerblue", "goldenrod2", "chartreuse4",
          "darkorange2", "rosybrown1", "darkred", "lightskyblue",
          "seagreen", "gold1", "olivedrab", "royalblue", "linen", "maroon4",
          "mediumturquoise", "plum2", "darkslateblue", "sienna", "grey70", "grey90")

theme_composition <- function(base_size = 14, base_family = "sans") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size = base_size, base_family = base_family)
        + theme(
            plot.title   = element_text(face = "bold", size = rel(1.0), hjust = 0.5),
            text         = element_text(),
            panel.background = element_rect(colour = NA, fill = NA),
            plot.background  = element_rect(colour = NA, fill = NA),
            panel.border     = element_rect(colour = NA),
            axis.title       = element_text(face = "bold", size = rel(1)),
            axis.title.y     = element_text(angle = 90, vjust = 2),
            axis.title.x     = element_text(vjust = -0.2),
            axis.text.x      = element_text(angle = 45, hjust = 1),
            axis.text        = element_text(),
            axis.line        = element_line(colour = "black"),
            axis.ticks       = element_line(),
            panel.grid.major = element_line(colour = "#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key       = element_rect(colour = NA),
            legend.position  = "right",
            legend.key.size  = unit(0.4, "cm"),
            legend.spacing   = unit(0, "cm"),
            legend.text      = element_text(size = rel(0.7)),
            plot.margin      = unit(c(10, 5, 5, 5), "mm"),
            strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
            strip.text       = element_text(face = "bold")
        ))
}

theme_Publication <- function(base_size = 14, base_family = "sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    (theme_foundation(base_size = base_size, base_family = base_family)
        + theme(
            plot.title   = element_text(face = "bold", size = rel(0.8), hjust = 0.5),
            text         = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background  = element_rect(colour = NA),
            panel.border     = element_rect(colour = NA),
            axis.title       = element_text(face = "bold", size = rel(0.8)),
            axis.title.y     = element_text(angle = 90, vjust = 2),
            axis.title.x     = element_text(vjust = -0.2),
            axis.text        = element_text(),
            axis.line        = element_line(colour = "black"),
            axis.ticks       = element_line(),
            panel.grid.major = element_line(colour = "#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key       = element_rect(colour = NA),
            legend.position  = "bottom",
            legend.key.size  = unit(0.2, "cm"),
            legend.spacing   = unit(0, "cm"),
            plot.margin      = unit(c(10, 5, 5, 5), "mm"),
            strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
            strip.text       = element_text(face = "bold")
        ))
}

## Open data
phy     <- readRDS("${phyloseq}")
tax     <- readRDS("${taxtable}")
postfix <- "${postfix}"

tab <- as(phy@otu_table, 'matrix')
counts <- sample_sums(phy@otu_table)
tab <- as.data.frame(t(tab) / counts * 100)
tab_comp <- tab
tab_comp\$Sample <- rownames(tab_comp)
rowSums(tab)   # samples should all sum up to 100%

#### Species-level ####
N <- 20
d <- tab_comp %>%
    pivot_longer(-Sample, names_to = 'ASV', values_to = 'Abundance') %>%
    mutate(sampleID = Sample)
d\$Tax     <- tax\$Tax[match(d\$ASV, tax\$ASV)]
d\$Species <- tax\$Species[match(d\$ASV, tax\$ASV)]

top_taxa <- d %>%
    filter(!is.na(Species)) %>%
    group_by(ASV, Sample) %>%
    summarise(Abundance = sum(Abundance)) %>%
    group_by(ASV) %>%
    summarise(Abundance = mean(Abundance)) %>%
    arrange(-Abundance) %>%
    dplyr::select(ASV) %>%
    head(N)
toptax <- d\$Tax[match(top_taxa\$ASV, d\$ASV)]

dx <- d %>%
    mutate(
        Species2 = case_when(
            is.na(Species)                             ~ "Unknown",
            Tax %in% toptax                            ~ paste(Tax),
            !(Tax %in% toptax) & !(is.na(Species))    ~ "Other species"
        ),
        Species2 = as.factor(Species2)
    ) %>%
    group_by(Species2, Sample) %>%
    summarise(Abundance = sum(Abundance)) %>%
    group_by(Species2) %>%
    summarise(Abundance = mean(Abundance)) %>%
    mutate(
        group    = "all samples",
        Species2 = fct_reorder(Species2, Abundance),
        Species2 = fct_relevel(Species2, "Other species", after = 0L),
        Species2 = fct_relevel(Species2, "Unknown", after = 0L)
    )

lev   <- levels(dx\$Species2)
colsp <- c(cols[1:(length(lev) - 2)], cols[21:22])

comp_species <- dx %>%
    ggplot(aes(x = group, y = Abundance, fill = Species2)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = rev(colsp), labels = lev) +
    guides(fill = guide_legend(ncol = 1)) +
    labs(y = "Composition (%)", x = "", title = "Species", fill = "") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_composition()
ggsave(plot = comp_species, paste0("composition_species_", postfix, ".pdf"))

#### Genus-level ####
N <- 20
d <- tab_comp %>%
    pivot_longer(-Sample, names_to = 'ASV', values_to = 'Abundance') %>%
    mutate(sampleID = Sample)
d\$Genus <- tax\$Genus[match(d\$ASV, tax\$ASV)]

top_taxa <- d %>%
    group_by(Genus, Sample) %>%
    summarise(Abundance = sum(Abundance)) %>%
    group_by(Genus) %>%
    summarise(Abundance = mean(Abundance)) %>%
    arrange(-Abundance) %>%
    dplyr::select(Genus) %>%
    filter(Genus != 'ambiguous') %>%
    head(N) %>%
    unlist()
top_taxa

dx <- d %>%
    mutate(
        Genus2 = case_when(
            Genus %in% top_taxa            ~ paste(Genus),
            is.na(Genus)                   ~ paste("Unknown"),
            !(Genus %in% top_taxa)         ~ paste("Other genera")
        ),
        Genus2 = as.factor(Genus2)
    ) %>%
    group_by(Genus2, Sample) %>%
    summarise(Abundance = sum(Abundance)) %>%
    group_by(Genus2) %>%
    summarise(Abundance = mean(Abundance)) %>%
    mutate(
        group  = "all samples",
        Genus2 = fct_reorder(Genus2, Abundance),
        Genus2 = fct_relevel(Genus2, "Other genera", after = 0L),
        Genus2 = fct_relevel(Genus2, "Unknown", after = 0L)
    )

lev <- levels(dx\$Genus2)

comp_genus <- dx %>%
    ggplot(aes(x = group, y = Abundance, fill = Genus2)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = rev(cols), labels = lev) +
    guides(fill = guide_legend(title = "Genus", ncol = 1)) +
    labs(y = "Composition (%)", x = "", title = "Genus") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_composition()
ggsave(plot = comp_genus, paste0("composition_genus_", postfix, ".pdf"))

#### Family level ####
N <- 20
d <- tab_comp %>%
    pivot_longer(-Sample, names_to = 'ASV', values_to = 'Abundance') %>%
    mutate(sampleID = Sample)
d\$Family <- tax\$Family[match(d\$ASV, tax\$ASV)]

top_taxa <- d %>%
    group_by(Family, Sample) %>%
    summarise(Abundance = sum(Abundance)) %>%
    group_by(Family) %>%
    summarise(Abundance = mean(Abundance)) %>%
    arrange(-Abundance) %>%
    dplyr::select(Family) %>%
    filter(Family != 'ambiguous') %>%
    head(N) %>%
    unlist()

dx <- d %>%
    mutate(
        Family2 = case_when(
            Family %in% top_taxa           ~ paste(Family),
            is.na(Family)                  ~ paste("Unknown"),
            !(Family %in% top_taxa)        ~ paste("Other families")
        ),
        Family2 = as.factor(Family2)
    ) %>%
    group_by(Family2, Sample) %>%
    summarise(Abundance = sum(Abundance)) %>%
    group_by(Family2) %>%
    summarise(Abundance = mean(Abundance)) %>%
    mutate(
        group   = "all samples",
        Family2 = fct_reorder(Family2, Abundance),
        Family2 = fct_relevel(Family2, "Other families", after = 0L),
        Family2 = fct_relevel(Family2, "Unknown", after = 0L)
    )

lev <- levels(dx\$Family2)

comp_family <- dx %>%
    ggplot(aes(x = group, y = Abundance, fill = Family2)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = rev(cols), labels = lev) +
    guides(fill = guide_legend(title = "Family", ncol = 1)) +
    labs(y = "Composition (%)", x = "", title = "Family") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_composition()
ggsave(plot = comp_family, paste0("composition_family_", postfix, ".pdf"))

#### Phylum level ####
N <- 6
d <- tab_comp %>%
    pivot_longer(-Sample, names_to = 'ASV', values_to = 'Abundance') %>%
    mutate(sampleID = Sample)
d\$Phylum <- tax\$Phylum[match(d\$ASV, tax\$ASV)]

top_taxa <- d %>%
    group_by(Phylum, Sample) %>%
    summarise(Abundance = sum(Abundance)) %>%
    group_by(Phylum) %>%
    summarise(Abundance = mean(Abundance)) %>%
    arrange(-Abundance) %>%
    dplyr::select(Phylum) %>%
    filter(Phylum != 'ambiguous') %>%
    head(N) %>%
    unlist()
top_taxa

dx <- d %>%
    mutate(
        Phylum2 = case_when(
            Phylum %in% top_taxa           ~ paste(Phylum),
            is.na(Phylum)                  ~ paste("Unknown"),
            !(Phylum %in% top_taxa)        ~ paste("Other phyla")
        ),
        Phylum2 = as.factor(Phylum2)
    ) %>%
    group_by(Phylum2, Sample) %>%
    summarise(Abundance = sum(Abundance)) %>%
    group_by(Phylum2) %>%
    summarise(Abundance = mean(Abundance)) %>%
    mutate(
        group   = "all samples",
        Phylum2 = fct_reorder(Phylum2, Abundance),
        Phylum2 = fct_relevel(Phylum2, "Other phyla", after = 0L),
        Phylum2 = fct_relevel(Phylum2, "Unknown", after = 0L)
    )

lev     <- levels(dx\$Phylum2)
colphyl <- c(cols[1:N], cols[21:22])

comp_phylum <- dx %>%
    ggplot(aes(x = group, y = Abundance, fill = Phylum2)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = rev(colphyl), labels = lev) +
    guides(fill = guide_legend(title = "Phylum", ncol = 1)) +
    labs(y = "Composition (%)", x = "", title = "Phylum") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_composition()
ggsave(plot = comp_phylum, paste0("composition_phylum_", postfix, ".pdf"))

## Diversity metrics
shannon <- vegan::diversity(tab, index = 'shannon')
shandf  <- as.data.frame(shannon)
shanpl  <- ggplot(shandf, aes(x = shannon)) +
    geom_histogram(color = "black", fill = "royalblue", alpha = 0.8) +
    theme_Publication() +
    xlab("Shannon index") +
    ggtitle("Shannon index")
ggsave(plot = shanpl, paste0("shannon_index_", postfix, ".pdf"))

specrich <- specnumber(tab)
dfspec   <- as.data.frame(specrich)
specpl   <- ggplot(dfspec, aes(x = specrich)) +
    geom_histogram(color = "black", fill = "darkgreen", alpha = 0.8) +
    theme_Publication() +
    xlab("Number of species") +
    ggtitle("Species richness")
ggsave(plot = specpl, paste0("species_richness_", postfix, ".pdf"))

report <- paste0(
    "This dataset has ", nsamples(phy), " samples and ", ntaxa(phy), " taxa.\n",
    "Species richness: mean ", mean(specrich), ", sd ", sd(specrich),
    ", min ", min(specrich), ", max ", max(specrich), "\n",
    "Shannon index: mean ", mean(shannon), " sd ", sd(shannon), "."
)
writeLines(report, paste0("metrics_overview_", postfix, ".txt"))

################################################
## VERSIONS                                   ##
################################################

r.version       <- paste(R.version[['major']], R.version[['minor']], sep = '.')
phyloseq.version <- as.character(packageVersion('phyloseq'))
vegan.version    <- as.character(packageVersion('vegan'))
ggplot2.version  <- as.character(packageVersion('ggplot2'))

writeLines(
    c('"${task.process}":',
      paste0('    r-base: ',   r.version),
      paste0('    phyloseq: ', phyloseq.version),
      paste0('    vegan: ',    vegan.version),
      paste0('    ggplot2: ',  ggplot2.version)
    ),
    'versions.yml'
)
