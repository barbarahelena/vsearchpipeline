process PHYLOSEQ_DECONTAM {
    label 'process_multi_low'
    label 'phyloseq'

    input:
    path phyloseq
    path samplesheet

    output:
    path "phyloseq_decontam.RDS"                                                                                                         , emit: phyloseq
    path "decontam_report.txt"                                                                                                           , emit: report
    path "decontam_contaminants.csv"                                                                                                     , emit: contaminants
    path "decontam_prev_plot.pdf"                                                                                                        , emit: prev_plot
    tuple val("${task.process}"), val('R'),        eval('Rscript -e "cat(paste(R.version[c(\'major\',\'minor\')], collapse=\'.\'))"')  , emit: versions_r,        topic: versions
    tuple val("${task.process}"), val('phyloseq'), eval('Rscript -e "cat(as.character(packageVersion(\'phyloseq\')))"')                , emit: versions_phyloseq, topic: versions
    tuple val("${task.process}"), val('decontam'), eval('Rscript -e "cat(as.character(packageVersion(\'decontam\')))"')                , emit: versions_decontam, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    Rscript - <<EOF
    library(phyloseq)
    library(decontam)
    library(ggplot2)

    ## -------------------------------------------------------
    ## Load phyloseq object
    ## -------------------------------------------------------
    ps <- readRDS("$phyloseq")
    report_lines <- c()
    report_lines <- c(report_lines,
        paste0("Phyloseq object loaded: ", nsamples(ps), " samples, ", ntaxa(ps), " taxa."))

    ## -------------------------------------------------------
    ## Load samplesheet to get nucl_acid_conc if present
    ## -------------------------------------------------------
    meta <- read.csv("$samplesheet", stringsAsFactors = FALSE)
    has_conc <- "nucl_acid_conc" %in% colnames(meta) &&
                any(!is.na(meta\$nucl_acid_conc) & meta\$nucl_acid_conc > 0)

    ## -------------------------------------------------------
    ## Detect negative controls (sample names containing NEGCON)
    ## -------------------------------------------------------
    neg_pattern <- "NEGCON"
    negctrlsvec <- sample_names(ps)[grepl(neg_pattern, sample_names(ps))]
    has_negctrls <- length(negctrlsvec) > 0

    report_lines <- c(report_lines,
        paste0("Negative controls found (pattern '", neg_pattern, "'): ", length(negctrlsvec)),
        paste0("DNA concentration column available: ", has_conc))

    if (!has_negctrls && !has_conc) {
        report_lines <- c(report_lines,
            "WARNING: No negative controls (NEGCON) and no nucl_acid_conc column found.",
            "Skipping decontam. Saving original phyloseq as output.")
        writeLines(report_lines, "decontam_report.txt")
        write.csv(data.frame(ASV = character(0), method = character(0)),
                  "decontam_contaminants.csv", row.names = FALSE)
        pdf("decontam_prev_plot.pdf"); dev.off()
        saveRDS(ps, "phyloseq_decontam.RDS")
    } else {

        ## -------------------------------------------------------
        ## Add control flag and concentration to sample_data
        ## -------------------------------------------------------
        df_sd <- data.frame(sample = sample_names(ps), stringsAsFactors = FALSE)
        df_sd\$Ctrl <- df_sd\$sample %in% negctrlsvec

        if (has_conc) {
            meta_sub <- meta[, c("sample", "nucl_acid_conc")]
            # match by sample name
            df_sd <- merge(df_sd, meta_sub, by = "sample", all.x = TRUE)
        }
        rownames(df_sd) <- df_sd\$sample
        df_sd\$sample <- NULL

        # Merge with existing sample_data
        existing_sd <- data.frame(sample_data(ps))
        for (col in colnames(df_sd)) {
            existing_sd[[col]] <- df_sd[rownames(existing_sd), col]
        }
        sample_data(ps) <- sample_data(existing_sd)

        cont_asvs <- c()

        ## -------------------------------------------------------
        ## Method 1: Frequency-based (requires nucl_acid_conc)
        ## -------------------------------------------------------
        if (has_conc) {
            ps_freq <- prune_samples(
                !is.na(sample_data(ps)\$nucl_acid_conc) & sample_data(ps)\$nucl_acid_conc > 0,
                ps)
            contam_freq <- isContaminant(ps_freq, method = "frequency",
                                         conc = "nucl_acid_conc")
            cont_freq <- rownames(contam_freq)[which(contam_freq\$contaminant == TRUE)]
            n_freq <- sum(contam_freq\$contaminant, na.rm = TRUE)
            report_lines <- c(report_lines,
                paste0("Frequency-based decontam: ", n_freq, " contaminant ASVs identified."))
            cont_asvs <- union(cont_asvs, cont_freq)
        } else {
            report_lines <- c(report_lines,
                "Frequency-based decontam skipped (no nucl_acid_conc column).")
        }

        ## -------------------------------------------------------
        ## Method 2: Prevalence-based (requires negative controls)
        ## -------------------------------------------------------
        if (has_negctrls) {
            contamdf_prev <- isContaminant(ps, method = "prevalence", neg = "Ctrl")
            cont_prev <- rownames(contamdf_prev)[which(contamdf_prev\$contaminant == TRUE)]
            n_prev <- sum(contamdf_prev\$contaminant, na.rm = TRUE)
            report_lines <- c(report_lines,
                paste0("Prevalence-based decontam: ", n_prev, " contaminant ASVs identified."))
            cont_asvs <- union(cont_asvs, cont_prev)

            ## Prevalence plot
            ps_pa <- transform_sample_counts(ps, function(abund) 1 * (abund > 0))
            ps_neg_pa <- prune_samples(sample_data(ps_pa)\$Ctrl == TRUE, ps_pa)
            ps_smp_pa <- prune_samples(sample_data(ps_pa)\$Ctrl == FALSE, ps_pa)
            dfpa <- data.frame(
                pssample    = taxa_sums(ps_smp_pa),
                psneg       = taxa_sums(ps_neg_pa),
                contaminant = contamdf_prev\$contaminant
            )
            prev_plot <- ggplot(data = dfpa,
                                aes(x = psneg, y = pssample, color = contaminant)) +
                geom_point(alpha = 0.6) +
                scale_color_manual(values = c("royalblue4", "firebrick2")) +
                xlab("Prevalence (Negative Controls)") +
                ylab("Prevalence (True Samples)") +
                ggtitle("Decontam: prevalence-based contamination") +
                theme_minimal()
            ggsave("decontam_prev_plot.pdf", plot = prev_plot)
        } else {
            report_lines <- c(report_lines,
                "Prevalence-based decontam skipped (no NEGCON samples found).")
            pdf("decontam_prev_plot.pdf"); dev.off()
        }

        ## -------------------------------------------------------
        ## Remove contaminants
        ## -------------------------------------------------------
        report_lines <- c(report_lines,
            paste0("Total unique contaminant ASVs to remove: ", length(cont_asvs)))
        ps_noncontam <- prune_taxa(!taxa_names(ps) %in% cont_asvs, ps)

        ## -------------------------------------------------------
        ## Remove negative control samples
        ## -------------------------------------------------------
        neg_and_sb <- sample_names(ps_noncontam)[
            grepl("NEGCON|^SB", sample_names(ps_noncontam))
        ]
        if (length(neg_and_sb) > 0) {
            ps_noncontam <- prune_samples(
                !sample_names(ps_noncontam) %in% neg_and_sb,
                ps_noncontam)
            report_lines <- c(report_lines,
                paste0("Removed ", length(neg_and_sb),
                       " negative/blank control samples: ",
                       paste(neg_and_sb, collapse = ", ")))
        }

        report_lines <- c(report_lines,
            paste0("Decontam output: ", nsamples(ps_noncontam),
                   " samples, ", ntaxa(ps_noncontam), " taxa."))

        ## -------------------------------------------------------
        ## Write outputs
        ## -------------------------------------------------------
        write.csv(data.frame(ASV = cont_asvs, stringsAsFactors = FALSE),
                  "decontam_contaminants.csv", row.names = FALSE)
        writeLines(report_lines, "decontam_report.txt")
        saveRDS(ps_noncontam, "phyloseq_decontam.RDS")
    }
    EOF
    """

    stub:
    """
    touch phyloseq_decontam.RDS
    touch decontam_report.txt
    touch decontam_contaminants.csv
    touch decontam_prev_plot.pdf
    """
}
