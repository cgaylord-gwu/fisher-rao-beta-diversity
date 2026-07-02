#!/usr/bin/env Rscript
################################################################################
# Three-Stream Community Composition Barplot — All Taxonomic Levels
#
# Clark Gaylord
# Computational Biology Institute, The George Washington University
#
# Developed with assistance from Claude (Anthropic).
#
# Generates stacked relative-abundance barplots for all four taxonomic
# levels (Phylum, Class, Family, Genus) from three analytical streams
# (Amplicon, MEGAHIT Contigs, Read-Based), organized by experimental
# cell (CT/EW x Day21/63).
#
# Output files (per level, e.g. for Phylum):
#   three_stream_experimental_design_barplot_phylum.pdf/.png
#   three_stream_experimental_design_barplot_phylum-notitle.pdf/.png
#
# Usage:
#   Source in RStudio; all four taxonomic levels are processed in sequence.
#   To restrict to a subset, edit tax_levels below.
################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(scales)

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a)) a else b

BASE <- "."

streams <- list(
  "Amplicon" = list(
    seqtab = file.path(BASE, "data/amplicon_native/16S/tables/seqtab.rds"),
    taxa   = file.path(BASE, "data/amplicon_native/16S/tables/taxa.rds")
  ),
  "MEGAHIT Contigs" = list(
    seqtab = file.path(BASE, "data/megahit_contig/tables/seqtab.rds"),
    taxa   = file.path(BASE, "data/megahit_contig/tables/taxa.rds")
  ),
  "Read-Based" = list(
    seqtab = file.path(BASE, "data/extracted_read/tables/seqtab.rds"),
    taxa   = file.path(BASE, "data/extracted_read/tables/taxa.rds")
  )
)

output_dir <- file.path(BASE, "figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

N_TOP <- 10  # number of top taxa to show before collapsing to "Other"

################################################################################
# Metadata
################################################################################

metadata <- data.frame(
  Sample    = c("ADN1",  "ADN2",  "ADN4",  "ADN5",  "ADN7",  "ADN8",
                "ADN10", "ADN11", "ADN13", "ADN14", "ADN16", "ADN17"),
  Treatment = c("EW","EW","CT","CT","EW","EW","CT","CT","EW","EW","CT","CT"),
  Time      = c(63, 21, 63, 21, 63, 21, 63, 21, 63, 21, 63, 21),
  stringsAsFactors = FALSE
)
metadata$Group <- paste0(metadata$Treatment, "_Day", metadata$Time)
group_order <- c("CT_Day21","CT_Day63","EW_Day21","EW_Day63")
metadata$Group <- factor(metadata$Group, levels = group_order)
metadata <- metadata[order(metadata$Group, metadata$Sample), ]
sample_order <- metadata$Sample

tax_levels <- c("Phylum", "Class", "Family", "Genus")

################################################################################
# Core function: load and aggregate to a given taxonomic level
################################################################################

load_stream <- function(stream_info, tax_level) {
  seqtab <- readRDS(stream_info$seqtab)
  taxa   <- readRDS(stream_info$taxa)
  
  # Normalize sample names to ADN format
  rn  <- rownames(seqtab)
  adn <- regmatches(rn, regexpr("ADN[0-9]+", rn))
  if (length(adn) == nrow(seqtab)) rownames(seqtab) <- adn
  
  # Extract the requested taxonomic level column
  if (!tax_level %in% colnames(taxa))
    stop(sprintf("Level '%s' not found in taxa table. Available: %s",
                 tax_level, paste(colnames(taxa), collapse = ", ")))
  
  tax_vec <- taxa[, tax_level]
  tax_vec[is.na(tax_vec)] <- "Unassigned"
  
  # Aggregate counts by taxon
  taxa_names <- unique(tax_vec)
  agg <- matrix(0, nrow = nrow(seqtab),
                ncol = length(taxa_names),
                dimnames = list(rownames(seqtab), taxa_names))
  for (tx in taxa_names) {
    idx <- which(tax_vec == tx)
    agg[, tx] <- rowSums(seqtab[, idx, drop = FALSE])
  }
  agg
}

################################################################################
# Build stacked data frame for one level
################################################################################

build_stack <- function(level_data, top_taxa, tax_level) {
  do.call(rbind, lapply(names(level_data), function(nm) {
    mat     <- level_data[[nm]]
    n_reads <- rowSums(mat)
    rel     <- mat / n_reads
    
    df   <- as.data.frame(rel)
    keep <- intersect(top_taxa, colnames(df))
    others <- setdiff(colnames(df), c(keep, "Unassigned"))
    
    df2 <- df[, keep, drop = FALSE]
    df2$Other <- if (length(others) > 0) rowSums(df[, others, drop = FALSE]) else 0
    if ("Unassigned" %in% colnames(df)) df2$Unassigned <- df[, "Unassigned"]
    
    df2$Sample  <- rownames(df2)
    df2$Method  <- nm
    df2$N_reads <- n_reads[rownames(df2)]
    
    pivot_longer(df2,
                 cols      = -c(Sample, Method, N_reads),
                 names_to  = tax_level,
                 values_to = "Abundance")
  }))
}

################################################################################
# Plot function for one level
################################################################################

make_barplot <- function(stack_data, tax_level, taxon_order, pal, title = TRUE) {
  stack_data[[tax_level]] <- factor(stack_data[[tax_level]], levels = rev(taxon_order))
  
  p <- ggplot(stack_data,
              aes(x = Sample, y = Abundance, fill = .data[[tax_level]])) +
    geom_bar(stat = "identity", position = "stack", width = 0.85) +
    facet_grid(Method ~ Group, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = pal, name = tax_level) +
    scale_y_continuous(labels = percent_format(), expand = c(0, 0),
                       limits = c(0, 1.05)) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y        = element_text(size = 8),
      strip.text.x       = element_text(face = "bold", size = 10),
      strip.text.y       = element_text(face = "bold", size = 9),
      legend.position    = "right",
      legend.key.size    = unit(0.6, "cm"),
      legend.text        = element_text(size = 11),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.spacing.x    = unit(0.25, "cm"),
      panel.spacing.y    = unit(0.4, "cm")
    ) +
    labs(x = NULL, y = "Relative Abundance")
  
  if (title) {
    p <- p + labs(
      title    = sprintf("%s-Level Community Composition by Extraction Method",
                         tax_level),
      subtitle = "CT = Control, EW = Earthworm"
    )
  }
  p
}

################################################################################
# Main: loop over all four taxonomic levels
################################################################################

for (tax_level in tax_levels) {
  cat(sprintf("\nProcessing: %s\n", tax_level))
  
  # Load and aggregate
  level_data <- lapply(streams, load_stream, tax_level = tax_level)
  
  # Common samples across streams
  common     <- Reduce(intersect, lapply(level_data, rownames))
  level_data <- lapply(level_data, function(m) m[common, , drop = FALSE])
  
  # Top N taxa (excluding Unassigned and Other) across all streams
  # Use a named union approach to handle different taxa sets per stream
  all_taxa <- Reduce(union, lapply(level_data, colnames))
  all_totals <- setNames(numeric(length(all_taxa)), all_taxa)
  for (m in level_data) {
    shared <- intersect(names(all_totals), colnames(m))
    all_totals[shared] <- all_totals[shared] + colSums(m)[shared]
  }
  all_totals <- all_totals[!names(all_totals) %in% c("Unassigned", "Other")]
  top_taxa <- names(sort(all_totals, decreasing = TRUE))[
    seq_len(min(N_TOP, sum(all_totals > 0)))]
  
  # Build stacked data frame
  stack_data <- build_stack(level_data, top_taxa, tax_level)
  stack_data <- left_join(stack_data,
                          metadata[, c("Sample","Group","Treatment","Time")],
                          by = "Sample")
  stack_data$Sample <- factor(stack_data$Sample, levels = sample_order)
  stack_data$Method <- factor(stack_data$Method,
                              levels = c("Amplicon","MEGAHIT Contigs","Read-Based"))
  stack_data$Group  <- factor(stack_data$Group, levels = group_order)
  
  # Build taxon_order ensuring no duplicates
  present <- unique(stack_data[[tax_level]])
  taxon_order <- c(
    top_taxa[top_taxa %in% present],
    intersect("Other", present),
    intersect("Unassigned", present)
  )
  taxon_order <- unique(taxon_order)  # guard against any remaining duplicates
  
  # Colors
  n_col <- length(taxon_order)
  pal   <- if (n_col <= 12) {
    brewer.pal(max(3, n_col), "Set3")[seq_len(n_col)]
  } else {
    colorRampPalette(brewer.pal(12, "Set3"))(n_col)
  }
  names(pal)        <- taxon_order
  pal["Other"]      <- "grey70"
  pal["Unassigned"] <- "grey90"
  
  # Base filename
  base <- file.path(output_dir,
                    sprintf("three_stream_experimental_design_barplot_%s",
                            tolower(tax_level)))
  
  # Titled version
  p_title <- make_barplot(stack_data, tax_level, taxon_order, pal, title = TRUE)
  ggsave(paste0(base, ".pdf"), p_title, width = 18, height = 11, dpi = 300)
  ggsave(paste0(base, ".png"), p_title, width = 18, height = 11, dpi = 300)
  cat(sprintf("  Saved: %s.pdf/.png\n", basename(base)))
  
  # No-title version for manuscript figures
  p_notitle <- make_barplot(stack_data, tax_level, taxon_order, pal, title = FALSE)
  ggsave(paste0(base, "-notitle.pdf"), p_notitle, width = 18, height = 11, dpi = 300)
  ggsave(paste0(base, "-notitle.png"), p_notitle, width = 18, height = 11, dpi = 300)
  cat(sprintf("  Saved: %s-notitle.pdf/.png\n", basename(base)))
}

cat("\nDone.\n")