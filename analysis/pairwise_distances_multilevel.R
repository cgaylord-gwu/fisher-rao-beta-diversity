#!/usr/bin/env Rscript
################################################################################
# Pairwise FR and BC Distance Comparison — Phylum, Family, Genus
#
# Runs the same analysis at three taxonomic levels in one pass.
################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)

################################################################################
# Configuration
################################################################################

BASE <- "."

streams <- list(
  Amplicon = list(
    seqtab = file.path(BASE, "data/amplicon_native/16S/tables/seqtab.rds"),
    taxa   = file.path(BASE, "data/amplicon_native/16S/tables/taxa.rds")
  ),
  `Read Based` = list(
    seqtab = file.path(BASE, "data/extracted_read/tables/seqtab.rds"),
    taxa   = file.path(BASE, "data/extracted_read/tables/taxa.rds")
  ),
  `MEGAHIT Contigs` = list(
    seqtab = file.path(BASE, "data/megahit_contig/tables/seqtab.rds"),
    taxa   = file.path(BASE, "data/megahit_contig/tables/taxa.rds")
  )
)

tax_levels <- c("Phylum", "Class", "Family", "Genus")

output_dir <- file.path(BASE, "figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# Distance functions
################################################################################

fisher_rao <- function(p, q) {
  eps <- 1e-10
  p <- (p + eps) / sum(p + eps)
  q <- (q + eps) / sum(q + eps)
  bc <- sum(sqrt(p * q))
  bc <- min(bc, 1.0)
  2 * acos(bc)
}

bray_curtis <- function(p, q) {
  sum(abs(p - q)) / sum(p + q)
}

################################################################################
# Load and aggregate
################################################################################

load_stream <- function(stream_info, tax_level) {
  seqtab <- readRDS(stream_info$seqtab)
  taxa   <- readRDS(stream_info$taxa)
  
  shared <- intersect(colnames(seqtab), rownames(taxa))
  seqtab <- seqtab[, shared, drop = FALSE]
  taxa   <- taxa[shared, , drop = FALSE]
  
  tax_assign <- taxa[, tax_level]
  tax_assign[is.na(tax_assign)] <- "Unassigned"
  
  agg <- matrix(0, nrow = nrow(seqtab), ncol = length(unique(tax_assign)))
  colnames(agg) <- sort(unique(tax_assign))
  rownames(agg) <- rownames(seqtab)
  
  for (taxon in colnames(agg)) {
    idx <- which(tax_assign == taxon)
    if (length(idx) == 1) {
      agg[, taxon] <- seqtab[, idx]
    } else {
      agg[, taxon] <- rowSums(seqtab[, idx, drop = FALSE])
    }
  }
  return(agg)
}

normalize_sample_names <- function(mat) {
  rn <- rownames(mat)
  adn <- regmatches(rn, regexpr("ADN[0-9]+", rn))
  if (length(adn) == nrow(mat)) {
    rownames(mat) <- adn
  }
  return(mat)
}

harmonize <- function(mat, all_cols) {
  missing <- setdiff(all_cols, colnames(mat))
  if (length(missing) > 0) {
    extra <- matrix(0, nrow = nrow(mat), ncol = length(missing))
    colnames(extra) <- missing
    mat <- cbind(mat, extra)
  }
  mat[, sort(all_cols), drop = FALSE]
}

################################################################################
# Main loop over taxonomic levels
################################################################################

all_results <- data.frame()

for (tax_level in tax_levels) {
  cat("\n========================================\n")
  cat("Taxonomic level:", tax_level, "\n")
  cat("========================================\n\n")
  
  # Load
  level_data <- list()
  for (nm in names(streams)) {
    mat <- load_stream(streams[[nm]], tax_level)
    mat <- normalize_sample_names(mat)
    level_data[[nm]] <- mat
    cat("  ", nm, ":", ncol(mat), "taxa,", nrow(mat), "samples\n")
  }
  
  common_samples <- Reduce(intersect, lapply(level_data, rownames))
  cat("  Common samples:", length(common_samples), "\n")
  
  for (nm in names(level_data)) {
    level_data[[nm]] <- level_data[[nm]][common_samples, , drop = FALSE]
  }
  
  # Harmonize
  all_taxa <- unique(unlist(lapply(level_data, colnames)))
  level_harmonized <- lapply(level_data, harmonize, all_taxa)
  cat("  Total taxa (union):", length(all_taxa), "\n\n")
  
  # Compute distances
  stream_names <- names(level_harmonized)
  
  for (i in 1:(length(stream_names) - 1)) {
    for (j in (i + 1):length(stream_names)) {
      s1 <- stream_names[i]
      s2 <- stream_names[j]
      pair_label <- paste(s1, "vs", s2)
      
      for (samp in common_samples) {
        p <- level_harmonized[[s1]][samp, ]
        q <- level_harmonized[[s2]][samp, ]
        
        fr <- fisher_rao(p, q)
        bc <- bray_curtis(p, q)
        
        all_results <- rbind(all_results, data.frame(
          Level = tax_level,
          Pair = pair_label,
          Sample = samp,
          Metric = "Fisher-Rao",
          Distance = fr
        ))
        all_results <- rbind(all_results, data.frame(
          Level = tax_level,
          Pair = pair_label,
          Sample = samp,
          Metric = "Bray-Curtis",
          Distance = bc
        ))
      }
      
      fr_vals <- all_results$Distance[all_results$Level == tax_level &
                                        all_results$Pair == pair_label &
                                        all_results$Metric == "Fisher-Rao"]
      bc_vals <- all_results$Distance[all_results$Level == tax_level &
                                        all_results$Pair == pair_label &
                                        all_results$Metric == "Bray-Curtis"]
      
      cat(sprintf("  %-35s  FR: %.3f (%.3f-%.3f)  BC: %.3f (%.3f-%.3f)\n",
                  pair_label, mean(fr_vals), min(fr_vals), max(fr_vals),
                  mean(bc_vals), min(bc_vals), max(bc_vals)))
    }
  }
  
  # CV comparison
  cat("\n  Variance (CV):\n")
  pair_order <- c("Amplicon vs MEGAHIT Contigs",
                  "Amplicon vs Read Based",
                  "Read Based vs MEGAHIT Contigs")
  
  for (pair in pair_order) {
    fr_vals <- all_results$Distance[all_results$Level == tax_level &
                                      all_results$Pair == pair &
                                      all_results$Metric == "Fisher-Rao"]
    bc_vals <- all_results$Distance[all_results$Level == tax_level &
                                      all_results$Pair == pair &
                                      all_results$Metric == "Bray-Curtis"]
    
    if (length(fr_vals) > 0 && length(bc_vals) > 0) {
      fr_cv <- sd(fr_vals) / mean(fr_vals)
      bc_cv <- sd(bc_vals) / mean(bc_vals)
      cat(sprintf("  %-35s  FR CV: %.3f  BC CV: %.3f  (ratio: %.2f)\n",
                  pair, fr_cv, bc_cv, fr_cv / bc_cv))
    }
  }
}

################################################################################
# Save full results
################################################################################

write.csv(all_results, file.path(output_dir, "pairwise_distances_all_levels.csv"),
          row.names = FALSE)

################################################################################
# Multi-panel figure: one row per taxonomic level
################################################################################

cat("\n\nCreating multi-level boxplot...\n")

pair_order <- c("Amplicon vs MEGAHIT Contigs",
                "Amplicon vs Read Based",
                "Read Based vs MEGAHIT Contigs")
pair_order <- pair_order[pair_order %in% unique(all_results$Pair)]
all_results$Pair <- factor(all_results$Pair, levels = pair_order)
all_results$Level <- factor(all_results$Level, levels = tax_levels)

# Short pair labels for readability
pair_labels <- c(
  "Amplicon vs MEGAHIT Contigs" = "Amp vs MEGAHIT",
  "Amplicon vs Read Based" = "Amp vs Reads",
  "Read Based vs MEGAHIT Contigs" = "Reads vs MEGAHIT"
)
all_results$Pair_Short <- factor(pair_labels[as.character(all_results$Pair)],
                                 levels = pair_labels[pair_order])


p_multi <- ggplot(all_results, aes(x = Pair_Short, y = Distance, fill = Pair)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.5) +
  # facet_grid(Level ~ Metric, scales = "free") +
  facet_wrap(~ Level + Metric, scales = "free", ncol = 2) +
  scale_fill_manual(values = c(
    "Amplicon vs MEGAHIT Contigs" = "#4DAF4A",
    "Amplicon vs Read Based" = "#E41A1C",
    "Read Based vs MEGAHIT Contigs" = "#377EB8"
  )) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 11)
  ) +
  labs(
    title = "Pairwise Distance Between Extraction Methods by Taxonomic Level",
    subtitle = "Per-sample paired distances (n = 12 per pair)",
    x = NULL,
    y = "Distance"
  )

ggsave(file.path(output_dir, "pairwise_distances_multilevel.png"),
       p_multi, width = 10, height = 10, dpi = 300)
ggsave(file.path(output_dir, "pairwise_distances_multilevel.pdf"),
       p_multi, width = 10, height = 10, dpi = 300)

p_multi <- ggplot(all_results, aes(x = Pair_Short, y = Distance, fill = Pair)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.5) +
  # facet_grid(Level ~ Metric, scales = "free") +
  facet_wrap(~ Level + Metric, scales = "free", ncol = 2) +
  scale_fill_manual(values = c(
    "Amplicon vs MEGAHIT Contigs" = "#4DAF4A",
    "Amplicon vs Read Based" = "#E41A1C",
    "Read Based vs MEGAHIT Contigs" = "#377EB8"
  )) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 11)
  ) +
  labs(
    # title = "Pairwise Distance Between Extraction Methods by Taxonomic Level",
    # subtitle = "Per-sample paired distances (n = 12 per pair)",
    x = NULL,
    y = "Distance"
  )

ggsave(file.path(output_dir, "pairwise_distances_multilevel-notitle.png"),
       p_multi, width = 10, height = 10, dpi = 300)
ggsave(file.path(output_dir, "pairwise_distances_multilevel-notitle.pdf"),
       p_multi, width = 10, height = 10, dpi = 300)

cat("Saved: pairwise_distances_multilevel.png/pdf\n")

cat("\n=== Done ===\n")
