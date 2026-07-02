#!/usr/bin/env Rscript
################################################################################
# PERMANOVA: Treatment × Time effects using Fisher-Rao and Bray-Curtis
# at four taxonomic levels: Phylum, Class, Family, Genus
#
# For each level, runs adonis2 on all three streams (Amplicon, MEGAHIT Contigs,
# Read Based) under both FR and BC distance matrices.
#
# Produces:
#   - permanova_fr_bc_multilevel.csv   : summary table (Level × Stream × Metric)
#   - pcoa_permanova_{level}.png/pdf   : 3-stream × 2-metric PCoA grid per level
#
# Experimental design:
#   Treatment: EW (earthworm) vs CT (control), n=6 each
#   Time: 21 days (active) vs 63 days (maturation), n=6 each
#   Reactors: 7,9,11 (EW) and 8,10,12 (CT), 3 replicates per cell
################################################################################

library(vegan)
library(ggplot2)
library(dplyr)
library(gridExtra)

################################################################################
# Configuration
################################################################################

# BASE <- "/GWSPH/groups/cbi/Users/cgaylord/research_data/genomics/vermiculture"
BASE <- "~/src/dissertation/vermiculture_aim1"

streams <- list(
  Amplicon = list(
    seqtab = file.path(BASE, "aim1/amplicon_native/16S/tables/seqtab.rds"),
    taxa   = file.path(BASE, "aim1/amplicon_native/16S/tables/taxa.rds")
  ),
  `MEGAHIT Contigs` = list(
    seqtab = file.path(BASE, "aim1/megahit_dada2/tables/seqtab.rds"),
    taxa   = file.path(BASE, "aim1/megahit_dada2/tables/taxa.rds")
  ),
  `Read Based` = list(
    seqtab = file.path(BASE, "aim1/dada2/all_extracted/tables/seqtab_r2_only.rds"),
    taxa   = file.path(BASE, "aim1/dada2/all_extracted/tables/taxa_r2_only.rds")
  )
)

tax_levels <- c("Phylum", "Class", "Family", "Genus")

output_dir <- file.path(BASE, "aim1/figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# Sample metadata
################################################################################

metadata <- data.frame(
  Sample    = c("ADN1",  "ADN2",  "ADN4",  "ADN5",  "ADN7",  "ADN8",
                "ADN10", "ADN11", "ADN13", "ADN14", "ADN16", "ADN17"),
  Reactor   = c(7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12),
  Age_days  = c(63, 21, 63, 21, 63, 21, 63, 21, 63, 21, 63, 21),
  Treatment = c("EW", "EW", "CT", "CT", "EW", "EW",
                "CT", "CT", "EW", "EW", "CT", "CT"),
  stringsAsFactors = FALSE
)
metadata$Treatment <- factor(metadata$Treatment, levels = c("CT", "EW"))
metadata$Time      <- factor(metadata$Age_days, levels = c(21, 63),
                             labels = c("Day21", "Day63"))
rownames(metadata) <- metadata$Sample

cat("Experimental design:\n")
print(table(metadata$Treatment, metadata$Time))
cat("\n")

################################################################################
# Helper functions
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
  
  # Normalize sample names to ADN format
  rn  <- rownames(agg)
  adn <- regmatches(rn, regexpr("ADN[0-9]+", rn))
  if (length(adn) == nrow(agg)) rownames(agg) <- adn
  
  return(agg)
}

fisher_rao_dist <- function(mat) {
  rel <- mat / rowSums(mat)
  n   <- nrow(rel)
  d   <- matrix(0, n, n)
  rownames(d) <- colnames(d) <- rownames(rel)
  eps <- 1e-10
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      p <- rel[i, ] + eps; p <- p / sum(p)
      q <- rel[j, ] + eps; q <- q / sum(q)
      bc_coef <- min(sum(sqrt(p * q)), 1.0)
      d[i, j] <- d[j, i] <- 2 * acos(bc_coef)
    }
  }
  as.dist(d)
}

################################################################################
# Main loop: taxonomic level → stream → metric
################################################################################

all_summary  <- list()   # collects one row per level × stream × metric
all_results  <- list()   # collects full adonis2 objects for downstream use

for (tax_level in tax_levels) {
  
  cat("\n")
  cat("################################################################\n")
  cat("Taxonomic level:", tax_level, "\n")
  cat("################################################################\n\n")
  
  level_results <- list()   # stream / metric → result list
  
  for (stream_name in names(streams)) {
    
    cat("  ----------------------------------------------------------------\n")
    cat("  Stream:", stream_name, "\n")
    cat("  ----------------------------------------------------------------\n\n")
    
    mat    <- load_stream(streams[[stream_name]], tax_level)
    common <- intersect(rownames(mat), rownames(metadata))
    mat    <- mat[common, , drop = FALSE]
    meta   <- metadata[common, ]
    
    cat("  Samples:", nrow(meta), "  |  ", tax_level, "taxa:", ncol(mat), "\n\n")
    
    for (metric in c("Fisher-Rao", "Bray-Curtis")) {
      
      cat("  ---", metric, "---\n")
      
      if (metric == "Fisher-Rao") {
        d <- fisher_rao_dist(mat)
      } else {
        rel <- mat / rowSums(mat)
        d   <- vegdist(rel, method = "bray")
      }
      
      # Full factorial PERMANOVA
      set.seed(42)
      perm_full <- adonis2(d ~ Treatment * Time, data = meta,
                           permutations = 9999, by = "terms")
      
      # Marginal models
      set.seed(42)
      perm_trt  <- adonis2(d ~ Treatment, data = meta, permutations = 9999)
      set.seed(42)
      perm_time <- adonis2(d ~ Time,      data = meta, permutations = 9999)
      
      cat(sprintf("  Treatment: R²=%.3f p=%.4f  |  Time: R²=%.3f p=%.4f  |  TxT: R²=%.3f p=%.4f\n\n",
                  perm_full["Treatment",    "R2"],
                  perm_full["Treatment",    "Pr(>F)"],
                  perm_full["Time",         "R2"],
                  perm_full["Time",         "Pr(>F)"],
                  perm_full["Treatment:Time","R2"],
                  perm_full["Treatment:Time","Pr(>F)"]))
      
      label <- paste(tax_level, stream_name, metric, sep = " / ")
      
      level_results[[label]] <- list(
        level   = tax_level,
        stream  = stream_name,
        metric  = metric,
        full    = perm_full,
        treatment = perm_trt,
        time    = perm_time,
        dist    = d,
        mat     = mat,
        meta    = meta
      )
      
      all_summary[[label]] <- data.frame(
        Level          = tax_level,
        Stream         = stream_name,
        Metric         = metric,
        Treatment_R2   = round(perm_full["Treatment",     "R2"],      4),
        Treatment_p    = round(perm_full["Treatment",     "Pr(>F)"],  4),
        Time_R2        = round(perm_full["Time",          "R2"],      4),
        Time_p         = round(perm_full["Time",          "Pr(>F)"],  4),
        Interaction_R2 = round(perm_full["Treatment:Time","R2"],      4),
        Interaction_p  = round(perm_full["Treatment:Time","Pr(>F)"],  4)
      )
    }
  }
  
  all_results[[tax_level]] <- level_results
  
  ##############################################################################
  # PCoA figure for this taxonomic level: 3 streams × 2 metrics = 6 panels
  ##############################################################################
  
  cat("  Creating PCoA figure for", tax_level, "...\n")
  
  pcoa_plots <- list()
  
  stream_order <- names(streams)          # Amplicon, MEGAHIT Contigs, Read Based
  metric_order <- c("Fisher-Rao", "Bray-Curtis")
  
  for (stream_name in stream_order) {
    for (metric in metric_order) {
      
      label <- paste(tax_level, stream_name, metric, sep = " / ")
      r     <- level_results[[label]]
      d     <- r$dist
      
      pc           <- cmdscale(d, k = 2, eig = TRUE)
      eig_pos      <- pc$eig[pc$eig > 0]
      var_explained <- round(100 * pc$eig[1:2] / sum(eig_pos), 1)
      
      pcoa_df <- data.frame(
        PC1       = pc$points[, 1],
        PC2       = pc$points[, 2],
        Treatment = metadata[rownames(pc$points), "Treatment"],
        Time      = metadata[rownames(pc$points), "Time"],
        Sample    = rownames(pc$points)
      )
      
      # Pull PERMANOVA p-values for subtitle annotation
      trt_p  <- r$full["Treatment",     "Pr(>F)"]
      time_p <- r$full["Time",          "Pr(>F)"]
      int_p  <- r$full["Treatment:Time","Pr(>F)"]
      
      subtitle_txt <- sprintf("Trt p=%.3f  |  Time p=%.3f  |  TxT p=%.3f",
                              trt_p, time_p, int_p)
      
      p <- ggplot(pcoa_df,
                  aes(x = PC1, y = PC2, color = Treatment, shape = Time)) +
        geom_point(size = 3.5, alpha = 0.85) +
        scale_color_manual(values = c("CT" = "#377EB8", "EW" = "#E41A1C")) +
        scale_shape_manual(values = c("Day21" = 16, "Day63" = 17)) +
        theme_bw(base_size = 10) +
        labs(
          title    = paste0(stream_name, "  —  ", metric),
          subtitle = subtitle_txt,
          x        = paste0("PCoA1 (", var_explained[1], "%)"),
          y        = paste0("PCoA2 (", var_explained[2], "%)")
        ) +
        theme(
          plot.title    = element_text(size = 10, face = "bold"),
          plot.subtitle = element_text(size = 8,  color = "grey40"),
          legend.position = "bottom",
          legend.box = "horizontal"
        )
      
      pcoa_plots[[paste(stream_name, metric)]] <- p
    }
  }
  
  # Arrange: rows = streams (Amplicon, MEGAHIT, Reads), cols = metrics (FR, BC)
  plot_order <- c()
  for (s in stream_order) {
    for (m in metric_order) {
      plot_order <- c(plot_order, paste(s, m))
    }
  }
  
  p_grid <- arrangeGrob(
    grobs = pcoa_plots[plot_order],
    ncol  = 2,
    top   = paste0("PCoA Ordination: Treatment × Time — ", tax_level, " Level")
  )
  
  outbase <- file.path(output_dir, paste0("pcoa_permanova_", tolower(tax_level)))
  ggsave(paste0(outbase, ".png"), p_grid, width = 12, height = 14, dpi = 300)
  ggsave(paste0(outbase, ".pdf"), p_grid, width = 12, height = 14)

  p_grid <- arrangeGrob(
    grobs = pcoa_plots[plot_order],
    ncol  = 2
    #,    top   = paste0("PCoA Ordination: Treatment × Time — ", tax_level, " Level")
  )
  
  outbase <- file.path(output_dir, paste0("pcoa_permanova_", tolower(tax_level), "-notitle"))
  ggsave(paste0(outbase, ".png"), p_grid, width = 12, height = 14, dpi = 300)
  ggsave(paste0(outbase, ".pdf"), p_grid, width = 12, height = 14)
  
  
  cat("  Saved:", paste0("pcoa_permanova_", tolower(tax_level), ".png/pdf\n\n"))
}

################################################################################
# Cross-level summary table
################################################################################

cat("\n################################################################\n")
cat("CROSS-LEVEL SUMMARY\n")
cat("################################################################\n\n")

summary_df <- do.call(rbind, all_summary)
rownames(summary_df) <- NULL

# Order for printing
summary_df$Level  <- factor(summary_df$Level,  levels = tax_levels)
summary_df$Stream <- factor(summary_df$Stream, levels = names(streams))
summary_df$Metric <- factor(summary_df$Metric, levels = c("Fisher-Rao", "Bray-Curtis"))
summary_df <- summary_df[order(summary_df$Level, summary_df$Stream, summary_df$Metric), ]

cat(sprintf("%-10s  %-20s  %-12s  Trt: R²    p      Time: R²   p      TxT: R²    p\n",
            "Level", "Stream", "Metric"))
cat(strrep("-", 95), "\n")

for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  cat(sprintf("%-10s  %-20s  %-12s  %.3f  %.4f    %.3f  %.4f    %.3f  %.4f\n",
              as.character(r$Level),
              as.character(r$Stream),
              as.character(r$Metric),
              r$Treatment_R2, r$Treatment_p,
              r$Time_R2,      r$Time_p,
              r$Interaction_R2, r$Interaction_p))
}

out_csv <- file.path(output_dir, "permanova_fr_bc_multilevel.csv")
write.csv(summary_df, out_csv, row.names = FALSE)
cat("\nSaved:", out_csv, "\n")

cat("\n=== Done ===\n")
