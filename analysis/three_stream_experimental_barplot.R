#!/usr/bin/env Rscript
################################################################################
# Three-Stream Experimental Design Barplot
# 3 rows (Amplicon, MEGAHIT Contigs, Read-Based) x
# 4 columns (CT_Day21, CT_Day63, EW_Day21, EW_Day63)
# 3 replicates per cell
# Includes read count annotation per sample
################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(scales)

# Define %||% if rlang not loaded
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

N_TOP_PHYLA <- 10

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

# Order samples within group
metadata <- metadata[order(metadata$Group, metadata$Sample), ]
sample_order <- metadata$Sample

################################################################################
# Helper: load and aggregate to phylum
################################################################################

load_phylum <- function(stream_info) {
  seqtab <- readRDS(stream_info$seqtab)
  taxa   <- readRDS(stream_info$taxa)
  
  # Normalize sample names to ADN format
  rn <- rownames(seqtab)
  adn <- regmatches(rn, regexpr("ADN[0-9]+", rn))
  if (length(adn) == nrow(seqtab)) rownames(seqtab) <- adn
  
  # Aggregate by phylum
  phyla <- taxa[, "Phylum"]
  phyla[is.na(phyla)] <- "Unassigned"
  agg <- matrix(0, nrow = nrow(seqtab),
                ncol = length(unique(phyla)),
                dimnames = list(rownames(seqtab), unique(phyla)))
  for (p in unique(phyla)) {
    idx <- which(phyla == p)
    agg[, p] <- rowSums(seqtab[, idx, drop = FALSE])
  }
  agg
}

################################################################################
# Load all three streams
################################################################################

phylum_data <- lapply(streams, load_phylum)

# Keep only common samples
common <- Reduce(intersect, lapply(phylum_data, rownames))
phylum_data <- lapply(phylum_data, function(m) m[common, , drop = FALSE])

################################################################################
# Top phyla across all streams
################################################################################

all_totals <- numeric()
for (nm in names(phylum_data)) {
  tots <- colSums(phylum_data[[nm]])
  for (t in names(tots)) {
    if (t %in% names(all_totals)) {
      all_totals[t] <- all_totals[t] + tots[t]
    } else {
      all_totals[t] <- tots[t]
    }
  }
}
all_totals <- all_totals[names(all_totals) != "Unassigned"]
top_phyla  <- names(sort(all_totals, decreasing = TRUE))[
  seq_len(min(N_TOP_PHYLA, length(all_totals)))]

################################################################################
# Build stacked data frame
################################################################################

prep_stack <- function(mat, method_name) {
  n_reads <- rowSums(mat)          # raw counts before normalization
  rel <- mat / n_reads
  df  <- as.data.frame(rel)
  
  # Collapse non-top to Other
  keep   <- intersect(top_phyla, colnames(df))
  others <- setdiff(colnames(df), c(keep, "Unassigned"))
  df2 <- df[, keep, drop = FALSE]
  df2$Other <- if (length(others) > 0) rowSums(df[, others, drop = FALSE]) else 0
  if ("Unassigned" %in% colnames(df)) df2$Unassigned <- df[, "Unassigned"]
  
  df2$Sample  <- rownames(df2)
  df2$Method  <- method_name
  df2$N_reads <- n_reads[rownames(df2)]   # total reads per sample
  
  pivot_longer(df2,
               cols      = -c(Sample, Method, N_reads),
               names_to  = "Phylum",
               values_to = "Abundance"
  )
}

stack_data <- do.call(rbind, lapply(names(phylum_data), function(nm) {
  prep_stack(phylum_data[[nm]], nm)
}))

# Join metadata
stack_data <- left_join(stack_data,
                        metadata[, c("Sample","Group","Treatment","Time")], by = "Sample")

# Factor ordering
stack_data$Sample <- factor(stack_data$Sample, levels = sample_order)
stack_data$Method <- factor(stack_data$Method,
                            levels = c("Amplicon","MEGAHIT Contigs","Read-Based"))
stack_data$Group  <- factor(stack_data$Group, levels = group_order)

taxon_order <- c(top_phyla, "Other", "Unassigned")
taxon_order <- taxon_order[taxon_order %in% unique(stack_data$Phylum)]
stack_data$Phylum <- factor(stack_data$Phylum, levels = rev(taxon_order))

################################################################################
# Colors
################################################################################

n_col <- length(taxon_order)
if (n_col <= 12) {
  phy_colors <- brewer.pal(max(3, n_col), "Set3")[seq_len(n_col)]
} else {
  phy_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_col)
}
names(phy_colors) <- taxon_order
phy_colors["Other"]      <- "grey70"
phy_colors["Unassigned"] <- "grey90"

################################################################################
# Plot
################################################################################

p <- ggplot(stack_data,
            aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", width = 0.85) +
  facet_grid(Method ~ Group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = phy_colors, name = "Phylum") +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0),
                     limits = c(0, 1.05)) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y      = element_text(size = 8),
    strip.text.x     = element_text(face = "bold", size = 10),
    strip.text.y     = element_text(face = "bold", size = 9),
    legend.position  = "right",
    legend.key.size  = unit(0.6, "cm"),
    legend.text      = element_text(size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.spacing.x    = unit(0.25, "cm"),
    panel.spacing.y    = unit(0.4, "cm")
  ) +
  labs(
    title    = "Phylum-Level Community Composition by Extraction Method",
    subtitle = "CT = Control, EW = Earthworm",
    x        = NULL,
    y        = "Relative Abundance"
  )

ggsave(file.path(output_dir, "three_stream_experimental_design_barplot.pdf"),
       p, width = 18, height = 11, dpi = 300)
ggsave(file.path(output_dir, "three_stream_experimental_design_barplot.png"),
       p, width = 18, height = 11, dpi = 300)

p <- ggplot(stack_data,
            aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", width = 0.85) +
  facet_grid(Method ~ Group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = phy_colors, name = "Phylum") +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0),
                     limits = c(0, 1.05)) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y      = element_text(size = 8),
    strip.text.x     = element_text(face = "bold", size = 10),
    strip.text.y     = element_text(face = "bold", size = 9),
    legend.position  = "right",
    legend.key.size  = unit(0.6, "cm"),
    legend.text      = element_text(size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.spacing.x    = unit(0.25, "cm"),
    panel.spacing.y    = unit(0.4, "cm")
  ) +
  labs(
    # title    = "",
    # subtitle = "",
    x        = NULL,
    y        = "Relative Abundance"
  )

ggsave(file.path(output_dir, "three_stream_experimental_design_barplot-notitle.pdf"),
       p, width = 18, height = 11, dpi = 300)
ggsave(file.path(output_dir, "three_stream_experimental_design_barplot-notitle.png"),
       p, width = 18, height = 11, dpi = 300)


cat("Saved: three_stream_experimental_design_barplot.pdf/png\n")
