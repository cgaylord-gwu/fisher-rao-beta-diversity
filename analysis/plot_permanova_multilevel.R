#!/usr/bin/env Rscript
################################################################################
# Cross-level PERMANOVA summary figures
#
# Reads permanova_fr_bc_multilevel.csv and produces:
#
#   1. permanova_crosslevel_pvalues.png/pdf
#      Tile heatmap: rows = Level × Stream, cols = Effect × Metric
#      Color = -log10(p), with significance threshold overlaid
#
#   2. permanova_crosslevel_r2.png/pdf
#      Dot + segment plot: R² by level, faceted by Effect,
#      colored by Metric, shaped by Stream
#
# Both figures are publication-quality at 300 dpi.
################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)

BASE <- "."

output_dir <- file.path(BASE, "figures")
csv_path   <- file.path(output_dir, "permanova_fr_bc_multilevel.csv")

df <- read.csv(csv_path, stringsAsFactors = FALSE)

# Factor ordering
df$Level  <- factor(df$Level,  levels = c("Phylum", "Class", "Family", "Genus"))
df$Stream <- factor(df$Stream, levels = c("Amplicon", "MEGAHIT Contigs", "Read Based"))
df$Metric <- factor(df$Metric, levels = c("Fisher-Rao", "Bray-Curtis"))

################################################################################
# Figure 1: p-value heatmap
################################################################################

# Reshape to long form: one row per Level × Stream × Metric × Effect
p_long <- df %>%
  select(Level, Stream, Metric,
         Treatment_p, Time_p, Interaction_p) %>%
  pivot_longer(
    cols      = c(Treatment_p, Time_p, Interaction_p),
    names_to  = "Effect",
    values_to = "p_value"
  ) %>%
  mutate(
    Effect = recode(Effect,
                    Treatment_p   = "Treatment",
                    Time_p        = "Time",
                    Interaction_p = "Treatment × Time"),
    Effect      = factor(Effect, levels = c("Treatment", "Time", "Treatment × Time")),
    neg_log10_p = -log10(p_value),
    sig_label   = case_when(
      p_value <= 0.001 ~ "***",
      p_value <= 0.01  ~ "**",
      p_value <= 0.05  ~ "*",
      p_value <= 0.10  ~ "†",
      TRUE             ~ ""
    ),
    # Row label: Level + Stream
    Row = paste0(Level, "\n", Stream)
  )

# Ordered row factor: Phylum→Genus top-to-bottom, streams within level
row_order <- p_long %>%
  arrange(Level, Stream) %>%
  distinct(Row) %>%
  pull(Row)
p_long$Row <- factor(p_long$Row, levels = rev(row_order))  # rev so Phylum on top

# Column label: Effect + Metric side by side
p_long$Col <- paste0(Effect = as.character(p_long$Effect), "\n", p_long$Metric)
col_order <- c()
for (eff in c("Treatment", "Time", "Treatment × Time")) {
  for (met in c("Fisher-Rao", "Bray-Curtis")) {
    col_order <- c(col_order, paste0(eff, "\n", met))
  }
}
p_long$Col <- factor(p_long$Col, levels = col_order)

p1 <- ggplot(p_long, aes(x = Col, y = Row, fill = neg_log10_p)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sig_label), size = 4.5, fontface = "bold",
            color = "white", vjust = 0.5) +
  scale_fill_gradientn(
    colours = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
    values  = scales::rescale(c(0, 1, 1.301, 2, 3.5)),
    limits  = c(0, NA),
    name    = expression(-log[10](p))
  ) +
  # Vertical separator between effect groups
  geom_vline(xintercept = c(2.5, 4.5), color = "grey50", linewidth = 0.8) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 0, hjust = 0.5, size = 9,
                                    lineheight = 1.1),
    axis.text.y      = element_text(size = 9, lineheight = 1.1),
    axis.title       = element_blank(),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "grey70"),
    legend.position  = "right",
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 10, color = "grey40")
  ) +
  labs(
    title    = "PERMANOVA: Treatment × Time Effects Across Taxonomic Levels",
    subtitle = "† p<0.10   * p<0.05   ** p<0.01   *** p<0.001"
  )

out1 <- file.path(output_dir, "permanova_crosslevel_pvalues")
ggsave(paste0(out1, ".png"), p1, width = 11, height = 7, dpi = 300)
ggsave(paste0(out1, ".pdf"), p1, width = 11, height = 7)

p1 <- ggplot(p_long, aes(x = Col, y = Row, fill = neg_log10_p)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sig_label), size = 4.5, fontface = "bold",
            color = "white", vjust = 0.5) +
  scale_fill_gradientn(
    colours = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
    values  = scales::rescale(c(0, 1, 1.301, 2, 3.5)),
    limits  = c(0, NA),
    name    = expression(-log[10](p))
  ) +
  # Vertical separator between effect groups
  geom_vline(xintercept = c(2.5, 4.5), color = "grey50", linewidth = 0.8) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 0, hjust = 0.5, size = 9,
                                    lineheight = 1.1),
    axis.text.y      = element_text(size = 9, lineheight = 1.1),
    axis.title       = element_blank(),
    panel.grid       = element_blank(),
    panel.border     = element_rect(color = "grey70"),
    legend.position  = "right",
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 10, color = "grey40")
  ) +
  labs(
    # title    = "PERMANOVA: Treatment × Time Effects Across Taxonomic Levels",
    # subtitle = "† p<0.10   * p<0.05   ** p<0.01   *** p<0.001"
  )

out1 <- file.path(output_dir, "permanova_crosslevel_pvalues-notitle")
ggsave(paste0(out1, ".png"), p1, width = 11, height = 7, dpi = 300)
ggsave(paste0(out1, ".pdf"), p1, width = 11, height = 7)


cat("Saved: permanova_crosslevel_pvalues.png/pdf\n")

################################################################################
# Figure 2: R² dot plot — how much variance explained, by level
################################################################################

r2_long <- df %>%
  select(Level, Stream, Metric,
         Treatment_R2, Time_R2, Interaction_R2) %>%
  pivot_longer(
    cols      = c(Treatment_R2, Time_R2, Interaction_R2),
    names_to  = "Effect",
    values_to = "R2"
  ) %>%
  mutate(
    Effect = recode(Effect,
                    Treatment_R2   = "Treatment",
                    Time_R2        = "Time",
                    Interaction_R2 = "Treatment × Time"),
    Effect = factor(Effect, levels = c("Treatment", "Time", "Treatment × Time"))
  )

# Join p-values for significance shading
p_vals <- p_long %>%
  select(Level, Stream, Metric, Effect, p_value) %>%
  mutate(Effect = as.character(Effect))
r2_long <- r2_long %>%
  mutate(Effect = as.character(Effect)) %>%
  left_join(p_vals, by = c("Level", "Stream", "Metric", "Effect")) %>%
  mutate(
    Effect    = factor(Effect, levels = c("Treatment", "Time", "Treatment × Time")),
    Metric    = factor(Metric, levels = c("Fisher-Rao", "Bray-Curtis")),
    Sig       = p_value < 0.05,
    alpha_val = ifelse(Sig, 1.0, 0.35)
  )

p2 <- ggplot(r2_long,
             aes(x     = Level,
                 y     = R2,
                 color = Metric,
                 shape = Stream,
                 alpha = alpha_val,
                 group = interaction(Stream, Metric))) +
  geom_line(aes(group = interaction(Stream, Metric)),
            linewidth = 0.5, linetype = "dashed",
            alpha = 0.25) +   # fixed low alpha for all lines
  geom_point(aes(alpha = alpha_val), size = 3.5) +
  scale_color_manual(values = c("Fisher-Rao"  = "#E41A1C",
                                "Bray-Curtis" = "#377EB8")) +
  scale_shape_manual(values = c("Amplicon"        = 16,
                                "MEGAHIT Contigs" = 17,
                                "Read Based"      = 15)) +
  scale_alpha_identity() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_wrap(~ Effect, ncol = 3, scales = "free_y") +
  theme_bw(base_size = 11) +
  theme(
    strip.text      = element_text(face = "bold", size = 11),
    axis.text.x     = element_text(angle = 30, hjust = 1),
    legend.position = "bottom",
    legend.box      = "vertical",
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 10, color = "grey40")
  ) +
  labs(
    title    = "PERMANOVA R² by Taxonomic Level",
    subtitle = "Faded points: p ≥ 0.05",
    x        = "Taxonomic Level",
    y        = "R²",
    color    = "Metric",
    shape    = "Stream"
  )

out2 <- file.path(output_dir, "permanova_crosslevel_r2")
ggsave(paste0(out2, ".png"), p2, width = 12, height = 6, dpi = 300)
ggsave(paste0(out2, ".pdf"), p2, width = 12, height = 6)
cat("Saved: permanova_crosslevel_r2.png/pdf\n")

p2 <- ggplot(r2_long,
             aes(x     = Level,
                 y     = R2,
                 color = Metric,
                 shape = Stream,
                 alpha = alpha_val,
                 group = interaction(Stream, Metric))) +
  geom_line(aes(group = interaction(Stream, Metric)),
            linewidth = 0.5, linetype = "dashed",
            alpha = 0.25) +   # fixed low alpha for all lines
  geom_point(aes(alpha = alpha_val), size = 3.5) +
  scale_color_manual(values = c("Fisher-Rao"  = "#E41A1C",
                                "Bray-Curtis" = "#377EB8")) +
  scale_shape_manual(values = c("Amplicon"        = 16,
                                "MEGAHIT Contigs" = 17,
                                "Read Based"      = 15)) +
  scale_alpha_identity() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_wrap(~ Effect, ncol = 3, scales = "free_y") +
  theme_bw(base_size = 11) +
  theme(
    strip.text      = element_text(face = "bold", size = 11),
    axis.text.x     = element_text(angle = 30, hjust = 1),
    legend.position = "bottom",
    legend.box      = "vertical",
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 10, color = "grey40")
  ) +
  labs(
    # title    = "PERMANOVA R² by Taxonomic Level",
    # subtitle = "Faded points: p ≥ 0.05",
    x        = "Taxonomic Level",
    y        = "R²",
    color    = "Metric",
    shape    = "Stream"
  )

out2 <- file.path(output_dir, "permanova_crosslevel_r2-notitle")
ggsave(paste0(out2, ".png"), p2, width = 12, height = 6, dpi = 300)
ggsave(paste0(out2, ".pdf"), p2, width = 12, height = 6)
cat("Saved: permanova_crosslevel_r2-notitle.png/pdf\n")

cat("\n=== Done ===\n")
