# ------------------------------------------------------------------------------
# Suppl. Fig. 5f - Dot plot, ATAC-seq (z-scored mean enrichment per genotype)
#
# Input:  ATAC_AvgExp.csv
# Output: R_doplot_ATAC.svg
#         ATAC_Averaged_Peaks_Per_Gene.csv
#         ATAC_ZScored_Average_Expression.csv
# ------------------------------------------------------------------------------

# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read input data
df <- read.csv("ATAC_AvgExp.csv")

# Average peak values per gene and genotype
df_avg <- df %>%
  pivot_longer(
    cols = c(AO, APP, OGG1, WT),
    names_to = "Genotype",
    values_to = "Value"
  ) %>%
  group_by(Gene, Genotype) %>%
  summarise(AverageExpression = mean(Value, na.rm = TRUE)) %>%
  ungroup()

# Save averaged values
write.csv(df_avg, "ATAC_Averaged_Peaks_Per_Gene.csv", row.names = FALSE)

# Select genes to plot
target_genes <- c("Alk", "Ppp2r2b", "Cbln2", "Firre", "Fmr1", "Hepacam", "Ntrk2", "Nwd2")
df_filtered <- df_avg %>%
  filter(Gene %in% target_genes)

# Calculate z-scores per gene
df_zscored <- df_filtered %>%
  group_by(Gene) %>%
  mutate(
    ZScore = scale(AverageExpression),
    ZScore = ifelse(is.na(ZScore), 0, ZScore)
  ) %>%
  ungroup() %>%
  filter(!is.na(AverageExpression))

# Order genes by z-score in APP (highest at top)
gene_order <- df_zscored %>%
  filter(Genotype == "APP") %>%
  arrange(desc(ZScore)) %>%
  pull(Gene) %>%
  unique()

# Set gene and genotype order
df_zscored$Gene <- factor(df_zscored$Gene, levels = rev(gene_order))
df_zscored$Genotype <- factor(
  df_zscored$Genotype,
  levels = c("WT", "OGG1", "APP", "AO")
)

# Create dot plot
p <- ggplot(df_zscored, aes(x = Genotype, y = Gene)) +
  geom_point(
    aes(fill = ZScore),
    size = 7,
    shape = 21,
    color = "black",
    stroke = 0.7
  ) +
  scale_fill_gradient2(
    name = "Normalized\n mean enrichment\n(z-score)",
    low = "#2C5AA0",
    mid = "#f7f7f7",
    high = "#E86B00",
    midpoint = 0,
    limits = c(-2, 2),
    na.value = "red",
    guide = guide_colorbar(
      barheight = unit(3, "cm"),
      barwidth = unit(0.5, "cm"),
      frame.colour = "black",
      frame.linewidth = 0.35,
      ticks.colour = "black",
      direction = "vertical",
      title.theme = element_text(size = 11, hjust = 0.5),
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.15, 0.15))) +
  scale_y_discrete(expand = expansion(mult = c(0.1, 0.1))) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                               size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.42),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 10, hjust = 0.5),
    legend.text = element_text(size = 9, hjust = 0.5),
    plot.margin = margin(1.5, 1.5, 1.5, 1.5, "cm"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(x = NULL, y = NULL) +
  coord_fixed(ratio = 0.92, clip = "off")

# Show plot
print(p)

# Save plot
ggsave("R_doplot_ATAC.svg", plot = p,
       width = 6,
       height = 5.5,
       bg = "white", device = "svg")

# Save z-scores in wide format
df_zscored_wide <- df_zscored %>%
  select(Gene, Genotype, ZScore) %>%
  pivot_wider(
    names_from = Genotype,
    values_from = ZScore
  )

write.csv(df_zscored_wide, "ATAC_ZScored_Average_Expression.csv", row.names = FALSE)
