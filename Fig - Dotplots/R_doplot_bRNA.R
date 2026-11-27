# ------------------------------------------------------------------------------
# Fig. 1i - Dot plot, bulk RNA-seq (z-scored mean expression per genotype)
#
# Input:  Average-FPKM_bRNA.csv
# Output: R_doplot_bRNA.svg
# ------------------------------------------------------------------------------

# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read input data
df <- read.csv("Average-FPKM_bRNA.csv")

# Reshape to long format
df_long <- df %>%
  pivot_longer(
    cols = c(AO, APP, OGG1, WT),
    names_to = "Genotype",
    values_to = "AverageExpression"
  ) %>%
  select(Gene, Genotype, AverageExpression)

# Select genes to plot
target_genes <- c("Cd48", "Cst7", "Ifi204", "Ifit1", "Irf7", "Nlrc5", "Nlrp1a", "Zbp1")

# Order genes by expression in APP (highest at top)
gene_order <- df %>%
  filter(Gene %in% target_genes) %>%
  arrange(desc(APP)) %>%
  pull(Gene)

df_filtered <- df_long %>%
  filter(Gene %in% target_genes) %>%
  mutate(Gene = factor(Gene, levels = rev(gene_order)))

# Calculate z-scores per gene
df_zscored <- df_filtered %>%
  group_by(Gene) %>%
  mutate(
    ZScore = scale(AverageExpression),
    ZScore = ifelse(is.na(ZScore), 0, ZScore)
  ) %>%
  ungroup() %>%
  filter(!is.na(AverageExpression))

# Set genotype order
df_zscored$Genotype <- factor(
  df_zscored$Genotype,
  levels = c("WT", "OGG1", "APP", "AO")
)

# Z-score range
zscore_range <- range(df_zscored$ZScore, na.rm = TRUE)

# Number of genes and genotypes
num_genes <- length(target_genes)
num_genotypes <- nlevels(df_zscored$Genotype)

# Create dot plot
p <- ggplot(df_zscored, aes(x = Genotype, y = Gene)) +
  geom_point(
    aes(fill = ZScore),
    size = 6,
    shape = 21,
    color = "black",
    stroke = 0.7
  ) +
  scale_fill_gradient2(
    name = "Normalized\n mean expression\n(z-score)",
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
      title.theme = element_text(size = 11, hjust = 0.5),
      direction = "vertical",
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
ggsave("R_doplot_bRNA.svg", plot = p,
       width = 6,
       height = 5.5,
       bg = "white", device = "svg")
