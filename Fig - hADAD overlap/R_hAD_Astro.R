# ------------------------------------------------------------------------------
# Fig. 5f - Heatmap, overlap with human AD and ADAD data, astrocytes
#
# Input:  Astrocytes_R.csv
# Output: hADAD_Astro.svg
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(reshape2)

# Read input data. Human logFC values (nebula, natural log) are converted to
# log2FC to match the mouse values (limma-voom, log2)
data <- read.csv("Astrocytes_R.csv") %>%
  mutate(across(c(sAD, hADAD), ~ .x / log(2)))

# Column order for the heatmap
logFC_columns <- c("sAD", "hADAD", "APPvsWT", "AOvsWT", "Ogg1vsWT")
pval_columns <- c("sAD_pval", "hADAD_pval","APPvsWT_pval", "AOvsWT_pval", "Ogg1vsWT_pval")

# Symmetric color scale limits (rounded up to 0.5)
fc_limit <- ceiling(max(abs(data[, logFC_columns]), na.rm = TRUE) * 2) / 2

# Order genes by APPvsWT log2FC (highest at top)
sorted_genes <- data %>%
  arrange(desc(APPvsWT)) %>%
  pull(Gene)

sorted_genes <- rev(sorted_genes)

# Position of each gene in the heatmap
gene_order <- data.frame(
  Gene = sorted_genes,
  Order = 1:length(sorted_genes)
)

# Genes shown in bold
bold_genes <- c("COL16A1", "PAQR4", "RLBP1", "IRAK4", "CHKB", "PPP1R26", "MAPK12", "CMTM5", "DERA", "PXMP4", "CRLF1", "STK32C", "COL25A1", "DOCK10", "CNIH3", "CACNA1H", "IL16", "FRMD5", "CDKL1", "TMEM266", "ZNHIT3", "FILIP1", "PTPRE", "SH3PXD2A", "AJAP1", "MMP24"
)

# Reshape to long format and add p-values
plot_data <- data %>%
  left_join(gene_order, by = "Gene") %>%
  dplyr::select(Gene, Order, all_of(logFC_columns)) %>%
  melt(id.vars = c("Gene", "Order"), variable.name = "Comparison", value.name = "logFC") %>%
  # Add p-values
  left_join(
    data %>%
      dplyr::select(Gene, all_of(pval_columns)) %>%
      melt(id.vars = "Gene", variable.name = "pval_var", value.name = "pval") %>%
      mutate(Comparison = str_remove(pval_var, "_pval")) %>%
      dplyr::select(-pval_var),
    by = c("Gene", "Comparison")
  ) %>%
  # Significance stars
  mutate(
    pval_stars = case_when(
      pval < 0.001 ~ "***",
      pval < 0.01 ~ "**",
      pval < 0.05 ~ "*",
      TRUE ~ ""
    ),
    is_significant = pval < 0.05
  ) %>%
  # Set column and gene order
  mutate(
    Comparison = factor(Comparison, levels = logFC_columns),
    Gene = factor(Gene, levels = sorted_genes)
  ) %>%
  arrange(Order, Comparison)

# Column labels
comparison_labels <- c("hAD", "hADAD", "mAPP/PS1", "mAO", "mOgg1-/-")
levels(plot_data$Comparison) <- comparison_labels

# Number of genes and comparisons
n_genes <- length(unique(plot_data$Gene))
n_comparisons <- length(unique(plot_data$Comparison))

# Cell size (inches)
cell_width <- 0.08
cell_height <- 0.08

# Plot size from number of cells
plot_width <- n_comparisons * cell_width + 2.2
plot_height <- n_genes * cell_height + 0.8

# Gene labels
gene_labels <- data.frame(
  Gene = sorted_genes,
  y_pos = 1:length(sorted_genes),
  is_bold = sorted_genes %in% bold_genes
)

# Create heatmap
p <- ggplot(plot_data, aes(x = Comparison, y = Gene)) +
  geom_tile(aes(fill = logFC), color = "white", linewidth = 0.1) +
  # Significance stars
  geom_text(aes(label = pval_stars, fontface = ifelse(is_significant, "bold", "plain")),
            size = 2.6, color = "black", vjust = 0.8) +
  # Gene labels left of the heatmap
  geom_text(data = gene_labels,
            aes(x = 0.4, y = y_pos, label = Gene,
                fontface = ifelse(is_bold, "bold", "plain")),
            hjust = 1,
            size = 2.2,
            color = "black",
            inherit.aes = FALSE) +
  scale_fill_gradient2(
    low = "#377B9F",
    mid = "#f5f5f5",
    high = "#d6604d",
    midpoint = 0,
    limits = c(-fc_limit, fc_limit),
    name = expression(log[2]*FC),
    guide = guide_colorbar(
      barwidth = 0.5,
      barheight = 5,
      title.position = "top",
      title.hjust = 0.5,
      frame.colour = "black",
      frame.linewidth = 0.2,
      ticks.colour = "black",
      title.theme = element_text(size = 7),
      label.theme = element_text(size = 6)
    )
  ) +
  scale_x_discrete(
    position = "top",
    expand = expansion(0)
  ) +
  scale_y_discrete(
    expand = expansion(0),
    labels = NULL
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    aspect.ratio = NULL,
    axis.text.x.top = element_text(
      angle = 35,
      hjust = 0,
      vjust = 0.1,
      size = 7,
      color = "black",
      margin = margin(b = 0, unit = "pt")
    ),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.2
    ),
    panel.background = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 8, color = "black"),
    legend.text = element_text(size = 6, color = "black"),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.key.size = unit(0.4, "cm"),
    plot.margin = margin(5, 5, 5, 60, unit = "pt"),  # room for gene labels
    panel.spacing = unit(0, "lines"),
    plot.background = element_blank()
  ) +
  coord_cartesian(
    expand = FALSE,
    xlim = c(0.5, n_comparisons + 0.5),
    clip = "off"  # allow gene labels outside the panel
  )

# Show plot
print(p)

# Save plot
ggsave("hADAD_Astro.svg", plot = p, width = plot_width, height = plot_height,
       units = "in", dpi = 1200, bg = "white")
