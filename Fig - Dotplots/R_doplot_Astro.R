# ------------------------------------------------------------------------------
# Fig. 5d - Dot plot, astrocytes (mean expression and fraction of cells)
#
# Input:  Astrocytes.csv
# Output: R_doplot_Astro.svg
# ------------------------------------------------------------------------------

# Load libraries
library(ggplot2)
library(dplyr)
library(svglite)

# Read input data
df <- read.csv("Astrocytes.csv")

# Order genes by expression in APP (highest at top)
app_expression <- df %>%
  filter(Genotype == "APP") %>%
  arrange(desc(AverageExpression))

gene_order <- app_expression$Gene

df$Gene <- factor(df$Gene, levels = rev(gene_order))

# Set genotype order
df$Genotype <- factor(df$Genotype, levels = c("WT", "Ogg1", "APP", "AO"))

# Symmetric color scale limits
expr_range <- range(df$AverageExpression, na.rm = TRUE)
max_abs <- max(abs(expr_range))
symmetric_limit <- ceiling(max_abs)  # round up to whole number

# Number of genes and genotypes
num_genes <- length(levels(df$Gene))
num_genotypes <- length(levels(df$Genotype))

# Create dot plot (values are plotted as given, no z-scoring)
p <- ggplot(df, aes(x = Genotype, y = Gene)) +
  geom_point(aes(fill = AverageExpression, size = Fraction),
             shape = 21, color = "black", stroke = 0.7) +
  scale_fill_gradient2(
    name = "Normalized\n mean expression",
    low = "#2C5AA0",
    mid = "#f7f7f7",
    high = "#E86B00",
    midpoint = 0,
    limits = c(-symmetric_limit, symmetric_limit),
    breaks = c(-symmetric_limit, -symmetric_limit/2, 0, symmetric_limit/2, symmetric_limit),
    na.value = "red"
  ) +
  scale_size(
    name = "Fraction of cells (%)",
    range = c(2, 8),
    limits = c(0, 100),
    breaks = c(5, 25, 60)
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
    legend.direction = "vertical",
    legend.box.just = "left",
    legend.spacing.y = unit(0.5, "cm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(1.5, 1.5, 1.5, 1.5, "cm"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(x = NULL, y = NULL) +
  coord_fixed(ratio = 0.92, clip = "off")

# Legend order: color bar first, dot size second
p <- p + guides(
  fill = guide_colorbar(
    order = 1,
    barheight = unit(3, "cm"),
    barwidth = unit(0.5, "cm"),
    frame.colour = "black",
    frame.linewidth = 0.35,
    ticks.colour = "black",
    direction = "vertical",
    title.position = "top",
    title.hjust = 0.5
  ),
  size = guide_legend(
    order = 2,
    override.aes = list(fill = "white"),
    direction = "vertical",
    title.position = "top",
    title.hjust = 0.5
  )
)

# Show plot
print(p)

# Save plot
ggsave("R_doplot_Astro.svg", plot = p,
       width = 6,
       height = 5.5,
       bg = "white", device = "svg")
