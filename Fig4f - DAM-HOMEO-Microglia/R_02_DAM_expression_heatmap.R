# ------------------------------------------------------------------------------
# Fig. 4f - Heatmap of DAM gene expression in microglia
#
# Input:  DAM_expression_results.csv
# Output: DAM_expression_heatmap.svg
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# Read input data and reshape to long format
df <- read.csv("DAM_expression_results.csv") %>%
  dplyr::select(-ENS_ID) %>%
  pivot_longer(
    cols = -Gene_Name,
    names_to = "Condition",
    values_to = "Expression"
  ) %>%
  # Split Condition into Genotype and Cluster
  separate(Condition, into = c("Genotype", "Cluster"), sep = "_") %>%
  rename(Gene = Gene_Name) %>%
  # Remove missing values
  filter(!is.na(Expression))

# Average duplicated genes
df_unique <- df %>%
  group_by(Gene, Condition = paste(Genotype, Cluster, sep = "_")) %>%
  summarise(AverageExpression = mean(Expression, na.rm = TRUE), .groups = "drop")

# Column order for the heatmap
condition_order <- c("WT_Microglia", "Ogg1_Microglia", "APP_Microglia", "AO_Microglia", "APP_Act1", "AO_Act1")

# Reshape to wide format
expression_wide <- df_unique %>%
  pivot_wider(
    names_from = Condition,
    values_from = AverageExpression,
    values_fill = 0
  )

# Convert to matrix
expression_matrix <- expression_wide %>%
  dplyr::select(Gene, all_of(condition_order)) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

# Calculate z-scores per gene
scaled_matrix <- t(scale(t(expression_matrix)))
scaled_matrix[is.na(scaled_matrix)] <- 0

# Order genes by hierarchical clustering
row_cluster <- hclust(dist(scaled_matrix))
gene_order <- row_cluster$order

scaled_matrix <- scaled_matrix[gene_order, ]

# Genes shown in bold
bold_genes <- c()

row_fontface <- ifelse(rownames(scaled_matrix) %in% bold_genes, "bold", "plain")

# Annotation data
annotation_df <- data.frame(
  Genotype = c("WT", "Ogg1", "APP/PS1", "AO", "APP/PS1", "AO"),
  Cluster = c("Microglia", "Microglia", "Microglia", "Microglia", "Act1", "Act1")
)
rownames(annotation_df) <- colnames(scaled_matrix)

# Set legend order
annotation_df$Genotype <- factor(annotation_df$Genotype,
                                 levels = c("WT", "Ogg1", "APP/PS1", "AO"))
annotation_df$Cluster <- factor(annotation_df$Cluster,
                                levels = c("Microglia", "Act1"))

# Annotation colors
annotation_colors <- list(
  Genotype = c("WT" = "#5e5e5e", "Ogg1" = "#543c6c", "APP/PS1" = "#056475", "AO" = "#F57A0C"),
  Cluster = c(Microglia = "#5DA8A2", Act1 = "#A85D89")
)

# Color scale, symmetric around zero
max_abs_z <- max(abs(scaled_matrix))
color_breaks <- seq(-max_abs_z, max_abs_z, length.out = 100)
col_fun <- colorRamp2(
  breaks = color_breaks,
  colors = colorRampPalette(c("#6A3D73", "#F2F1EE", "#E6B450"))(100)
)

# Create column annotation
ha_col <- HeatmapAnnotation(
  df = annotation_df,
  col = annotation_colors,
  annotation_name_side = "left",
  show_legend = TRUE,
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Genotype = list(
      title = "Genotype",
      title_gp = gpar(fontsize = 12, fontface = "plain"),
      labels_gp = gpar(fontsize = 11, fontface = "plain"),
      at = c("WT", "Ogg1", "APP/PS1", "AO")
    ),
    Cluster = list(
      title = "Cluster",
      title_gp = gpar(fontsize = 12, fontface = "plain"),
      labels_gp = gpar(fontsize = 11, fontface = "plain"),
      at = c("Microglia", "Act1")
    )
  )
)

# Create heatmap
cell_width <- unit(0.4, "inches")
cell_height <- unit(0.15, "inches")

ht <- Heatmap(
  scaled_matrix,
  name = "Normalized mean expression (z-score)",
  col = col_fun,
  border = TRUE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 11, fontface = row_fontface),
  row_title = NULL,
  show_column_names = FALSE,
  top_annotation = ha_col,
  row_gap = unit(2, "mm"),
  width = ncol(scaled_matrix) * cell_width,
  height = nrow(scaled_matrix) * cell_height,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 11, fontface = "plain"),
    labels_gp = gpar(fontsize = 8),
    legend_height = unit(4, "cm"),
    legend_width = unit(1, "cm"),
    border = "black",
    title_position = "leftcenter-rot"
  )
)

# Save heatmap
svg("DAM_expression_heatmap.svg", width = 6, height = 5, bg = "white")
draw(ht, heatmap_legend_side = "right")
dev.off()
