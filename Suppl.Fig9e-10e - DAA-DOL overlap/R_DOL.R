# ------------------------------------------------------------------------------
# Suppl. Fig. 10e - Heatmap of DOL gene expression in oligodendrocytes
#
# Input:  03_Oligo.csv
# Output: DOL.svg
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# Read input data
df <- read.csv("03_Oligo.csv") %>%
  mutate(
    # Rename APP to APP/PS1 for display
    Genotype = case_when(
      Genotype == "APP" ~ "APP/PS1",
      TRUE ~ Genotype
    ),
    Condition = Genotype
  )

# Average duplicated genes
df_unique <- df %>%
  group_by(Gene, Condition) %>%
  summarise(AverageExpression = mean(AverageExpression, na.rm = TRUE), .groups = "drop")

# Column order for the heatmap
condition_order <- c("WT", "Ogg1", "APP/PS1", "AO")

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

# Order genes by hierarchical clustering (values used as given, no z-scoring)
row_cluster <- hclust(dist(expression_matrix))
gene_order <- row_cluster$order

expression_matrix_ordered <- expression_matrix[gene_order, ]

# Genes shown in bold
bold_genes <- c("Taf13", "Cd9", "Lap3", "Emb", "Rtn4", "Apod", "Pde4b", "Sox2ot", "Ptgds", "Erbin", "Bfsp2", "C4b", "Rcan2", "Stmn2")

row_fontface <- ifelse(rownames(expression_matrix_ordered) %in% bold_genes, "bold", "plain")

# Annotation data
annotation_df <- data.frame(
  Genotype = c("WT", "Ogg1", "APP/PS1", "AO")
)
rownames(annotation_df) <- colnames(expression_matrix_ordered)

# Set legend order
annotation_df$Genotype <- factor(annotation_df$Genotype,
                                 levels = c("WT", "Ogg1", "APP/PS1", "AO"))

# Annotation colors
annotation_colors <- list(
  Genotype = c("WT" = "#5e5e5e", "Ogg1" = "#543c6c", "APP/PS1" = "#056475", "AO" = "#F57A0C")
)

# Color scale covering the data range
expr_range <- range(expression_matrix_ordered)
color_breaks <- seq(expr_range[1], expr_range[2], length.out = 100)
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
    )
  )
)

# Create heatmap
cell_width <- unit(0.4, "inches")
cell_height <- unit(0.15, "inches")

ht <- Heatmap(
  expression_matrix_ordered,
  name = "Normalized mean expression",
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
  width = ncol(expression_matrix_ordered) * cell_width,
  height = nrow(expression_matrix_ordered) * cell_height,
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
svg("DOL.svg", width = 5, height = 5, bg = "white")
draw(ht, heatmap_legend_side = "right")
dev.off()
