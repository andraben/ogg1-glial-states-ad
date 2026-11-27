# ------------------------------------------------------------------------------
# Fig. 6a - Concordance between ATAC-seq, bulk RNA-seq and snRNA-seq
#
# Input:  prevDEG_<celltype>_LogFCs.csv              (DEG-only concordance)
#         <genotype>_<celltype>_allgenes_LogFCs.csv  (genome-wide concordance)
# Output: Figure_1_spearman_rho_<approach>_<celltype>.png / .svg
#         Figure_3_median_abs_log2FC_<approach>_<celltype>.png / .svg
#
# Both approaches are run for all three cell types (astro, micro, oligo):
# 2 approaches x 3 cell types x 2 figures x 2 formats = 24 files
#
# Figure 1 shows Spearman's rho between assay pairs, Figure 3 the median
# |log2FC| per assay. The Deming slope (Figure 2) is still calculated but
# not plotted
#
# The snRNA-seq columns are named "<genotype>_scRNA" in the input files,
# only the figure labels say "snRNA"
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(svglite)

# Settings
OUTPUT_SIZE_IN <- 6  # size of the square canvas (inches)

offwhite <- "#F7F7F4"
# Colors of the Deming slope figure (not plotted)
div_low  <- "#2166AC"
div_high <- "#B35838"

# Figure 1 (Spearman's rho): diverging palette, used as continuous scale
RHO_PALETTE <- c("#2166AC", "#67A4C5", "#D1E1E8",
                 "#F7F7F7",
                 "#F1D0BC", "#D98765", "#B35838")

# Figure 3 (median |log2FC|): sequential palette, used as continuous scale
MAG_PALETTE <- c("#313695", "#4575B4", "#74ADD1",
                 "#ABD9C6", "#D0E5B5", "#E5E9A0", "#D6D88A")

# Luminance threshold (0-255) for the tile labels: black text on light tiles,
# white text on dark tiles (see tile_text_color() below)
TEXT_LUMINANCE_THRESHOLD <- 150

# Plot titles
FIG1_TITLE <- "Spearman Correlation (\u03c1)"
FIG3_TITLE <- "Median |log2FC|"

# Fixed font so the SVG text looks the same on every machine
FONT_FAMILY <- "Arial"

genotypes      <- c("Ogg1", "APP", "AO")
# Genotype labels as plotmath expressions (superscript for Ogg1). "-/-" is
# quoted so it is treated as text
genotype_labels <- c(
  Ogg1 = 'Ogg1^{"-/-"}',
  APP  = '"APP/PS1"',
  AO   = '"AO"'
)
celltypes      <- c("astro", "micro", "oligo")
# Assay pairs. The first assay of each pair is the x (reference) axis for the
# Deming slope
pair_list      <- list(c("ATAC", "bRNA"), c("bRNA", "scRNA"), c("ATAC", "scRNA"))
pair_order     <- map_chr(pair_list, ~ paste0(.x[1], "-", .x[2]))
# Order and labels of the assay pairs on the x-axis of Figure 1
pair_display_limits <- c("ATAC-scRNA", "ATAC-bRNA", "bRNA-scRNA")
pair_display_labels <- c("ATAC-snRNA", "ATAC-bRNA", "bRNA-snRNA")
assays_c_order <- c("ATAC", "bRNA", "scRNA")  # column order in Figure 3
# Assay labels on the x-axis of Figure 3
assay_display_labels <- c(ATAC = "ATAC", bRNA = "bRNA", scRNA = "snRNA")
RHO_STABLE     <- 0.3  # slope is not estimated if |rho| is below this value
# Color scale limits and breaks for rho and slope, set slightly above the
# observed range of the data. Breaks are given explicitly so that the scale
# limits are labeled in the legend
RHO_LIMIT      <- 0.7  # observed rho: -0.34 to 0.66
SLOPE_LIMIT    <- 7  # observed slope - 1: -6.53 to 4.04
RHO_BREAKS     <- seq(-RHO_LIMIT, RHO_LIMIT, by = 0.35)
SLOPE_BREAKS   <- seq(-SLOPE_LIMIT, SLOPE_LIMIT, by = 3.5)
# Color scale of Figure 3, shared by both approaches. The genome-wide values
# only reach about 0.27 and therefore stay at the lower end of the scale
MAG_LIMIT      <- 1.5
MAG_BREAKS     <- seq(0, MAG_LIMIT, by = 0.5)
BORDER_LWD     <- 0.5  # line width of tile and legend borders

# Font sizes (pt) and spacing (pt)
AXIS_FONT_SIZE     <- 13  # axis labels
LEGEND_FONT_SIZE   <- 13  # legend title and labels
HEADLINE_FONT_SIZE <- 14  # plot title
TILE_LABEL_SIZE    <- 5.2  # numbers inside the tiles (mm)
LEGEND_TITLE_GAP <- 16  # between legend title and color bar
HEADLINE_GAP     <- 6  # between plot title and plot
AXIS_LABEL_GAP   <- 2  # between axis labels and plot

# Input files use ";" as delimiter, "," as decimal mark and "#N/A" for
# missing values
eu_locale <- locale(decimal_mark = ",")
read_eu   <- function(path) {
  read_delim(path, delim = ";", locale = eu_locale,
             na = c("#N/A", "NA", ""), show_col_types = FALSE)
}

# Read prevDEG file (one file per cell type, all genotypes in one table)
load_prevdeg <- function(celltype) {
  df <- read_eu(paste0("prevDEG_", celltype, "_LogFCs.csv"))
  names(df)[1] <- "Gene"
  df
}

# Read allgenes files (one file per genotype and cell type) and join the three
# genotypes by gene. The snRNA-seq column (e.g. "AO_astro") is renamed to
# "<genotype>_scRNA"
load_allgenes <- function(celltype) {
  per_genotype <- map(genotypes, function(g) {
    path <- paste0(g, "_", celltype, "_allgenes_LogFCs.csv")
    df <- read_eu(path)
    names(df)[1] <- "Gene"
    sc_col_old <- paste0(g, "_", celltype)
    sc_col_new <- paste0(g, "_scRNA")
    df %>% rename(!!sc_col_new := all_of(sc_col_old))
  })
  reduce(per_genotype, full_join, by = "Gene")
}

# Deming (major-axis) regression slope
deming_slope <- function(x, y, delta = 1) {
  xbar <- mean(x); ybar <- mean(y)
  sxx <- mean((x - xbar)^2)
  syy <- mean((y - ybar)^2)
  sxy <- mean((x - xbar) * (y - ybar))
  (syy - delta * sxx + sqrt((syy - delta * sxx)^2 + 4 * delta * sxy^2)) / (2 * sxy)
}

# Black or white tile label, depending on the luminance of the tile color
# (works for sequential and diverging palettes)
tile_text_color <- function(values, limits, palette) {
  frac <- (values - limits[1]) / (limits[2] - limits[1])
  # NA (n.d.) is set to the midpoint
  frac <- ifelse(is.na(frac), 0.5, pmin(pmax(frac, 0), 1))
  rgb_mat <- grDevices::colorRamp(palette, space = "Lab")(frac)
  luminance <- 0.2126 * rgb_mat[, 1] + 0.7152 * rgb_mat[, 2] + 0.0722 * rgb_mat[, 3]
  ifelse(luminance > TEXT_LUMINANCE_THRESHOLD, "black", "white")
}

# Spearman's rho and Deming slope for one genotype and assay pair
compute_pair_stats <- function(df, genotype, a, b) {
  ca <- paste0(genotype, "_", a)
  cb <- paste0(genotype, "_", b)
  sub <- df %>% select(all_of(c(ca, cb))) %>% drop_na()
  rho <- suppressWarnings(cor(sub[[ca]], sub[[cb]], method = "spearman"))
  slope <- if (!is.na(rho) && abs(rho) >= RHO_STABLE) {
    deming_slope(sub[[ca]], sub[[cb]])
  } else {
    NA_real_
  }
  tibble(genotype = genotype, pair_label = paste0(a, "-", b), rho = rho, slope = slope)
}

# Plot theme. No axis titles and no panel border (the tile borders form the
# frame)
base_theme <- theme_minimal(base_size = AXIS_FONT_SIZE) +
  theme(
    text             = element_text(family = FONT_FAMILY),
    axis.text        = element_text(color = "black", face = "plain", size = AXIS_FONT_SIZE,
                                    margin = margin(t = AXIS_LABEL_GAP)),
    axis.text.y      = element_text(color = "black", face = "plain", size = AXIS_FONT_SIZE,
                                    margin = margin(r = AXIS_LABEL_GAP)),
    axis.title       = element_blank(),
    legend.title     = element_text(color = "black", face = "plain", size = LEGEND_FONT_SIZE),
    legend.text      = element_text(color = "black", face = "plain", size = LEGEND_FONT_SIZE),
    panel.grid       = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    plot.title       = element_text(color = "black", face = "plain", size = HEADLINE_FONT_SIZE,
                                    hjust = 0.5, margin = margin(b = HEADLINE_GAP))
  )

# Color bar with black frame and ticks
black_colorbar <- guide_colorbar(
  frame.colour = "black", frame.linewidth = BORDER_LWD,
  ticks.colour = "black", ticks.linewidth = BORDER_LWD,
  title.theme = element_text(color = "black", size = LEGEND_FONT_SIZE,
                             margin = margin(b = LEGEND_TITLE_GAP))
)

# Save plot as PNG and SVG
save_fig <- function(plot, name, size = OUTPUT_SIZE_IN) {
  ggsave(paste0(name, ".png"), plot, width = size, height = size, dpi = 300, bg = "white")
  ggsave(paste0(name, ".svg"), plot, width = size, height = size, device = svglite, bg = "white")
}

# Create and save Figure 1 and Figure 3 for one approach and cell type
build_figs <- function(df, approach, celltype) {

  results <- map_dfr(genotypes, function(g) {
    map_dfr(pair_list, function(p) compute_pair_stats(df, g, p[1], p[2]))
  }) %>%
    mutate(
      genotype   = factor(genotype, levels = genotypes),
      pair_label = factor(pair_label, levels = pair_order)
    )

  rho_df   <- results %>%
    select(genotype, pair_label, rho) %>%
    mutate(text_color = tile_text_color(rho, c(-RHO_LIMIT, RHO_LIMIT), RHO_PALETTE))

  mag_df <- map_dfr(genotypes, function(g) {
    cols <- paste0(g, "_", assays_c_order)
    sub  <- df %>% select(all_of(cols)) %>% drop_na()
    meds <- vapply(sub, function(v) median(abs(v)), numeric(1))
    tibble(genotype = g, assay = assays_c_order, med_abs_logfc = meds)
  }) %>%
    mutate(
      genotype = factor(genotype, levels = genotypes),
      assay    = factor(assay, levels = assays_c_order)
    )

  # y-axis: Ogg1, APP/PS1, AO from top to bottom, labels parsed as plotmath
  y_scale <- scale_y_discrete(limits = rev(genotypes),
                              labels = function(b) parse(text = genotype_labels[b]))
  # x-axis of Figure 1: assay pairs, labels tilted to avoid overlap
  x_scale <- scale_x_discrete(limits = pair_display_limits, labels = pair_display_labels)
  tilted_x <- theme(axis.text.x = element_text(angle = 33, hjust = 1, vjust = 1))

  # Figure 1: Spearman's rho
  p1 <- ggplot(rho_df, aes(x = pair_label, y = genotype, fill = rho)) +
    geom_tile(color = "black", linewidth = BORDER_LWD) +
    geom_text(aes(label = ifelse(is.na(rho), "n.d.", sprintf("%.2f", rho)), color = text_color),
              size = TILE_LABEL_SIZE, family = FONT_FAMILY) +
    scale_color_identity() +
    scale_fill_gradientn(
      colors = RHO_PALETTE,
      limits = c(-RHO_LIMIT, RHO_LIMIT), breaks = RHO_BREAKS, na.value = "grey90",
      name = "\u03c1", guide = black_colorbar
    ) +
    x_scale +
    y_scale +
    coord_fixed() +
    base_theme +
    tilted_x +
    labs(title = FIG1_TITLE)

  # Figure 3: median |log2FC|
  mag_df <- mag_df %>%
    mutate(text_color = tile_text_color(med_abs_logfc, c(0, MAG_LIMIT), MAG_PALETTE))

  p3 <- ggplot(mag_df, aes(x = assay, y = genotype, fill = med_abs_logfc)) +
    geom_tile(color = "black", linewidth = BORDER_LWD) +
    geom_text(aes(label = sprintf("%.2f", med_abs_logfc), color = text_color), size = TILE_LABEL_SIZE, family = FONT_FAMILY) +
    scale_color_identity() +
    scale_fill_gradientn(colors = MAG_PALETTE,
                         limits = c(0, MAG_LIMIT), breaks = MAG_BREAKS,
                         name = "Median\n|log2FC|", guide = black_colorbar) +
    scale_x_discrete(labels = assay_display_labels) +
    y_scale +
    coord_fixed() +
    base_theme +
    labs(title = FIG3_TITLE)

  save_fig(p1, paste0("Figure_1_spearman_rho_", approach, "_", celltype))
  save_fig(p3, paste0("Figure_3_median_abs_log2FC_", approach, "_", celltype))

  message(sprintf(" [%s / %s] n_genes (after per-genotype pair drop_na varies): %d rows in merged table",
                  approach, celltype, nrow(df)))
}

# Run both approaches for all three cell types
for (ct in celltypes) {
  build_figs(load_prevdeg(ct),   approach = "prevDEG",  celltype = ct)
  build_figs(load_allgenes(ct),  approach = "allgenes",  celltype = ct)
}

message("Done. Wrote 24 files (2 approaches x 3 celltypes x 2 figures x [png+svg]).")
