# ------------------------------------------------------------------------------
# Fig. 2g - Overlap of bulk RNA-seq and ATAC-seq gene lists per genotype
#
# Input:  Ogg1_UP-DOWN_bRNA-epigenetics.csv
#         APP_UP-DOWN_bRNA-epigenetics.csv
#         AO_UP-DOWN_bRNA-epigenetics.csv
# Output: ATAC_bRNA_intersection_all_genotypes.pdf / .png / .svg
#         ATAC_bRNA_overlapping_genes.xlsx
# ------------------------------------------------------------------------------

# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(openxlsx)

# Input files (one per genotype)
input_files <- list(
  "Ogg1"    = "Ogg1_UP-DOWN_bRNA-epigenetics.csv",
  "APP/PS1" = "APP_UP-DOWN_bRNA-epigenetics.csv",
  "AO"      = "AO_UP-DOWN_bRNA-epigenetics.csv"
)

# Helper functions
clean_gene_column <- function(column_data) {
  genes <- column_data[!is.na(column_data) & column_data != ""]
  genes <- genes[!genes %in% c("bRNA_UP", "bRNA_DOWN", "ATAC_UP", "ATAC_DOWN",
                               "H3K9ac_UP", "H3K9ac_DOWN", "H3K27ac_UP", "H3K27ac_DOWN")]
  return(genes)
}

count_intersect <- function(set1, set2) {
  sum(set1 %in% set2)
}

# Overlap between bulk RNA-seq and ATAC-seq gene lists per genotype
results    <- data.frame()
gene_lists <- list()  # gene names per genotype and category

for (geno in names(input_files)) {

  dat <- read.csv(input_files[[geno]], stringsAsFactors = FALSE)

  bRNA_UP   <- clean_gene_column(dat$bRNA_UP)
  bRNA_DOWN <- clean_gene_column(dat$bRNA_DOWN)
  ATAC_UP   <- clean_gene_column(dat$ATAC_UP)
  ATAC_DOWN <- clean_gene_column(dat$ATAC_DOWN)

  # Overlapping genes (unique)
  genes_up_up     <- sort(intersect(bRNA_UP,   ATAC_UP))
  genes_down_down <- sort(intersect(bRNA_DOWN, ATAC_DOWN))
  genes_up_down   <- sort(intersect(bRNA_UP,   ATAC_DOWN))
  genes_down_up   <- sort(intersect(bRNA_DOWN, ATAC_UP))

  gene_lists[[geno]] <- list(
    UP                = genes_up_up,
    DOWN              = genes_down_down,
    bRNA_UP_ATAC_DOWN = genes_up_down,
    bRNA_DOWN_ATAC_UP = genes_down_up
  )

  cat(sprintf("\n=== %s ===\n", geno))
  cat(sprintf("  UP: %d | DOWN: %d | bRNA\u2191/ATAC\u2193: %d | bRNA\u2193/ATAC\u2191: %d\n",
              length(genes_up_up), length(genes_down_down),
              length(genes_up_down), length(genes_down_up)))

  results <- rbind(results, data.frame(
    Genotype  = geno,
    Category  = c("UP", "DOWN", "bRNA\u2191/ATAC\u2193", "bRNA\u2193/ATAC\u2191"),
    Count     = c(length(genes_up_up), length(genes_down_down),
                  length(genes_up_down), length(genes_down_up)),
    stringsAsFactors = FALSE
  ))
}

# Set category and genotype order
results$Category <- factor(results$Category,
                           levels = c("UP", "DOWN", "bRNA\u2191/ATAC\u2193", "bRNA\u2193/ATAC\u2191"))

results$Genotype <- factor(results$Genotype,
                           levels = c("Ogg1", "APP/PS1", "AO"))

# Same y-axis maximum for all genotype panels
y_axis_max   <- 175

# Genotype colors
genotype_colours <- c(
  "Ogg1"    = "#543c6d",
  "APP/PS1" = "#056475",
  "AO"      = "#f57a0c"
)

# Lighter tint of a color (used for the discordant categories)
tint <- function(hex, fac = 0.45) {
  rgb_vals <- col2rgb(hex) / 255
  tinted   <- rgb_vals + (1 - rgb_vals) * fac
  rgb(tinted[1], tinted[2], tinted[3])
}

# Bar colors: genotype color for UP/DOWN, tint for the discordant categories
results <- results %>%
  mutate(
    label     = ifelse(Count > 0, as.character(Count), ""),
    base_col  = genotype_colours[as.character(Genotype)],
    bar_fill  = ifelse(
      Category %in% c("UP", "DOWN"),
      base_col,
      mapply(tint, base_col)
    )
  )

# Create bar plot
p <- ggplot(results, aes(x = Category, y = Count, fill = bar_fill)) +

  geom_col(
    width     = 0.62,
    colour    = "black",
    linewidth = 0.25
  ) +

  # Count labels above the bars, same size as the 12 pt axis text
  geom_text(
    aes(label = label),
    vjust  = -0.45,
    size   = 12 / .pt,
    colour = "black",
    family = "Helvetica"
  ) +

  scale_fill_identity() +

  scale_y_continuous(
    limits = c(0, y_axis_max * 1.12),  # headroom for count labels
    breaks = c(0, 25, 50, 75, 100, 125, 150, 175),
    expand = expansion(mult = c(0.015, 0))  # small gap below 0
  ) +

  facet_wrap(~Genotype, nrow = 1, axes = "all") +  # axes on every panel

  labs(
    x = NULL,
    y = "Overlapping DEGs / DARs"
  ) +

  theme_classic(base_size = 13, base_family = "Helvetica") +

  theme(
    aspect.ratio      = 1,  # square panels

    strip.text        = element_text(size = 13, face = "bold.italic",
                                     colour = "black", margin = margin(t = 3, b = 3)),
    strip.background  = element_rect(fill = "grey92", colour = "black", linewidth = 0.45),
    strip.clip        = "on",

    axis.text.x       = element_text(size = 12, colour = "black",
                                     angle = 38, hjust = 1, vjust = 1),
    axis.text.y       = element_text(size = 12, colour = "black"),
    axis.title.y      = element_text(size = 13, colour = "black",
                                     margin = margin(r = 5)),
    axis.line         = element_blank(),
    axis.ticks        = element_line(colour = "black", linewidth = 0.35),
    axis.ticks.length = unit(2.5, "pt"),

    panel.border      = element_rect(fill = NA, colour = "black", linewidth = 0.5),
    panel.grid        = element_blank(),
    panel.spacing     = unit(10, "pt"),

    plot.margin       = margin(t = 8, r = 8, b = 6, l = 6, unit = "pt"),
    plot.background   = element_rect(fill = "white", colour = NA)
  )

# Show plot
print(p)

# Save plot (PDF, PNG and SVG)
ggsave("ATAC_bRNA_intersection_all_genotypes.pdf", p,
       width = 9.6, height = 4.0, device = cairo_pdf)

ggsave("ATAC_bRNA_intersection_all_genotypes.png", p,
       width = 9.6, height = 4.0, dpi = 600, bg = "white")

ggsave("ATAC_bRNA_intersection_all_genotypes.svg", p,
       width = 9.6, height = 4.0, device = svglite::svglite, bg = "white")

cat("\nSaved: ATAC_bRNA_intersection_all_genotypes.pdf / .png / .svg\n")

# Print summary table
cat("\n=== SUMMARY TABLE ===\n")
print(
  results %>%
    select(Genotype, Category, Count) %>%
    arrange(Genotype, Category),
  row.names = FALSE
)

# Excel export: summary sheet and one sheet per genotype with gene names
cat_labels <- c(
  UP                = "UP (bRNA-up / ATAC-up)",
  DOWN              = "DOWN (bRNA-down / ATAC-down)",
  bRNA_UP_ATAC_DOWN = "Discordant (bRNA-up / ATAC-down)",
  bRNA_DOWN_ATAC_UP = "Discordant (bRNA-down / ATAC-up)"
)
cat_keys <- names(cat_labels)

wb_out <- createWorkbook()

hdr_style  <- createStyle(halign = "CENTER", valign = "center", wrapText = TRUE, textDecoration = "bold")
bold_style <- createStyle(halign = "LEFT",   valign = "center", textDecoration = "bold")
cnt_style  <- createStyle(halign = "CENTER", valign = "center", textDecoration = "bold")
gene_style <- createStyle(halign = "LEFT",   valign = "center")
note_style <- createStyle(halign = "LEFT")

# Summary sheet
addWorksheet(wb_out, "Summary")

col_names_sum <- c("Genotype (vs. WT)",
                   "UP (bRNA-up / ATAC-up)",
                   "DOWN (bRNA-down / ATAC-down)",
                   "Discordant (bRNA-up / ATAC-down)",
                   "Discordant (bRNA-down / ATAC-up)",
                   "Total")
writeData(wb_out, "Summary", t(col_names_sum),
          startRow = 1, startCol = 1, colNames = FALSE)
addStyle(wb_out, "Summary", hdr_style, rows = 1, cols = 1:6, gridExpand = TRUE)
setRowHeights(wb_out, "Summary", rows = 1, heights = 38)

for (i in seq_along(names(gene_lists))) {
  geno   <- names(gene_lists)[i]
  gl     <- gene_lists[[geno]]
  counts <- sapply(cat_keys, function(k) length(gl[[k]]))
  row_df <- data.frame(as.list(c(Genotype = geno, counts, Total = sum(counts))))
  writeData(wb_out, "Summary", row_df, startRow = 1 + i, startCol = 1, colNames = FALSE)
  addStyle(wb_out, "Summary", bold_style, rows = 1 + i, cols = 1)
  addStyle(wb_out, "Summary", cnt_style,  rows = 1 + i, cols = 2:6, gridExpand = TRUE)
  setRowHeights(wb_out, "Summary", rows = 1 + i, heights = 18)
}

setColWidths(wb_out, "Summary", cols = 1,   widths = 18)
setColWidths(wb_out, "Summary", cols = 2:6, widths = 26)

note_row <- length(gene_lists) + 3
writeData(wb_out, "Summary",
          "Note: Intersections computed as unique gene set overlap. All comparisons vs. WT.",
          startRow = note_row, startCol = 1)
addStyle(wb_out, "Summary", note_style, rows = note_row, cols = 1)

# One sheet per genotype
sheet_names <- c("Ogg1" = "Ogg1", "APP/PS1" = "APP-PS1", "AO" = "AO")

for (geno in names(gene_lists)) {
  sname <- sheet_names[geno]
  gl    <- gene_lists[[geno]]

  addWorksheet(wb_out, sname)

  # Column headers
  writeData(wb_out, sname, t(unname(cat_labels)),
            startRow = 1, startCol = 1, colNames = FALSE)
  addStyle(wb_out, sname, hdr_style, rows = 1, cols = 1:4, gridExpand = TRUE)
  setRowHeights(wb_out, sname, rows = 1, heights = 38)

  # Count row
  for (ci in 1:4) {
    writeData(wb_out, sname, paste0("n = ", length(gl[[cat_keys[ci]]])),
              startRow = 2, startCol = ci, colNames = FALSE)
    addStyle(wb_out, sname, cnt_style, rows = 2, cols = ci)
  }
  setRowHeights(wb_out, sname, rows = 2, heights = 18)

  # Gene lists
  max_genes <- max(sapply(cat_keys, function(k) length(gl[[k]])), 0)
  for (ri in seq_len(max_genes)) {
    for (ci in 1:4) {
      gene_vec <- gl[[cat_keys[ci]]]
      val      <- if (ri <= length(gene_vec)) gene_vec[ri] else ""
      writeData(wb_out, sname, val, startRow = 2 + ri, startCol = ci, colNames = FALSE)
      addStyle(wb_out, sname, gene_style, rows = 2 + ri, cols = ci)
    }
    setRowHeights(wb_out, sname, rows = 2 + ri, heights = 16)
  }

  setColWidths(wb_out, sname, cols = 1:4, widths = 30)
}

xlsx_path <- "ATAC_bRNA_overlapping_genes.xlsx"
saveWorkbook(wb_out, xlsx_path, overwrite = TRUE)
cat(sprintf("\nExcel saved: %s\n", xlsx_path))
