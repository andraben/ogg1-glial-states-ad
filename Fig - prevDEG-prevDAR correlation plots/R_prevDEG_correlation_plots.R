# ------------------------------------------------------------------------------
# Fig. 1i, 4g, 5d, 5j and Suppl. Fig. 5f - Correlation plots of prevented DEGs
# (prevDEGs) and prevented DARs (prevDARs)
#
# Each plot shows log2FC(AO vs APP/PS1) against log2FC(APP/PS1 vs WT) for the
# genes (or ATAC peaks) whose APP/PS1-induced change is counteracted in AO mice.
# Points are colored by the compensation score
#   Score = -log2FC(AO vs APP/PS1) / log2FC(APP/PS1 vs WT) * 100
# and only genes with a positive score are shown. The grey line is y = -x,
# i.e. 100 % compensation. Selected genes are circled and labeled.
#
# Plots written by Teri Sakshaug (teri.sakshaug@ntnu.no)
#
# Input:  prevDEGs_bRNA.txt               (Fig. 1i, bulk RNA-seq)
#         prevDEGs_scMicroglia.txt        (Fig. 4g, snRNA-seq microglia)
#         prevDEGs_scAstrocytes.txt       (Fig. 5d, snRNA-seq astrocytes)
#         prevDEGs_scOligodendrocytes.txt (Fig. 5j, snRNA-seq oligodendrocytes)
#         prevDARs_ATAC.txt               (Suppl. Fig. 5f, ATAC-seq, per peak)
# Output: CorPlots_prevDEGs_bRNA_BIG.svg
#         CorPlots_prevDEGs_Microglia_BIG.svg
#         CorPlots_prevDEGs_Astrocytes_BIG.svg
#         CorPlots_prevDEGs_Oligodendrocytes_BIG.svg
#         CorPlots_prevDARs_ATAC_BIG.svg
# ------------------------------------------------------------------------------

# Load libraries
library(dplyr)
library(ggplot2)
library(ggrepel)
library(svglite)

# Input files are tab-separated with a header. The RNA-seq files have the gene
# name as row name and two columns, APPvsWT and AOvsAPP. The ATAC-seq file has
# one row per peak: the gene name in the first column (Gene), followed by
# APPvsWT, AOvsAPP and the peak ID. A gene with several prevented peaks is
# therefore plotted and labeled once per peak.

# Read one result file and calculate the compensation score
read_results <- function(file, gene_column = FALSE) {
  if (gene_column) {
    df <- read.table(file, header = TRUE, sep = "\t", row.names = NULL)
    df$Score <- (-df[, 3] / df[, 2]) * 100
  } else {
    df <- read.table(file, header = TRUE, sep = "\t", row.names = 1)
    df$Gene <- rownames(df)
    df$Score <- (-df[, 2] / df[, 1]) * 100
  }
  df %>% filter(Score > 0)
}

# Correlation plot with selected genes labeled
correlation_plot <- function(df, genes_to_label) {
  df_labels <- df[df$Gene %in% genes_to_label, ]

  # Color scale limits, rounded to the nearest 10
  min_val <- floor(min(df$Score) / 10) * 10
  max_val <- ceiling(max(df$Score) / 10) * 10

  # Symmetric axis limits
  lims <- ceiling(max(abs(c(df$AOvsAPP, df$APPvsWT))))

  ggplot(df, aes(x = AOvsAPP, y = APPvsWT)) +
    geom_abline(slope = -1, intercept = 0, color = "grey50", linewidth = 0.42, linetype = "solid") +
    geom_point(aes(color = Score), shape = 16, size = 4.5, alpha = 0.8) +
    geom_point(data = df_labels, aes(x = AOvsAPP, y = APPvsWT), shape = 21, color = "black", size = 5, stroke = 1.5) +
    geom_text_repel(data = df_labels, aes(label = Gene), size = 6, max.overlaps = Inf, box.padding = 0.6, force = 3) +
    scale_color_gradient(
      low = "#F2B705", high = "#7A1FA2",
      limits = c(min_val, max_val), breaks = seq(min_val, max_val, 50),
      name = "Compensation\n(%)\n",
      guide = guide_colorbar(
        frame.colour = "black", frame.linewidth = .35,
        ticks.colour = "black", ticks.linewidth = .35,
        barwidth = unit(.5, "cm"), barheight = unit(3, "cm"),
        title.theme = element_text(size = 16),
        label.theme = element_text(size = 14))) +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text = element_text(size = 16, color = "black"),
      axis.title = element_text(size = 18),
      axis.line = element_blank(),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(0.2, "cm")) +
    labs(x = "log2FC (AOvsAPP/PS1)", y = "log2FC (APP/PS1vsWT)") +
    xlim(-lims, lims) +
    ylim(-lims, lims)
}

# Fig. 1i - bulk RNA-seq
rna <- read_results("prevDEGs_bRNA.txt")
p <- correlation_plot(rna, c("Nlrp1a", "Ifit1", "Irf7", "Nlrc5", "Ifi204", "Zbp1", "Cd48", "Cst7"))
ggsave("CorPlots_prevDEGs_bRNA_BIG.svg", p, width = 7.5, height = 6, device = "svg")

# Fig. 4g - snRNA-seq microglia
micro <- read_results("prevDEGs_scMicroglia.txt")
p <- correlation_plot(micro, c("Lacc1", "Phf11c", "Card6", "Tmem71", "Elf4", "Axl", "Clec7a", "Cd83", "Sash3", "Cst7"))
ggsave("CorPlots_prevDEGs_Microglia_BIG.svg", p, width = 7.5, height = 6, device = "svg")

# Fig. 5d - snRNA-seq astrocytes
astro <- read_results("prevDEGs_scAstrocytes.txt")
p <- correlation_plot(astro, c("Gfap", "Clu", "Aqp4", "Cd63", "Col25a1", "Ntng1", "Plppr1", "Rasgrf2", "Nlrc5", "Bsn"))
ggsave("CorPlots_prevDEGs_Astrocytes_BIG.svg", p, width = 7.5, height = 6, device = "svg")

# Fig. 5j - snRNA-seq oligodendrocytes
oligo <- read_results("prevDEGs_scOligodendrocytes.txt")
p <- correlation_plot(oligo, c("Cd9", "C4b", "Jup", "Snap25", "Mpzl1", "Shank1", "Mobp", "Hepacam", "Sema3d", "Ppp2r2b"))
ggsave("CorPlots_prevDEGs_Oligodendrocytes_BIG.svg", p, width = 7.5, height = 6, device = "svg")

# Suppl. Fig. 5f - ATAC-seq (prevDARs)
atac <- read_results("prevDARs_ATAC.txt", gene_column = TRUE)
p <- correlation_plot(atac, c("Firre", "Nwd2", "Ntrk2", "Fmr1", "Cbln2", "Ppp2r2b", "Alk", "Hepacam"))
ggsave("CorPlots_prevDARs_ATAC_BIG.svg", p, width = 7.5, height = 6, device = "svg")
