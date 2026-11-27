# ------------------------------------------------------------------------------
# Suppl. Fig. 9e - Heatmap of DAA gene expression in astrocytes
#
# Input:  02_Astrocytes.csv
# Output: DAA.svg
# ------------------------------------------------------------------------------

# Load libraries
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# Read input data
df <- read.csv("02_Astrocytes.csv") %>%
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
bold_genes <- c("ENSMUSG00000121395", "Hormad2", "Gm43154", "Rarres2", "Car12", "Galnt14", "Iqcj", "Gm21972", "Trabd2b", "Serpinf1", "1700094M23Rik", "Trim56", "Dph3", "Ccbe1", "Nrg3os", "Ypel4", "Nr4a3", "Gm11681", "Col25a1", "Zfp804b", "Gm9962", "Cdh13", "Filip1", "Gm10605", "Cyp4v3", "Thbs4", "ENSMUSG00000120798", "ENSMUSG00000121275", "Galnt4", "Serinc2", "A830018L16Rik", "Kcnk13", "A730036I17Rik", "Atp2b3", "ENSMUSG00000120042", "Ppm1m", "Ntng1", "Sult4a1", "Stk17b", "Rap1gap2", "Gm42744", "Mmp24", "Slc26a4", "Unc5d", "Actn1", "Snn", "Myo5b", "Frmd5", "Gm12339", "Cacng3", "Gprasp2", "Ppp1r3d", "Agtrap", "Nectin1", "4930517G24Rik", "Sh3glb2", "Chkb", "Gm28050", "Gng2", "Plppr1", "Ccdc120", "Imp4", "Rasgrf2", "Gm42549", "Gja1", "Cntnap5a", "Rnaset2a", "Id4", "Gm31508", "Cfap20dc", "Tmem117", "Lcat", "Cntn5", "Gm49776", "Gm20754", "ENSMUSG00000119995", "Bhlhe22", "Ttll11", "Agpat4", "Adarb2", "Nkain2", "Chrm1", "Phf24", "Gm49445", "Cmss1", "Ksr2", "Slc2a6", "Fmnl3", "Gpr37l1", "Zfp831", "1700111E14Rik", "Cd63", "Il16", "Scn9a", "Ntrk3", "Kank2", "Gm9916", "Gm38353", "Nsg1", "Gpc1", "Ddr1", "ENSMUSG00000120976", "Prmt8", "Ppp1r26", "Vipr1", "A930029G22Rik", "Neurl1b", "Eif2ak2", "Gm35552", "Gm30092", "Arhgap10", "Tmem47", "Mtln", "Gm30655", "Rerg", "Dera", "Pdia5", "Tep1", "Kcnh1", "Ppp2r2c", "Clstn2", "Ankrd33b", "Neurl1a", "Col16a1", "Celf3", "Pak6", "Sncb", "Hs3st1", "Trim12c", "Slc35f3", "Prkcz", "Susd1", "Ptchd1", "Cfap97d2", "Rph3a", "Trank1", "Pip5k1b", "4930532I03Rik", "Gm13031", "St8sia4", "Chst9", "Cdkl1", "Spock3", "Spr", "B930018H19Rik", "Mfsd13b", "Plxna4", "Tmem266", "Mctp1", "Ptpre", "Htr4", "Adamts17", "Gm15624", "Cacna1c", "ENSMUSG00000121191", "Gm43666", "Nmnat2", "Cacnb3", "Dtx3l", "Susd4", "Slit3", "Fndc1", "Cdk5r1", "Fgf9", "Slc2a10", "Kcnj6", "Rbfox3", "Inka2", "4732463B04Rik", "Gm32786", "Arl6ip1", "Pcgf6", "Gdpd2", "Arhgap20", "Lgals3bp", "Slc7a14", "Lrrtm4", "Gm44724", "Ptk2b", "Bsn", "Tspan4", "Ajap1", "Thpo", "Sez6l", "Pxmp4", "Chil1", "Vstm2b", "Cabyr", "Syt7", "Cnih3", "Akain1", "Synpr", "Gm10561", "Dbpht2", "Clu", "S100a6", "Stk26", "Mapk12", "A730091E23Rik", "Lingo2", "Mgat5b", "Cdh9", "Irak4", "Rhoj", "Asap3", "Slc16a7", "Kctd12", "Trpc6", "Cmtm5", "Gm20743", "Sh3pxd2a", "Ppfia2", "Kcnj4", "Znhit3", "Lratd1", "Parp9", "Rimbp2", "Gsg1l", "Cemip", "Pde1b", "Paqr9", "Igsf10", "Cst3", "Kcng2", "Cyp11a1", "Mpp3", "Gm6999", "Ablim3", "Oasl2", "Fkbp5", "Slc29a4", "Aqp4", "Mfsd4a", "Tafa2", "Smad3", "Rnf213", "Dock10", "Rundc3a", "Stk32b", "Lrp10", "Trim34a", "Ptpro", "Ntf3", "Pgbd5", "Shisa9", "Napepld", "Ehd2", "Egln2", "Gm14455", "Rbp1", "ENSMUSG00000120362", "Rlbp1", "Tspan18", "Gnal", "Cracdl", "Smoc2", "Gpnmb", "Pcsk1n", "Mx2", "S100b", "Stard10", "Tiam1", "Prox1os", "Plekha4", "Ptms", "Paqr4", "Stk32c", "Gm21954", "Ifitm7", "Rad51b", "Gm30524", "Rigi", "Cd9", "Sv2c", "Gadl1", "Gm4804", "Osmr", "Gbp7", "Cacna1h", "Htra4", "Serpina3n", "Gfap", "Cd274", "Sema5a", "Tnfaip6", "Ifih1", "Mkx", "Tunar", "Crlf1", "Slc14a1", "Hs3st5", "Mr1", "Ddx60", "Gbp3", "Nlrc5", "C4b", "Tlr4", "Gm42047", "Ggta1", "Ip6k3", "Gm4951", "Bcl3", "Colec12")

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
svg("DAA.svg", width = 5, height = 5, bg = "white")
draw(ht, heatmap_legend_side = "right")
dev.off()
