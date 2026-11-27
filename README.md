# Code depository

Code for the figures and the sequencing data processing of the manuscript:
Ogg1 regulates context-dependent glial states in Alzheimer's disease
(bioRxiv preprint: https://www.biorxiv.org/content/10.64898/2026.11.06.764857v1).

## Contributors

- Andreas Abentung (andreas.abentung@ahus.no), Department of Pathology,
  Akershus University Hospital, Lørenskog, Norway: figure scripts in the
  folders below
- Teri Sakshaug (teri.sakshaug@ntnu.no), Department of Neuromedicine and Movement Science,
  Norwegian University of Science and Technology (NTNU), Trondheim, Norway: bulk RNA-seq and ATAC-seq data
  processing (`Data processing/`) and the prevDEG / prevDAR correlation plots
  (Fig. 1i, 4g, 5d, 5j and Suppl. Fig. 5f)

## Figure scripts

Each folder contains the R script(s) for one or more figures together with the
input files. Run a script with its own folder as working directory. Output
files are written to the same folder.

| Figure | Folder | Script | Input | Output |
|---|---|---|---|---|
| Fig. 1i | Fig - Dotplots | `R_doplot_bRNA.R` | `Average-FPKM_bRNA.csv` | `R_doplot_bRNA.svg` |
| Fig. 1i | Fig - prevDEG-prevDAR correlation plots | `R_prevDEG_correlation_plots.R` | `prevDEGs_bRNA.txt` | `CorPlots_prevDEGs_bRNA_BIG.svg` |
| Fig. 2g | Fig2g - Intersection bRNA-ATAC | `R_intersection_bRNA_ATAC.R` | `Ogg1_UP-DOWN_bRNA-epigenetics.csv`, `APP_UP-DOWN_bRNA-epigenetics.csv`, `AO_UP-DOWN_bRNA-epigenetics.csv` | `ATAC_bRNA_intersection_all_genotypes.pdf / .png / .svg`, `ATAC_bRNA_overlapping_genes.xlsx` |
| Fig. 4f | Fig4f - DAM-HOMEO-Microglia | `R_02_DAM_expression_heatmap.R` | `DAM_expression_results.csv` | `DAM_expression_heatmap.svg` |
| Fig. 4f | Fig4f - DAM-HOMEO-Microglia | `R_02_HOMEO_expression_heatmap.R` | `HOMEO_expression_results.csv` | `HOMEO_expression_heatmap.svg` |
| Fig. 4g | Fig - Dotplots | `R_doplot_Microglia.R` | `Microglia.csv` | `R_doplot_Microglia.svg` |
| Fig. 4g | Fig - prevDEG-prevDAR correlation plots | `R_prevDEG_correlation_plots.R` | `prevDEGs_scMicroglia.txt` | `CorPlots_prevDEGs_Microglia_BIG.svg` |
| Fig. 4i | Fig - hADAD overlap | `R_hAD_Micro.R` | `Microglia_R.csv` | `hADAD_Microglia.svg` |
| Fig. 5d | Fig - Dotplots | `R_doplot_Astro.R` | `Astrocytes.csv` | `R_doplot_Astro.svg` |
| Fig. 5d | Fig - prevDEG-prevDAR correlation plots | `R_prevDEG_correlation_plots.R` | `prevDEGs_scAstrocytes.txt` | `CorPlots_prevDEGs_Astrocytes_BIG.svg` |
| Fig. 5f | Fig - hADAD overlap | `R_hAD_Astro.R` | `Astrocytes_R.csv` | `hADAD_Astro.svg` |
| Fig. 5j | Fig - Dotplots | `R_doplot_Oligo.R` | `Oligodendrocytes.csv` | `R_doplot_Oligo.svg` |
| Fig. 5j | Fig - prevDEG-prevDAR correlation plots | `R_prevDEG_correlation_plots.R` | `prevDEGs_scOligodendrocytes.txt` | `CorPlots_prevDEGs_Oligodendrocytes_BIG.svg` |
| Fig. 5l | Fig - hADAD overlap | `R_hAD_Oligo_hADAD_top26-pvalinhAD.R` | `Oligo_R_hADAD-top26.csv` | `hADAD_Oligo.svg` |
| Fig. 6a | Fig6a - Concordance heatmaps | `R_concordance_heatmaps.R` | `prevDEG_<celltype>_LogFCs.csv`, `<genotype>_<celltype>_allgenes_LogFCs.csv` | `Figure_1_spearman_rho_<approach>_<celltype>.png / .svg`, `Figure_3_median_abs_log2FC_<approach>_<celltype>.png / .svg` |
| Fig. 6b | Fig6b - Chord diagrams | `R_Chord_diagrams.R` | `GO_terms.xlsx` | `Chord_diagrams_all.svg` |
| Suppl. Fig. 5f | Fig - Dotplots | `R_doplot_ATAC.R` | `ATAC_AvgExp.csv` | `R_doplot_ATAC.svg`, `ATAC_Averaged_Peaks_Per_Gene.csv`, `ATAC_ZScored_Average_Expression.csv` |
| Suppl. Fig. 5f | Fig - prevDEG-prevDAR correlation plots | `R_prevDEG_correlation_plots.R` | `prevDARs_ATAC.txt` | `CorPlots_prevDARs_ATAC_BIG.svg` |
| Suppl. Fig. 7b | Suppl.Fig7b - Propeller - all proportions | `R_propeller_analysis.R` | `cell_type_proportions.csv` | `propeller_anova.csv`, `propeller_pairwise.csv`, `propeller_statistics.txt`, `cell_type_proportions.pdf` |
| Suppl. Fig. 9e | Suppl.Fig9e-10e - DAA-DOL overlap | `R_DAA.R` | `02_Astrocytes.csv` | `DAA.svg` |
| Suppl. Fig. 10e | Suppl.Fig9e-10e - DAA-DOL overlap | `R_DOL.R` | `03_Oligo.csv` | `DOL.svg` |

`R_prevDEG_correlation_plots.R` makes all five correlation plots in one run.

## Data processing

The folder `Data processing` documents how the sequencing data were processed
from raw reads to the tables the figure scripts start from, as commented
shell and R code in the order it was run:

| File | Content |
|---|---|
| `bRNA-seq_processing.md` | Bulk RNA-seq: QC, HISAT2 alignment, featureCounts, edgeR differential expression, PCA, DEG heatmap, prevDEGs |
| `ATAC-seq_processing.md` | ATAC-seq: trimming, Bowtie 2 alignment, MACS3 peak calling, peak annotation, coverage tracks and TSS profile, featureCounts, edgeR differential accessibility, PCA, prevDARs, Firre genome track, HOMER motifs |
| `sessionInfo_bRNA-Seq-DE.txt`, `sessionInfo_ATAC-Seq-DE.txt` | R and package versions of the two differential analyses |

Each document starts with a table saying which output is used for which
figure. The snRNA-seq data were analyzed in Trailmaker, as described in the
Materials and Methods of the paper; the per-gene expression tables used by the
figure scripts are exports from there.

## Notes

- Fig. 4f: `DAM_expression_results.csv` and `HOMEO_expression_results.csv`
  contain the average expression per genotype and cluster, calculated from the
  microglia Seurat object in a previous step that is not part of this
  repository.
- Fig. 2g and Fig. 6a use the fonts Helvetica and Arial.
- Suppl. Fig. 7b: speckle and limma are Bioconductor packages
  (`BiocManager::install(c("speckle", "limma"))`).

## R packages

circlize, ComplexHeatmap (Bioconductor), dplyr, ggplot2, ggrepel, limma
(Bioconductor), openxlsx, readxl, reshape2, scales, speckle (Bioconductor),
svglite, tidyr, tidyverse

Fig. 2g requires ggplot2 3.5.0 or newer. Packages for the data processing are
listed in the two documents in `Data processing`.
