# Code depository

Code for the figures and the sequencing data processing of the manuscript:
**Ogg1 regulates context-dependent glial states in Alzheimer's disease**
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
files are written to the same folder. Input and output file names are listed
in the header of each script.

| Figure | Folder | Script |
|---|---|---|
| Fig. 1i | Fig - Dotplots | `R_doplot_bRNA.R` |
| Fig. 1i | Fig - prevDEG-prevDAR correlation plots | `R_prevDEG_correlation_plots.R` |
| Fig. 2g | Fig2g - Intersection bRNA-ATAC | `R_intersection_bRNA_ATAC.R` |
| Fig. 4f | Fig4f - DAM-HOMEO-Microglia | `R_02_DAM_expression_heatmap.R` |
| Fig. 4f | Fig4f - DAM-HOMEO-Microglia | `R_02_HOMEO_expression_heatmap.R` |
| Fig. 4g | Fig - Dotplots | `R_doplot_Microglia.R` |
| Fig. 4g | Fig - prevDEG-prevDAR correlation plots | `R_prevDEG_correlation_plots.R` |
| Fig. 4i | Fig - hADAD overlap | `R_hAD_Micro.R` |
| Fig. 5d | Fig - Dotplots | `R_doplot_Astro.R` |
| Fig. 5d | Fig - prevDEG-prevDAR correlation plots | `R_prevDEG_correlation_plots.R` |
| Fig. 5f | Fig - hADAD overlap | `R_hAD_Astro.R` |
| Fig. 5j | Fig - Dotplots | `R_doplot_Oligo.R` |
| Fig. 5j | Fig - prevDEG-prevDAR correlation plots | `R_prevDEG_correlation_plots.R` |
| Fig. 5l | Fig - hADAD overlap | `R_hAD_Oligo_hADAD_top26-pvalinhAD.R` |
| Fig. 6a | Fig6a - Concordance heatmaps | `R_concordance_heatmaps.R` |
| Fig. 6b | Fig6b - Chord diagrams | `R_Chord_diagrams.R` |
| Suppl. Fig. 5f | Fig - Dotplots | `R_doplot_ATAC.R` |
| Suppl. Fig. 5f | Fig - prevDEG-prevDAR correlation plots | `R_prevDEG_correlation_plots.R` |
| Suppl. Fig. 7b | Suppl.Fig7b - Propeller - all proportions | `R_propeller_analysis.R` |
| Suppl. Fig. 9e | Suppl.Fig9e-10e - DAA-DOL overlap | `R_DAA.R` |
| Suppl. Fig. 10e | Suppl.Fig9e-10e - DAA-DOL overlap | `R_DOL.R` |

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
