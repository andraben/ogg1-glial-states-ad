# Bulk RNA-seq processing

Processing of the bulk RNA-seq data from raw reads to count tables,
differential expression, PCA, DEG heatmap and the prevented DEGs (prevDEGs).

Code written by Teri Sakshaug (teri.sakshaug@ntnu.no), NTNU.

The commands are given in the order they were run. Shell blocks were run on
Ubuntu 20.04, R blocks in one R session (R 4.3.2; the full package list is in
`sessionInfo_bRNA-Seq-DE.txt`). Paths starting with `/path/to/` must be
adapted.

## What is used where

| Section | Output | Used for |
|---|---|---|
| Quality & alignment | sorted, indexed BAM files; MultiQC report | all downstream steps |
| Differential expression | `Counts_Rsubread_GRCm38.102GTF_bRNA.txt`, `EdgeR_featureCounts_DEresults.xlsx`, `EdgeR_featureCounts_DEresults_FilterP01FC1.xlsx` | DEG numbers in the text; inputs of Fig. 2g and Fig. 6a (see `../README.md`) |
| PCA | `RNA_PCAtools_log2FPKM.svg` | Suppl. Fig. 4a |
| Heatmap | `RNA_DEG_Heatmap.svg` | Fig. 1e |
| prevDEGs | `resDEGs.xlsx` (= prevDEGs, Supplementary Table 1), `bRNA_resDEG_CorrelationPlot.svg` | Fig. 1i (final panel made with `../Fig - prevDEG-prevDAR correlation plots/R_prevDEG_correlation_plots.R`) |

## Software

| Tool | Version |
|---|---|
| FastQC | 0.11.9 |
| FastQ Screen | 0.15.2 |
| HISAT2 | 2.2.1 (index: GRCm38 genome_tran from http://daehwankimlab.github.io/hisat2/download/) |
| samtools | 1.10 (htslib 1.10.2) |
| MultiQC | ? |
| Picard | ? |
| R | 4.3.2 with Rsubread 2.16.1, edgeR 3.42.4, limma 3.56.2, PCAtools, countToFPKM, ComplexHeatmap, dendextend, openxlsx 4.2.5.2, dplyr 1.1.3, ggplot2 3.5.0, ggrepel 0.9.3, svglite |
| Annotation | Ensembl GRCm38 release 102 (`Mus_musculus.GRCm38.102.gtf`) |

## Folder map

```
GenomicRefFiles/
  hisat2_index/                            ${HISAT_INDEX}
  Mus_musculus.GRCm38.102.gtf
  Mus_musculus.GRCm38.102.genelengths.txt  (mean transcript length per gene)

RNA_Analysis/
  QC_Files/                                ${QC_DIR}
    FastQScreen/                           ${FASTSCREEN_DIR}
    FastQC/                                ${FASTQC_DIR}
    Hisat2_SummaryOutput/                  ${HISATOUT_DIR}
    FastQC_AlignedBAMs/                    ${FASTQCB_DIR}
  fastq/                                   ${INPUT_DIR}
    BAMs/                                  ${BAM_DIR}
  Analysis/
    DiffEnrich/                            ${DE_DIR}
      picard_metrics/                      ${PICARD_DIR}
      Counts_Rsubread_GRCm38.102GTF_bRNA.txt
      EdgeR_featureCounts_DEresults.xlsx
      EdgeR_featureCounts_DEresults_FilterP01FC1.xlsx
      sessionInfo_bRNA-Seq-DE.txt
      resDEGs.xlsx
      RNA_PCAtools_log2FPKM.svg
      bRNA_resDEG_CorrelationPlot.svg
      RNA_DEG_Heatmap.svg
```

## Quality and alignment

```sh
# Raw fastq files are in the fastq directory
mkdir fastq
mkdir QC_Files
mkdir Analysis

THREADS=$(nproc)
INPUT_DIR="fastq"
QC_DIR="QC_Files"

# FastQC (v0.11.9) and FastQ Screen (v0.15.2)
FASTQC_DIR="${QC_DIR}/FastQC"
FASTSCREEN_DIR="${QC_DIR}/FastQScreen"

mkdir -p "${FASTQC_DIR}"
mkdir -p "${FASTSCREEN_DIR}"

fastqc \
    --threads "${THREADS}" \
    --outdir "${FASTQC_DIR}" \
    "${INPUT_DIR}"/*.fq.gz

find "${INPUT_DIR}" -name "*.fq.gz" \
    -exec fastq_screen \
    --threads "${THREADS}" \
    --outdir "${FASTSCREEN_DIR}" \
    {} \;

# Alignment with HISAT2 (v2.2.1), GRCm38 genome_tran index
HISAT_INDEX="/path/to/Hisat2/grcm38/genome_tran"
BAM_DIR="${INPUT_DIR}/BAMs"
HISATOUT_DIR="${QC_DIR}/Hisat2_SummaryOutput"

mkdir -p "${BAM_DIR}"
mkdir -p "${HISATOUT_DIR}"

for f in "${INPUT_DIR}"/*_1.fq.gz
do
    prefix=$(basename "$f" _1.fq.gz)
    hisat2 \
        -p "${THREADS}" \
        -x "${HISAT_INDEX}" \
        -1 "${INPUT_DIR}/${prefix}_1.fq.gz" \
        -2 "${INPUT_DIR}/${prefix}_2.fq.gz" \
        2> "${HISATOUT_DIR}/${prefix}.txt" \
        -S "${BAM_DIR}/${prefix}.sam"
done
```

### Convert to BAM, sort and index

```sh
# samtools 1.10
for f in $(ls "${BAM_DIR}"/*.sam)
do
    prefix=${f%.*}
    samtools view -h -b "$f" > "${prefix}.bam"
done

for f in $(ls "${BAM_DIR}"/*.bam)
do
    prefix=${f%.*}
    samtools sort "$f" -o "${prefix}_s.bam"
done

for f in $(ls "${BAM_DIR}"/*_s.bam)
do
    samtools index "${f}"
done
```

### Quality of the final BAM files

```sh
FASTQCB_DIR="${QC_DIR}/FastQC_AlignedBAMs"
mkdir -p "${FASTQCB_DIR}"

fastqc \
    --threads "${THREADS}" \
    --outdir "${FASTQCB_DIR}" \
    "${BAM_DIR}"/*_s.bam

multiqc "${QC_DIR}" --data-dir "${QC_DIR}"

ANALYSIS_DIR="Analysis"
DE_DIR="${ANALYSIS_DIR}/DiffEnrich"
mkdir -p "${DE_DIR}"
```

### Fragment size metrics

Mean fragment length per sample, needed for the FPKM calculation below.

```sh
PICARD_DIR="${DE_DIR}/picard_metrics"
mkdir -p "${PICARD_DIR}"

for f in "${BAM_DIR}"/*_s.bam
do
    prefix=$(basename "$f" .bam)
    java -jar picard.jar CollectInsertSizeMetrics \
        -I "$f" \
        -O "${PICARD_DIR}/${prefix}.txt" \
        -H "${PICARD_DIR}/${prefix}.pdf" \
        -M 0.5
done
```

## Differential expression

Counts with Rsubread featureCounts on protein-coding, non-mitochondrial genes
(Ensembl GRCm38.102), differential expression with edgeR quasi-likelihood
F-tests. Five contrasts are tested: AO vs WT, APP/PS1 vs WT, Ogg1-/- vs WT,
AO vs APP/PS1 and AO vs Ogg1-/-. Significant genes (`_FilterP01FC1`) are
defined as P <= 0.01 and |log2FC| >= 1.

```r
# R version 4.3.2
library(Rsubread)
library(openxlsx)
library(edgeR)
library(dplyr)
library(rtracklayer)

BAM_DIR <- "/path/to/RNA_Analysis/fastq/BAMs"
DE_DIR  <- "/path/to/RNA_Analysis/Analysis/DiffEnrich"
QC_DIR  <- "/path/to/RNA_Analysis/QC_Files"
gtf     <- "/path/to/Mus_musculus.GRCm38.102.gtf"

# Protein-coding genes without the mitochondrial genome
gtf_data <- import(gtf, format = "gtf")
gtf102 <- gtf_data[gtf_data$gene_biotype == "protein_coding" & as.character(seqnames(gtf_data)) != "MT"]
export(gtf102, file.path(DE_DIR, "GRCm38.102_AO_bRNA.gtf"), format = "gtf")

files <- list.files(BAM_DIR, pattern='_s\\.bam$', full.names = TRUE)

fct <- featureCounts(files, annot.ext= file.path(DE_DIR, "GRCm38.102_AO_bRNA.gtf"), isGTFAnnotationFile=TRUE, GTF.featureType="exon", GTF.attrType="gene_name", chrAliases=NULL, useMetaFeatures=TRUE, allowMultiOverlap=FALSE,readExtension5=0, readExtension3=0, read2pos=NULL, countMultiMappingReads=TRUE, fraction=FALSE, minMQS=20, splitOnly=FALSE, nonSplitOnly=FALSE, primaryOnly=FALSE, ignoreDup=TRUE, strandSpecific=0, juncCounts=FALSE, genome=NULL, isPairedEnd=TRUE, requireBothEndsMapped=TRUE, autosort=TRUE, nthreads=10)

write.table(fct$counts, file.path(DE_DIR, "Counts_Rsubread_GRCm38.102GTF_bRNA.txt"), sep="\t",quote=FALSE,col.names=TRUE, row.names=TRUE)
write.table(fct$stat, file.path(QC_DIR, "Stats_Rsubread_GRCm38.102GTF_bRNA.txt"), sep="\t",quote=FALSE,col.names=TRUE, row.names=TRUE)

# Samples are in alphabetical order: AO1-3, APP1-3, Ogg1-3, WT1-3
Genotype <- c(rep("AO", each=3), rep("APP/PS1", each=3), rep("Ogg1-/-", each=3), rep("WT", each=3))
counts <- read.delim(file.path(DE_DIR, "Counts_Rsubread_GRCm38.102GTF_bRNA.txt"), sep='\t', header=TRUE, row.names=1)
d <- DGEList(counts=counts,lib.size=colSums(counts), group=factor(Genotype))

design <- model.matrix(~0+group, data=d$samples)
colnames(design) <- levels(d$samples$group)

drop <- which(apply(cpm(d), 1, max) < .09)
dge <- d[-drop,]
dge <- calcNormFactors(dge, method="TMM")
dge <- estimateDisp(dge,design)
fit <- glmQLFit(dge, design)

# Design columns (alphabetical factor levels): AO, APP/PS1, Ogg1-/-, WT
contrasts <- list(
  AO   = c(1, 0, 0,-1),
  APP  = c(0, 1, 0,-1),
  OGG  = c(0, 0, 1,-1),
  AOAP = c(1,-1, 0, 0),
  AOOG = c(1, 0,-1, 0))

results <- lapply(contrasts, function(x) glmQLFTest(fit, contrast = x))
results <- lapply(results, function(x) topTags(x, n = Inf)$table)
results <- lapply(results, as.data.frame)
results <- lapply(results, function(x) x[order(x$logFC), ])

sig_results <- lapply(results, function(x) {x %>% filter(PValue <= 0.01) %>% filter(abs(logFC) >= 1)})

compnames <- c("AOvsWT", "APPvsWT", "OGGvsWT", "AOvsAPP", "AOvsOGG")
names(results) <- compnames
names(sig_results) <- compnames

write.xlsx(results, file.path(DE_DIR, "EdgeR_featureCounts_DEresults.xlsx"), rowNames=TRUE, colNames=TRUE)
write.xlsx(sig_results, file.path(DE_DIR, "EdgeR_featureCounts_DEresults_FilterP01FC1.xlsx"), rowNames=TRUE, colNames=TRUE)

sink(file.path(DE_DIR, "sessionInfo_bRNA-Seq-DE.txt"), append = FALSE)
print(sessioninfo::session_info())
sink()
```

## PCA (Suppl. Fig. 4a)

FPKM values from the counts, gene lengths and the Picard mean fragment
lengths; PCA on log2(FPKM + 1).

```r
library(PCAtools)
library(edgeR)
library(svglite)
library(countToFPKM)

DE_DIR  <- "/path/to/RNA_Analysis/Analysis/DiffEnrich"
lengths <- "/path/to/Mus_musculus.GRCm38.102.genelengths.txt"
countData <- read.csv(file.path(DE_DIR, "Counts_Rsubread_GRCm38.102GTF_bRNA.txt"), header=TRUE, sep="\t", row.names=1)
PICARD_DIR <- file.path(DE_DIR, "picard_metrics")
gene.annotations <- read.csv(lengths, header=TRUE, sep="\t", row.names=1)

# Mean insert size per sample (same alphabetical sample order as the counts)
picard_files <- list.files(PICARD_DIR, pattern = "_s\\.txt$", full.names = TRUE)
mean_lengths <- lapply(picard_files, function(file) {
  x <- readLines(file)
  header_line <- grep("^MEDIAN_INSERT_SIZE", x)
  metrics <- read.delim(file, skip = header_line - 1, nrows = 1)
  data.frame(sample = basename(file), meanFragmentLength = metrics$MEAN_INSERT_SIZE)})

meanFragmentLength <- do.call(rbind, mean_lengths)
m <- match(rownames(countData), gene.annotations$GeneName)
length <- gene.annotations$mean[m]
sum(is.na(length))
complete_cases <- complete.cases(length)
countData <- countData[complete_cases, ]
length <- length[complete_cases]

meanFragmentLength <- as.numeric(meanFragmentLength[,"meanFragmentLength"])
fpkms <- fpkm(countData, length, meanFragmentLength)

Genotype <- c(rep("AO", each=3), rep("APP/PS1", each=3), rep("Ogg1-/-", each=3), rep("WT", each=3))
Replicate <- rep(1:3, times = 4)
meta = data.frame(Genotype, Replicate)
rownames(meta) <- colnames(countData)

logfpkm <- log2(fpkms+1)
p <- pca(logfpkm, metadata = meta, removeVar = 0.1)

svglite(file.path(DE_DIR, "RNA_PCAtools_log2FPKM.svg"), width = 8, height = 8)
biplot(p, x = "PC1", y = "PC2", pointSize=5, colby='Genotype', colkey = c('WT' = '#5e5e5e', 'APP/PS1' = '#056475', 'Ogg1-/-' = '#543c6d', 'AO' = '#f57a0c'), lab=p$metadata[,2], legendPosition='right', title='PCA of RNA-sequencing', caption='log2 FPKM')
dev.off()
```

## Heatmap (Fig. 1e)

Z-scored log2 FPKM of all genes with P <= 0.05 and |log2FC| >= 0.5 in at
least one of the three comparisons against WT, clustered into four groups.
Continues the PCA session (`logfpkm`). The count columns must be named
`WT1 ... AO3` for `sample_order` to match.

```r
library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(grid)
library(openxlsx)
library(dplyr)
library(svglite)

DE_DIR <- "/path/to/RNA_Analysis/Analysis/DiffEnrich"
DE_file <- file.path(DE_DIR, "EdgeR_featureCounts_DEresults.xlsx")

sheet_names <- getSheetNames(DE_file)
results <- lapply(sheet_names, function(sheet) {read.xlsx(DE_file, sheet = sheet, rowNames=TRUE)})

sig_results <- lapply(results, function(x) {x %>% filter(PValue <= 0.05) %>% filter(abs(logFC) >= .5)})
names(sig_results) <- sheet_names

genes <- unique(unlist(lapply(sig_results[1:3], rownames)))
DEGfpkm <- logfpkm[rownames(logfpkm) %in% genes, ]
DEGfpkm_zscore <- t(scale(t(DEGfpkm)))
sample_order <- c("WT1", "WT2", "WT3", "APP1", "APP2", "APP3", "Ogg1", "Ogg2", "Ogg3","AO1", "AO2", "AO3")
DEGfpkm_zscore <- DEGfpkm_zscore[, sample_order]
dim(DEGfpkm_zscore)

hc <- hclust(dist(DEGfpkm_zscore))
clusters <- cutree(hc, k = 4)

cluster_colors <- c("1" = "#E69F00", "2" = "#E69F00", "3" = "#56B4E9", "4" = "#009E73")
dend <- as.dendrogram(hc)
dend_clusters <- clusters[labels(dend)]

dend <- color_branches(dend, clusters = dend_clusters, col = cluster_colors[as.character(sort(unique(dend_clusters)))])
dend <- set(dend, "branches_lwd", 0.6)

annotation_row <- rowAnnotation(Cluster = factor(clusters[rownames(DEGfpkm_zscore)], levels = 1:4), col = list(Cluster = cluster_colors))

col_fun <- colorRamp2(c(-2, 0, 2), c("#0068ba", "#000000", "#fd7d00"))

RNAheat <- Heatmap(DEGfpkm_zscore,
		name = "Z-score",
		col = col_fun,
		cluster_rows = dend,
  		row_dend_reorder = FALSE,
  		row_dend_width = unit(25, "mm"),
		cluster_columns = FALSE,
  		show_column_names = TRUE,
  		show_row_names = FALSE,
		border = FALSE,
  		column_title = "Gene expression (Z-score)" )

svglite(file.path(DE_DIR, "RNA_DEG_Heatmap.svg"), width = 8, height = 8)
RNAheat
dev.off()
```

## Prevented DEGs (prevDEGs; Fig. 1i, Supplementary Table 1)

Genes significant in both APP/PS1 vs WT and AO vs APP/PS1 (P <= 0.01,
|log2FC| >= 1). The compensation is -log2FC(AO vs APP/PS1) / log2FC(APP/PS1
vs WT); genes with a positive value (change in APP/PS1 counteracted in AO)
are the prevDEGs. `resDEGs.xlsx` is written before this filter and contains
all overlapping genes. The correlation plot drawn here is the draft of
Fig. 1i; the published panel is made from `prevDEGs_bRNA.txt` with
`R_prevDEG_correlation_plots.R` in this repository.

```r
library(dplyr)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(svglite)

DE_DIR <- "/path/to/RNA_Analysis/Analysis/DiffEnrich"
DE_file <- file.path(DE_DIR, "EdgeR_featureCounts_DEresults_FilterP01FC1.xlsx")

sheet_names <- getSheetNames(DE_file)
results <- lapply(sheet_names, function(sheet) {read.xlsx(DE_file, sheet = sheet, rowNames=TRUE)})

cols <- c("logFC", "PValue", "FDR")
appWT <- results[[2]][, cols]
aoAPP <- results[[4]][, cols]

resDEGs <- merge(appWT, aoAPP, by = "row.names", suffix = c("_APP-vs-WT", "_AO-vs-APP"))
resDEGs$PercRescue <- -resDEGs[,5]/resDEGs[,2]
write.xlsx(resDEGs, file = file.path(DE_DIR, "resDEGs.xlsx"), rowNames = FALSE)
resDEGs <- resDEGs %>% filter(PercRescue>0)

resDEGs$Score <- (resDEGs$PercRescue*100)
resDEGs <- resDEGs %>% mutate(LogScore = log10(Score + 1e-3))

genes_to_label <- c( "Nlrp1a","Ifit1","Irf7","Nlrc5","Ifi204","Zbp1","Cd48","Cst7")
deg_labels <- resDEGs[resDEGs$Row.names %in% genes_to_label, ]

min_val <- floor(min(resDEGs$Score) / 10) * 10
max_val <- ceiling(max(resDEGs$Score) / 10) * 10
mid_val <- (min_val + max_val) / 2
lims <- ceiling(max(abs(c(resDEGs[,5], resDEGs[,2]))))

rna <- ggplot(resDEGs, aes(x = resDEGs[,5], y = resDEGs[,2])) +
  geom_abline(slope = -1, intercept = 0, color = "grey50", linewidth = 0.42, linetype = "solid") +
  geom_point(aes(color = Score), shape = 16, size = 4.5, alpha = 0.8) +
  geom_point(data = deg_labels, aes(x = `logFC_AO-vs-APP`, y = `logFC_APP-vs-WT`), shape = 21, color = "black", size = 5, stroke = 1.5) +
  geom_text_repel(data = deg_labels,aes(x = `logFC_AO-vs-APP`, y = `logFC_APP-vs-WT`, label = `Row.names`), size = 6, max.overlaps = Inf, box.padding = 0.6, force = 3) +
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
  labs(x = "log2FC(AOvsAPP)", y = "log2FC(APPvsWT)") +
  xlim(-lims, lims) +
  ylim(-lims, lims)

svglite(file.path(DE_DIR, "bRNA_resDEG_CorrelationPlot.svg"), width = 8, height = 8)
rna
dev.off()
```
