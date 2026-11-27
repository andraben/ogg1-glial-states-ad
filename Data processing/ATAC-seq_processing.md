# ATAC-seq processing

Processing of the ATAC-seq data from raw reads to peaks, count tables,
differentially accessible regions (DARs), PCA, coverage tracks, the prevented
DARs (prevDARs) and HOMER motif analysis.

Code written by Teri Sakshaug (teri.sakshaug@ntnu.no), NTNU.

The commands are given in the order they were run. Shell blocks were run on
Ubuntu 20.04, R blocks in R 4.3.2 (the full package list is in
`sessionInfo_ATAC-Seq-DE.txt`). Paths starting with `/path/to/` must be
adapted.

## What is used where

| Section | Output | Used for |
|---|---|---|
| Alignment and quality control | lane-merged BAM per sample, replicate-merged BAM per genotype | all downstream steps |
| Peak calling | `Consensus_Peaks.bed`, `*_Annotated_Genes.tsv`, `Genes_of_Consensus_Peaks.bed` | counting, profile plot, gene lists |
| BigWigs and profile plot | `SE_Profile-over_ConsensusPeakGenes_blacklist_refpoint.svg` | Suppl. Fig. 5b |
| Counts and DARs | `ATAC_Counts.txt`, `EdgeRqlf_featureCounts_DEresults.xlsx`, `EdgeRqlf_featureCounts_DEresults_FilterP05FC1.xlsx` | DAR numbers in the text; inputs of Fig. 2g and Fig. 6a (see `../README.md`) |
| PCA (inside the DAR block) | `ATAC_PCAtools_EdgeRcpm.svg` | Suppl. Fig. 5a |
| prevDARs | `resDARs.xlsx` (= prevDARs, Supplementary Table 1), `resPEAKs.bed`, `ATAC_resDAR_CorrelationPlot.svg` | Suppl. Fig. 5f (final panel made with `../Fig - prevDEG-prevDAR correlation plots/R_prevDEG_correlation_plots.R`) |
| Genome tracks | `ATAC_Firre_GenomeTrack_1XCov_UCSC.svg` | Fig. 2f |
| HOMER motifs | `HOMER_Unique/Motifs/*` | Fig. 2d, Suppl. Fig. 5d |

## Software

| Tool | Version |
|---|---|
| cutadapt | 2.8 |
| FastQC | 0.11.9 |
| Bowtie 2 | 2.3.5.1 (GRCm38 index) |
| samtools | 1.10 (htslib 1.10.2) |
| MACS3 | 3.0.0a7 |
| bedtools | ? |
| deepTools | 3.5.2 |
| featureCounts (Subread) | 2.0.0 |
| HOMER | ? |
| R | 4.3.2 with edgeR 3.42.4, limma 3.56.2, ChIPseeker 1.36.0, TxDb.Mmusculus.UCSC.mm10.knownGene 3.10.0, org.Mm.eg.db 3.17.0, GenomicRanges 1.52.1, rtracklayer 1.60.1, PCAtools, Gviz, openxlsx 4.2.5.2, dplyr 1.1.3, ggplot2 3.5.0, ggrepel 0.9.3, svglite |
| Blacklist | mm10-blacklist.v2 (Boyle lab, https://github.com/Boyle-Lab/Blacklist) |

## Folder map

```
GenomicRefFiles/
  adapters/                                ${ADAPTER_DIR}
  bowtie_index/                            ${BOWTIE2_INDEX}
  mm10-blacklist.v2.bed.gz                 ${BLACKLIST}

ATAC_Analysis/
  QC_Files/                                ${QC_DIR}
    CutAdapt_Summary/                      ${CASUMMARY_DIR}
    FastQC/                                ${FASTQC_DIR}
    Bowtie2_SummaryOutput/                 ${SUMMARY_DIR}
  fastq/                                   ${INPUT_DIR}
    trimmed_fastq/                         ${TRIMMED_DIR}
      BAMs/                                ${BAM_DIR}
        LanesMerged_BAMs/                  ${LANEMERGED_DIR}
        ReplicatesMerged_BAMs/             ${REPLICATEMERGED_DIR}
          bamCoverage_RPKM/                ${BW_DIR}
            *.bw
            computeMatrix_ConsensusPeakGenes_blacklist_refpoint.mat.gz
            SE_Profile-over_ConsensusPeakGenes_blacklist_refpoint.svg
          bamCoverage_1xCov/               ${BW1x_DIR}
            *.bw
            ATAC_Firre_GenomeTrack_1XCov_UCSC.svg
  Analysis/
    PeakCalling_macs3/                     ${PEAK_DIR}
      Consensus_Peaks.bed
      Consensus_Peaks.saf
      *_Annotated_Genes.tsv
      Genes_of_Consensus_Peaks.bed
      resPEAKs.bed
    DiffEnrich/                            ${DE_DIR}
      ATAC_Counts.txt
      EdgeRqlf_featureCounts_DEresults.xlsx
      EdgeRqlf_featureCounts_DEresults_FilterP05FC1.xlsx
      sessionInfo_ATAC-Seq-DE.txt
      resDARs.xlsx
      ATAC_PCAtools_EdgeRcpm.svg
      ATAC_resDAR_CorrelationPlot.svg
      HOMER_All/
      HOMER_Unique/                        ${HOMER_DIR}
        ATAC_Background.bed
        AOvsWT.bed, APPvsWT.bed, OGGvsWT.bed
        Motifs/
```

## Alignment and quality control

```sh
# Raw fastq files are in the fastq directory
mkdir fastq
mkdir QC_Files
mkdir Analysis

# Adapter trimming with cutadapt (v2.8)
ADAPTER_DIR="/path/to/adapters"
INPUT_DIR="fastq"
QC_DIR="QC_Files"
TRIMMED_DIR="${INPUT_DIR}/trimmed_fastq"
CASUMMARY_DIR="${QC_DIR}/CutAdapt_Summary"
THREADS=$(nproc)

mkdir -p "${TRIMMED_DIR}"
mkdir -p "${CASUMMARY_DIR}"

for f in $(ls "${INPUT_DIR}"/*_1.fq.gz | sed 's/_1.fq.gz//' | xargs -n1 basename | sort -u)
do
    cutadapt \
        -j ${THREADS} -m 95 -l 95 \
        -a "file:${ADAPTER_DIR}/NexteraPE-PE.fa" \
        -A "file:${ADAPTER_DIR}/NexteraPE-PE.fa" \
        -o "${TRIMMED_DIR}/${f}_1_trim.fq.gz" \
        -p "${TRIMMED_DIR}/${f}_2_trim.fq.gz" \
        "${INPUT_DIR}/${f}_1.fq.gz" \
        "${INPUT_DIR}/${f}_2.fq.gz" \
        --info-file "${CASUMMARY_DIR}/CutAdapt-Nextera_${f}.txt"
done

# FastQC (v0.11.9)
FASTQC_DIR="${QC_DIR}/FastQC"
mkdir -p "${FASTQC_DIR}"

fastqc \
    --threads "${THREADS}" \
    --outdir "${FASTQC_DIR}" \
    "${TRIMMED_DIR}"/*.fq.gz

# Alignment with Bowtie 2 (v2.3.5.1)
BOWTIE2_INDEX="/path/to/GRCm38/GRCm38"
BAM_DIR="${TRIMMED_DIR}/BAMs"
SUMMARY_DIR="${QC_DIR}/Bowtie2_SummaryOutput"

mkdir -p "${BAM_DIR}"
mkdir -p "${SUMMARY_DIR}"

for f in "${TRIMMED_DIR}"/*_1_trim.fq.gz
do
    prefix=$(basename "$f" _1_trim.fq.gz)
    bowtie2 \
        -p "${THREADS}" \
        -x "${BOWTIE2_INDEX}" \
        -1 "${TRIMMED_DIR}/${prefix}_1_trim.fq.gz" \
        -2 "${TRIMMED_DIR}/${prefix}_2_trim.fq.gz" \
        2> "${SUMMARY_DIR}/${prefix}.txt" \
        -S "${BAM_DIR}/${prefix}.sam"
done

# Sort, fix mates, mark duplicates (samtools 1.10). Duplicates are marked, not
# removed; they are ignored in featureCounts and bamCoverage below
for f in $(ls "${BAM_DIR}"/*.sam)
do
    prefix=${f%.*}
    samtools sort -n "$f" -o "${prefix}_sname.bam" || exit 1
done
rm "${BAM_DIR}"/*.sam

for f in $(ls "${BAM_DIR}"/*_sname.bam)
do
    prefix=${f%_sname.*}
    samtools fixmate -m "$f" "${prefix}_fixmate.bam" || exit 1
done
rm "${BAM_DIR}"/*_sname.bam

for f in $(ls "${BAM_DIR}"/*_fixmate.bam)
do
    prefix=${f%_fixmate.*}
    samtools sort "$f" -o "${prefix}_sorted.bam" || exit 1
done
rm "${BAM_DIR}"/*_fixmate.bam

for f in $(ls "${BAM_DIR}"/*_sorted.bam)
do
    prefix=${f%_sorted.*}
    samtools markdup "$f" "${prefix}_markdup.bam" || exit 1
done
rm "${BAM_DIR}"/*_sorted.bam

for f in $(ls "${BAM_DIR}"/*_markdup.bam)
do
    prefix=${f%_markdup.*}
    samtools sort "$f" -o "${prefix}_final.bam" || exit 1
done
rm "${BAM_DIR}"/*_markdup.bam

# Merge the two lanes per sample, and the three replicates per genotype
LANEMERGED_DIR="${BAM_DIR}/LanesMerged_BAMs"
REPLICATEMERGED_DIR="${BAM_DIR}/ReplicatesMerged_BAMs"

mkdir -p "${LANEMERGED_DIR}"
mkdir -p "${REPLICATEMERGED_DIR}"

for sample in AO1 AO2 AO3 APP1 APP2 APP3 Ogg1 Ogg2 Ogg3 WT1 WT2 WT3
do
    samtools merge \
        "${LANEMERGED_DIR}/${sample}.bam" \
        "${BAM_DIR}/${sample}_L1_final.bam" \
        "${BAM_DIR}/${sample}_L2_final.bam"
done

for group in AO APP Ogg WT
do
    samtools merge \
        "${REPLICATEMERGED_DIR}/${group}.bam" \
        "${LANEMERGED_DIR}/${group}1.bam" \
        "${LANEMERGED_DIR}/${group}2.bam" \
        "${LANEMERGED_DIR}/${group}3.bam"
done

for f in "${LANEMERGED_DIR}"/*.bam
do
    samtools index ${f}
done

for f in "${REPLICATEMERGED_DIR}"/*.bam
do
    samtools index ${f}
done
```

## Peak calling

Peaks are called per sample with MACS3, merged per genotype and then across
genotypes into one consensus peak set.

```sh
# MACS3 (v3.0.0a7), bedtools
ANALYSIS_DIR="Analysis"
PEAK_DIR="${ANALYSIS_DIR}/PeakCalling_macs3"

mkdir -p "${PEAK_DIR}"

for f in "${LANEMERGED_DIR}"/*.bam
do
    prefix=$(basename "$f" .bam)
    macs3 callpeak -g mm -t "$f" -n "${PEAK_DIR}/${prefix}" --nomodel -f BAMPE -q 0.05
done

for group in AO APP Ogg WT
do
    cat \
        "${PEAK_DIR}/${group}1_peaks.narrowPeak" \
        "${PEAK_DIR}/${group}2_peaks.narrowPeak" \
        "${PEAK_DIR}/${group}3_peaks.narrowPeak" \
        > "${PEAK_DIR}/${group}Peaks.bed"
done

for group in AO APP Ogg WT
do
    bedtools sort \
        -i "${PEAK_DIR}/${group}Peaks.bed" \
        > "${PEAK_DIR}/${group}_sorted.bed"

    bedtools merge \
        -i "${PEAK_DIR}/${group}_sorted.bed" \
        > "${PEAK_DIR}/${group}_Peaks.bed"
done

cat "${PEAK_DIR}"/*_Peaks.bed \
    > "${PEAK_DIR}/AllGroupsPeaks.bed"

bedtools sort \
    -i "${PEAK_DIR}/AllGroupsPeaks.bed" \
    > "${PEAK_DIR}/AllGroups_sorted.bed"

bedtools merge \
    -i "${PEAK_DIR}/AllGroups_sorted.bed" \
    > "${PEAK_DIR}/Consensus_Peaks.bed"
```

### Peak annotation

Peaks are annotated to the nearest gene (TSS region +/- 2 kb) with ChIPseeker.
The Bowtie 2 index uses Ensembl chromosome names, so "chr" is added before
annotation against the UCSC mm10 TxDb.

```r
# R version 4.3.2
library(dplyr)
library(ChIPseeker)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicRanges)
library(openxlsx)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
PEAK_DIR <- "/path/to/ATAC_Analysis/Analysis/PeakCalling_macs3"

To_GRange <- function(x) {
  makeGRangesFromDataFrame(x,
    ignore.strand = TRUE,
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE)}

# Gene lists per genotype peak set and for the consensus peaks
annotate_peaks <- function(bedfile) {
  bed <- read.table(bedfile, header = FALSE)
  colnames(bed)[1:3] <- c("seqnames", "start", "end")
  bed$seqnames <- paste0("chr", bed$seqnames)
  gr <- To_GRange(bed)
  anno <- annotatePeak(
    gr,
    tssRegion = c(-2000, 2000),
    TxDb = txdb,
    annoDb = "org.Mm.eg.db")
  anno_df <- as.data.frame(anno)
  genes <- anno_df[, c(1, 8, 9, 16, 10, 11)]
  genes$geneStrand[genes$geneStrand == 1] <- "+"
  genes$geneStrand[genes$geneStrand == 2] <- "-"
  genes <- genes[!duplicated(genes), ]
  outfile <- sub("\\.bed$", "_Annotated_Genes.tsv", bedfile)
  write.table(genes,
    file = outfile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE)}

peak_files <- list.files(PEAK_DIR, pattern = "_Peaks\\.bed$", full.names = TRUE)
lapply(peak_files, annotate_peaks)

# BED file of the genes of the consensus peaks (used by computeMatrix below)
promoter <- getPromoters(TxDb=txdb, upstream=2000, downstream=2000)
bed <- read.table(file.path(PEAK_DIR,"Consensus_Peaks.bed"))
colnames(bed) <- c("seqnames","start","end")
bed$seqnames <- paste0("chr", bed$seqnames)

ranged <- makeGRangesFromDataFrame(bed, ignore.strand = T, seqnames.field = "seqnames", start.field = "start", end.field = "end")

anno <- annotatePeak(ranged, tssRegion=c(-2000, 2000),TxDb=txdb, annoDb="org.Mm.eg.db")
file <- as.data.frame(anno)

# chrom, start, end, name, score, strand
genes <- file[,c(7,8,9,16,10,11)]
genes$geneStrand[genes$geneStrand == 2] <- c('-')
genes$geneStrand[genes$geneStrand == 1] <- c('+')
gene <- genes[!duplicated(genes),]
write.table(gene, file.path(PEAK_DIR, "Genes_of_Consensus_Peaks.bed"), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
```

## BigWigs and profile plot (Suppl. Fig. 5b)

```sh
# deepTools (v3.5.2)
BW_DIR="${REPLICATEMERGED_DIR}/bamCoverage_RPKM"
BLACKLIST="/path/to/blacklist/mm10-blacklist.v2.bed.gz"   # Boyle lab: https://github.com/Boyle-Lab/Blacklist

mkdir -p "${BW_DIR}"

# RPKM-normalized coverage per genotype, for the TSS profile plot
for f in "${REPLICATEMERGED_DIR}"/*.bam
do
    prefix=$(basename "$f" .bam)
    bamCoverage -p "$THREADS" --effectiveGenomeSize 2652783500 -bs 50 --smoothLength 150 --normalizeUsing RPKM --ignoreDuplicates -b ${f} -o "$BW_DIR/${prefix}.bw"
done

BWs=("${BW_DIR}"/*.bw)

computeMatrix reference-point -p "$THREADS" -S "${BWs[@]}" -R "$PEAK_DIR/Genes_of_Consensus_Peaks.bed" --blackListFileName "$BLACKLIST" -b 1000 -a 1000 --missingDataAsZero --skipZeros -o "$BW_DIR/computeMatrix_ConsensusPeakGenes_blacklist_refpoint.mat.gz"

plotProfile -m "$BW_DIR/computeMatrix_ConsensusPeakGenes_blacklist_refpoint.mat.gz" --averageType mean --plotType se --numPlotsPerRow 2 --perGroup --colors '#f57a0c' '#056475' '#543c6d' '#5e5e5e' --plotTitle "ATAC-seq profile over TSS" --regionsLabel " " --yAxisLabel "Mean Read Coverage" --plotWidth 15 --plotHeight 15 -out "$BW_DIR/SE_Profile-over_ConsensusPeakGenes_blacklist_refpoint.svg"

# 1x-normalized (RPGC) coverage per genotype, for the genome tracks
BW1x_DIR="${REPLICATEMERGED_DIR}/bamCoverage_1xCov"
mkdir -p "${BW1x_DIR}"

for f in "${REPLICATEMERGED_DIR}"/*.bam
do
    prefix=$(basename "$f" .bam)
    bamCoverage -p "$THREADS" --effectiveGenomeSize 2652783500 -bs 200 --smoothLength 600 --normalizeUsing RPGC --ignoreDuplicates -b ${f} -o "$BW1x_DIR/${prefix}.bw"
done
```

## Counts

Reads per consensus peak and sample (lane-merged BAMs).

```sh
# featureCounts (Subread v2.0.0)
DE_DIR="${ANALYSIS_DIR}/DiffEnrich"
mkdir -p "${DE_DIR}"

echo -e "GeneID\tChr\tStart\tEnd\tStrand" > "${PEAK_DIR}/Consensus_Peaks.saf"
awk 'OFS="\t" {print $1"_"$2"_"$3,$1,$2,$3,"+"}' "${PEAK_DIR}/Consensus_Peaks.bed" >> "${PEAK_DIR}/Consensus_Peaks.saf"

featureCounts -T 4 -F SAF -M -Q 30 --ignoreDup -p -a "${PEAK_DIR}/Consensus_Peaks.saf" -o "${DE_DIR}/ATAC_Counts.txt" "${LANEMERGED_DIR}"/*.bam
```

## Differentially accessible regions (DARs) and PCA (Suppl. Fig. 5a)

edgeR quasi-likelihood F-tests per consensus peak, same five contrasts as for
the bulk RNA-seq. Significant peaks (`_FilterP05FC1`) are defined as
P <= 0.05 and |log2FC| >= 1. Each peak is annotated to its nearest gene with
ChIPseeker. The PCA uses edgeR log2 CPM of the filtered peaks.

```r
library(limma)
library(edgeR)
library(dplyr)
library(ChIPseeker)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicRanges)
library(openxlsx)
library(svglite)
library(PCAtools)

DE_DIR <- "/path/to/ATAC_Analysis/Analysis/DiffEnrich"
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Split the peak ID (chr_start_end) back into coordinates
colsplit <- function(x){
  x1 = mutate(x, chr = sapply(strsplit(as.character(rownames(x)),'_'), "[", 1))
  x3 = mutate(x1, start = sapply(strsplit(as.character(rownames(x)),'_'), "[", 2))
  x4 = mutate(x3, end = sapply(strsplit(as.character(rownames(x)),'_'), "[", 3))
  return(x4)}

To_GRange <- function(x){
  makeGRangesFromDataFrame(x, ignore.strand = T, seqnames.field = "chr", start.field = "start", end.field = "end",keep.extra.columns=TRUE)}

# Columns 7-18 are the 12 samples in alphabetical order: AO1-3, APP1-3, Ogg1-3, WT1-3
countData <- read.csv(file.path(DE_DIR, "ATAC_Counts.txt"), header=TRUE, sep="\t")[,c(1,7:18)]
rownames(countData) <- countData$Geneid
countData<-countData[,-1]

Genotype <- c(rep("AO", each=3), rep("APP/PS1", each=3), rep("Ogg1-/-", each=3), rep("WT", each=3))
Replicate <- rep(1:3, times = 4)

d <- DGEList(counts=countData,lib.size=colSums(countData), group=factor(Genotype))
design <- model.matrix(~0+group, data=d$samples)
colnames(design) <- levels(d$samples$group)

keep <- filterByExpr(d, factor(Genotype))
dge <- d[keep,]
dge <- calcNormFactors(dge, method="TMM")
dge <- estimateDisp(dge,design)
fit <- glmQLFit(dge, design)

# PCA
cpms<- cpm(dge, log=TRUE)
meta = data.frame(Genotype, Replicate)
rownames(meta) <- colnames(countData)

p <- pca(cpms, metadata = meta, removeVar = 0.1)
svglite(file.path(DE_DIR, "ATAC_PCAtools_EdgeRcpm.svg"), width = 8, height = 6)
biplot(p, pointSize=5, colby='Genotype', colkey = c('WT' = '#5e5e5e', 'APP/PS1' = '#056475', 'Ogg1-/-' = '#543c6d', 'AO' = '#f57a0c'), lab=p$metadata[,2], legendPosition='right', title='PCA of ATAC-sequencing', caption='EdgeR log2 CPM')
dev.off()

# Design columns (alphabetical factor levels): AO, APP/PS1, Ogg1-/-, WT
contrasts <- list(
  AO   = c(1, 0, 0,-1),
  APP  = c(0, 1, 0,-1),
  OGG  = c(0, 0, 1,-1),
  AOAP = c(1,-1, 0, 0),
  AOOG = c(1, 0,-1, 0))

results <- lapply(contrasts, function(x) glmQLFTest(fit, contrast = x))
results <- lapply(results, function(x) topTags(x, n = Inf)$table)
results <- lapply(results, colsplit)
results <- lapply(results, To_GRange)
results <- lapply(results, function(x) annotatePeak(x, tssRegion = c(-2000, 2000), TxDb = txdb, annoDb = "org.Mm.eg.db"))

results <- lapply(results, as.data.frame)
results <- lapply(results, function(x) x[order(x$logFC), ])

sig_results <- lapply(results, function(x) {x %>% filter(PValue <= 0.05) %>% filter(abs(logFC) >= 1)})

compnames <- c("AOvsWT", "APPvsWT", "OGGvsWT", "AOvsAPP", "AOvsOGG")
names(results) <- compnames
names(sig_results) <- compnames

write.xlsx(results, file.path(DE_DIR, "EdgeRqlf_featureCounts_DEresults.xlsx"), rowNames=TRUE, colNames=TRUE)
write.xlsx(sig_results, file.path(DE_DIR, "EdgeRqlf_featureCounts_DEresults_FilterP05FC1.xlsx"), rowNames=TRUE, colNames=TRUE)

sink(file.path(DE_DIR, "sessionInfo_ATAC-Seq-DE.txt"), append = FALSE)
print(sessioninfo::session_info())
sink()
```

## Prevented DARs (prevDARs; Suppl. Fig. 5f, Supplementary Table 1)

Peaks significant in both APP/PS1 vs WT and AO vs APP/PS1 (P <= 0.05,
|log2FC| >= 1). The compensation is -log2FC(AO vs APP/PS1) / log2FC(APP/PS1
vs WT); peaks with a positive value are the prevDARs. `resDARs.xlsx` is
written before this filter. The correlation plot drawn here is the draft of
Suppl. Fig. 5f; the published panel is made from `prevDARs_ATAC.txt` with
`R_prevDEG_correlation_plots.R` in this repository.

```r
library(dplyr)
library(openxlsx)
library(ggplot2)
library(ggrepel)
library(svglite)

PEAK_DIR <- "/path/to/ATAC_Analysis/Analysis/PeakCalling_macs3"
DE_DIR <- "/path/to/ATAC_Analysis/Analysis/DiffEnrich"

DE_file <- file.path(DE_DIR, "EdgeRqlf_featureCounts_DEresults_FilterP05FC1.xlsx")
sheet_names <- getSheetNames(DE_file)
DE_results <- lapply(sheet_names, function(sheet) {read.xlsx(DE_file, sheet = sheet)})
names(DE_results) <- sheet_names

cols <- c("seqnames", "start", "end", "SYMBOL", "logFC", "PValue", "FDR")
appWT <- DE_results[[2]][, cols]
aoAPP <- DE_results[[4]][, cols]

resDARs <- inner_join(appWT, aoAPP, by = c("seqnames", "start", "end"), suffix = c("_APP-vs-WT", "_AO-vs-APP"))
resDARs$PercRescue <- -resDARs[,9]/resDARs[,5]
write.xlsx(resDARs, file = file.path(DE_DIR, "resDARs.xlsx"), rowNames = FALSE)
resDARs <- resDARs %>% filter(PercRescue>0)
resDARs$strand <- "*"
write.table(resDARs[,c(1:4,12:13)], file = file.path(PEAK_DIR, "resPEAKs.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Correlation plot
resDARs$Score <- resDARs$PercRescue*100

genes_to_label <- c( "Firre","Nwd2","Ntrk2","Fmr1","Cbln2","Ppp2r2b","Alk","Hepacam")
df_labels <- resDARs[resDARs[,4] %in% genes_to_label, ]

min_val <- floor(min(resDARs$Score) / 10) * 10
max_val <- ceiling(max(resDARs$Score) / 10) * 10
mid_val <- (min_val + max_val) / 2
lims <- ceiling(max(abs(c(resDARs[,9], resDARs[,5]))))

atac <- ggplot(resDARs, aes(x = resDARs[,9], y = resDARs[,5])) +
  geom_abline(slope = -1, intercept = 0, color = "grey50", linewidth = 0.42, linetype = "solid") +
  geom_point(aes(color = Score), shape = 16, size = 4.5, alpha = 0.8) +
  geom_point(data = df_labels, aes(x = `logFC_AO-vs-APP`, y = `logFC_APP-vs-WT`), inherit.aes = FALSE, shape = 21, color = "black", size = 5, stroke = 1.5) +
  geom_text_repel(data = df_labels, aes(x = `logFC_AO-vs-APP`, y = `logFC_APP-vs-WT`, label = `SYMBOL_APP-vs-WT`),
  inherit.aes = FALSE, size = 6, max.overlaps = Inf, box.padding = 0.6, force = 3) +
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

svglite(file.path(DE_DIR, "ATAC_resDAR_CorrelationPlot.svg"), width = 8, height = 8)
atac
dev.off()
```

## Genome tracks (Fig. 2f)

Coverage at the Firre locus (chrX:50,555,000-50,636,000, mm10) per genotype,
with the consensus peaks and prevDAR peaks. `UcscTrack()` downloads the gene
model from UCSC and needs internet access. The bigWigs use Ensembl chromosome
names ("X"), the UCSC gene track uses "chrX".

```r
library(Gviz)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(openxlsx)
library(svglite)
library(grid)

PEAK_DIR <- "/path/to/ATAC_Analysis/Analysis/PeakCalling_macs3"
BW_DIR <- "/path/to/ATAC_Analysis/fastq/trimmed_fastq/BAMs/ReplicatesMerged_BAMs/bamCoverage_1xCov"

ConsPeak_file <- file.path(PEAK_DIR, "Consensus_Peaks.bed")
ResPeak_file <- file.path(PEAK_DIR, "resPEAKs.bed")

# Firre: chrX 50555744-50635261
chromosome <- "chrX"
genome <- "mm10"
start <- 50555000
end <- 50636000

geneTrack <- UcscTrack(genome=genome, chromosome=chromosome, track="NCBI RefSeq",
                       table="refGene", from=start, to=end, trackType="GeneRegionTrack",
                       rstarts="exonStarts", rends="exonEnds",
                       gene="name", symbol="name2", transcript=NULL,
                       name="Genes")

gr <- ranges(geneTrack)
gr_firre <- gr[symbol(geneTrack) == "Firre"]
Firre <- GeneRegionTrack(gr_firre, genome=genome, chromosome=chromosome,
                         name="Firre", background.title="transparent",
                         col.title="black",cex.title = 1.1, fill="#095eb4", col="#095eb4",
                         col.line="black",collapseTranscripts="longest")

save(Firre, file="Firre_track.RData")

peaks <- import(ConsPeak_file, format="BED")
chrx <- peaks[seqnames(peaks) == "chrX"]
Peaks <- AnnotationTrack(chrx, genome=genome, chromosome=chromosome,
                         fill="#1E905F", col="#1E905F",
                         name="Peaks", background.title="transparent",
                         col.title="black", cex.title = 1.1)

respeaks <- import(ResPeak_file, format="BED")
resPeaks <- AnnotationTrack(respeaks, genome=genome, chromosome=chromosome,
                          fill="#1E709F", col="#1E709F",
                          name="resPeaks", background.title="transparent",
                          col.title="black", cex.title = 1.1)

load("Firre_track.RData")
chromosome <- "X"
genome <- "mm10"
start <- 50555000
end <- 50636000
low <- 0
lim <- 20

gtrack <- GenomeAxisTrack(range = GRanges(chromosome, IRanges(start, end)),
                          name = "", fill="#000000", col="#000000",
                          fontcolor="#000000", fontsize=15)

AO_bw <- import(file.path(BW_DIR, "AO.bw"), which=GRanges(chromosome, IRanges(start, end)))
AOTrack <- DataTrack(range=AO_bw, genome=genome, chromosome=chromosome,
                     name="AO", type="histogram", col.histogram="#f57a0c",
                     fill.histogram="#f57a0c", ylim=c(0, lim), showAxis=FALSE, cex.title = 1,
                     background.title="transparent", col.title="black", col.axis="black")

APP_bw <- import(file.path(BW_DIR,"APP.bw"), which=GRanges(chromosome, IRanges(start, end)))
APPTrack <- DataTrack(range=APP_bw, genome=genome, chromosome=chromosome,
                     name="APP/PS1", type="histogram", col.histogram="#056475",
                     fill.histogram="#056475", ylim=c(0, lim), showAxis=FALSE, cex.title = 1,
                     background.title="transparent", col.title="black", col.axis="black")

WT_bw <- import(file.path(BW_DIR,"WT.bw"), which=GRanges(chromosome, IRanges(start, end)))
WTTrack <- DataTrack(range=WT_bw, genome=genome, chromosome=chromosome,
                     name="WT", feature ="[0-20]", type="histogram", col.histogram="#5e5e5e",
                     fill.histogram="#5e5e5e",ylim=c(0, lim), showAxis=FALSE, cex.title = 1,
                     background.title="transparent", col.title="black", col.axis="black")

Ogg_bw <- import(file.path(BW_DIR,"Ogg.bw"), which=GRanges(chromosome, IRanges(start, end)))
OggTrack <- DataTrack(range=Ogg_bw, genome=genome, chromosome=chromosome,
                     name="Ogg1-/-", type="histogram", col.histogram="#543c6d",
                     fill.histogram="#543c6d",ylim=c(0, lim), showAxis=FALSE, cex.title = 1,
                     background.title="transparent", col.title="black", col.axis="black")

svglite(file.path(BW_DIR, "ATAC_Firre_GenomeTrack_1XCov_UCSC.svg"), width = 8, height = 8)
plotTracks(list(Firre, gtrack, Peaks, resPeaks, WTTrack, OggTrack, APPTrack, AOTrack), from=start, to=end, chromosome=chromosome)
grid.text(label = "Chr X",
          x = unit(0.08, "npc"),
          y = unit(0.98, "npc"),
          just = c("left", "top"),
          gp = gpar(col = "black", fontsize = 10))

grid.text(label = paste0("[", low, "-", lim, "]"),
          x = unit(0.05, "npc"),
          y = unit(0.75, "npc"),
          just = c("left", "top"),
          gp = gpar(col = "black", fontsize = 9))
grid.text(label = paste0("[", low, "-", lim, "]"),
          x = unit(0.05, "npc"),
          y = unit(0.55, "npc"),
          just = c("left", "top"),
          gp = gpar(col = "black", fontsize = 9))
grid.text(label = paste0("[", low, "-", lim, "]"),
          x = unit(0.05, "npc"),
          y = unit(0.35, "npc"),
          just = c("left", "top"),
          gp = gpar(col = "black", fontsize = 9))
grid.text(label = paste0("[", low, "-", lim, "]"),
          x = unit(0.05, "npc"),
          y = unit(0.15, "npc"),
          just = c("left", "top"),
          gp = gpar(col = "black", fontsize = 9))
dev.off()
```

## HOMER motifs (Fig. 2d, Suppl. Fig. 5d)

BED files of the significant peaks per genotype vs WT (all, and unique to
that genotype), with all tested peaks as background. Motif enrichment with
HOMER on the unique peaks.

```r
library(openxlsx)

DE_DIR <- "/path/to/ATAC_Analysis/Analysis/DiffEnrich"
DE_file <- file.path(DE_DIR, "EdgeRqlf_featureCounts_DEresults_FilterP05FC1.xlsx")
background <- file.path(DE_DIR, "EdgeRqlf_featureCounts_DEresults.xlsx")
OUT_DIR <- file.path(DE_DIR, "HOMER_All")
OUTUNIQUE_DIR <- file.path(DE_DIR, "HOMER_Unique")
dir.create(OUT_DIR, recursive = TRUE)
dir.create(OUTUNIQUE_DIR, recursive = TRUE)

cols <- c("seqnames", "start", "end", "SYMBOL", "logCPM", "strand")
df <- read.xlsx(background, sheet = 1)
df_filtered <- df[, cols]
write.table(df_filtered, file = file.path(OUT_DIR, "ATAC_Background.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(df_filtered, file = file.path(OUTUNIQUE_DIR, "ATAC_Background.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

sheet_names <- getSheetNames(DE_file)
DE_results <- lapply(sheet_names, function(sheet) {read.xlsx(DE_file, sheet = sheet)})
names(DE_results) <- sheet_names
DE_filtered <- lapply(DE_results[1:3], function(df) {df[, cols]})

lapply(names(DE_filtered), function(nm) {
  write.table(DE_filtered[[nm]],
    file = file.path(OUT_DIR, paste0(nm, ".bed")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)})

# Peaks significant in only one of the three genotypes
DE_unique <- lapply(seq_along(DE_filtered), function(i) {
  other_dfs <- do.call(rbind, DE_filtered[-i])
  df <- DE_filtered[[i]]
  df[!paste(df[, 1], df[, 2], df[, 3]) %in%
       paste(other_dfs[, 1], other_dfs[, 2], other_dfs[, 3]), ]})

names(DE_unique) <- names(DE_filtered)
lapply(names(DE_unique), function(nm) {
  write.table(DE_unique[[nm]],
    file = file.path(OUTUNIQUE_DIR, paste0(nm, ".bed")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)})
```

```sh
# HOMER
HOMER_DIR="${DE_DIR}/HOMER_Unique"

for f in "${HOMER_DIR}"/*vsWT.bed
do
    prefix=$(basename "$f" .bed)
    perl findMotifsGenome.pl "$f" mm10 "${HOMER_DIR}/Motifs/${prefix}" -p 10 -bg "${HOMER_DIR}/ATAC_Background.bed" -S 10
done
```
