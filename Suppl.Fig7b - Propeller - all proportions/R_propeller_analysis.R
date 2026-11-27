# ------------------------------------------------------------------------------
# Suppl. Fig. 7b - Cell type composition analysis with propeller (speckle)
#
# Tests whether the proportion of each cell type in the snRNA-seq data differs
# between the four genotypes (WT, Ogg1, APP/PS1, AO), overall and pairwise.
#
# Input:  cell_type_proportions.csv (rows = cell types, columns = samples,
#                                    values = % of cells per sample)
# Output: propeller_anova.csv       (overall test across the 4 genotypes)
#         propeller_pairwise.csv    (all 6 genotype pairs, per cell type)
#         propeller_statistics.txt  (tab-delimited statistics summary: test
#                                    details, test statistics, degrees of
#                                    freedom, P values, FDR, data summary)
#         cell_type_proportions.pdf (% of cells per genotype, one panel
#                                    per cell type)
#
# speckle is a Bioconductor package:
#   install.packages("BiocManager")
#   BiocManager::install(c("speckle", "limma"))
# ------------------------------------------------------------------------------

# Load libraries
library(speckle)
library(limma)

# Settings
input_file <- "cell_type_proportions.csv"
transform  <- "asin"   # "asin" (arcsine square root) or "logit"
robust     <- TRUE     # robust empirical Bayes variance moderation
trend      <- FALSE

# Approximate number of cells per sample, only used to back-calculate counts.
# It does not change the results with transform = "asin"; it does with "logit",
# so set the real value (one number or one value per sample) if you use logit
total_cells <- 5000

# Text written to the statistics file
figure_label <- "Supplementary Fig. 7b - snRNA-seq cell type proportions"
sample_note  <- "3 independent biological replicates (each pooled from the hippocampi of 3 mice)"

# Genotype colors
genotype_colours <- c(WT = "#8A8A8A", Ogg1 = "#543c6c", APPPS1 = "#056475", AO = "#F57A0C")

# Display names used in tables and plots (internal names must stay free of "/")
genotype_labels <- c(WT = "WT", Ogg1 = "Ogg1", APPPS1 = "APP/PS1", AO = "AO")

# Read data
pct <- as.matrix(read.csv(input_file, row.names = 1, check.names = FALSE))
cat("Column sums (should be ~100):\n"); print(round(colSums(pct), 2))

# Convert % to proportions (0-1)
props <- sweep(pct, 2, colSums(pct), "/")

# Genotype = sample name without the replicate suffix (WT_1 -> WT)
group <- factor(sub("_[0-9]+$", "", colnames(props)),
                levels = c("WT", "Ogg1", "APPPS1", "AO"))
print(data.frame(sample = colnames(props), group = group))
stopifnot(!anyNA(group))

# propeller
prop.list <- convertDataToList(props, data.type = "proportions",
                               transform = transform, scale.fac = total_cells)

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Overall test: does the proportion differ between any of the 4 genotypes?
res_anova <- propeller.anova(prop.list, design = design, coef = seq_len(ncol(design)),
                             robust = robust, trend = trend, sort = TRUE)

# Degrees of freedom of the overall test. propeller does not return them, so
# the same limma model that propeller.anova fits internally is fitted here
# (intercept design, F-test on the 3 genotype coefficients) and the moderated
# degrees of freedom are read from it. The check stops the script if the
# F statistics of this fit and propeller disagree.
design_int <- design
design_int[, 1] <- 1
fit_F <- eBayes(lmFit(prop.list$TransformedProps, design_int)[, 2:ncol(design_int)],
                robust = robust, trend = trend)
idx <- match(rownames(res_anova), rownames(prop.list$TransformedProps))
stopifnot(!anyNA(idx),
          isTRUE(all.equal(as.numeric(fit_F$F)[idx], res_anova$Fstatistic)))
res_anova$DFn <- ncol(design) - 1
res_anova$DFd <- as.numeric(fit_F$df.total)[idx]
res_anova <- res_anova[, c(grep("^PropMean", names(res_anova), value = TRUE),
                           "Fstatistic", "DFn", "DFd", "P.Value", "FDR")]

cat("\n=== Overall test (4 genotypes) ===\n"); print(res_anova)
write.csv(res_anova, "propeller_anova.csv")

# Pairwise tests: every genotype pair
pairs <- combn(levels(group), 2, simplify = FALSE)
res_pairs <- do.call(rbind, lapply(pairs, function(p) {
  contr <- makeContrasts(contrasts = paste(p[1], "-", p[2]), levels = design)
  res <- propeller.ttest(prop.list, design = design, contrasts = contr[, 1],
                         robust = robust, trend = trend, sort = FALSE)

  # Degrees of freedom: same limma fit as inside propeller.ttest, with a check
  # that the t statistics agree
  fit_t <- eBayes(contrasts.fit(lmFit(prop.list$TransformedProps, design),
                                contrasts = contr[, 1]),
                  robust = robust, trend = trend)
  stopifnot(isTRUE(all.equal(as.numeric(fit_t$t[, 1]), res$Tstatistic)))

  ct <- rownames(res)
  data.frame(comparison   = paste(genotype_labels[p[1]], "vs", genotype_labels[p[2]]),
             cell_type    = ct,
             mean_pct_1   = rowMeans(pct[ct, group == p[1], drop = FALSE]),
             mean_pct_2   = rowMeans(pct[ct, group == p[2], drop = FALSE]),
             PropRatio    = res$PropRatio,
             Tstatistic   = res$Tstatistic,
             DF           = as.numeric(fit_t$df.total),
             P.Value      = res$P.Value,
             FDR          = res$FDR,
             row.names = NULL)
}))

# The FDR from propeller corrects across cell types within one comparison.
# This adds a stricter correction across all cell types x all comparisons
res_pairs$FDR_all_comparisons <- p.adjust(res_pairs$P.Value, method = "BH")

# Keep the comparison order for the statistics file, then sort by P value
comparison_order <- unique(res_pairs$comparison)
res_pairs <- res_pairs[order(res_pairs$P.Value), ]

cat("\n=== Pairwise tests (top 15 by P value) ===\n")
print(head(res_pairs, 15), row.names = FALSE)
write.csv(res_pairs, "propeller_pairwise.csv", row.names = FALSE)

# Statistics summary (tab-delimited text, decimal point)
f3    <- function(x) formatC(x, format = "f", digits = 3)
f2    <- function(x) formatC(x, format = "f", digits = 2)
fp    <- function(x) format(signif(x, 4), scientific = FALSE, trim = TRUE)
stars <- function(p) ifelse(p < 0.0001, "****", ifelse(p < 0.001, "***",
                                                       ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
yesno <- function(x) ifelse(x, "Yes", "No")
tab   <- function(...) paste(..., sep = "\t")

lev    <- levels(group)
n_ct   <- nrow(pct)
n_test <- nrow(res_pairs)
out    <- c()

out <- c(out,
         "Cell type composition analysis (propeller)",
         figure_label,
         "",
         tab("Table analyzed", paste0(input_file, " (% of nuclei per cell type and sample)")),
         tab("Data sets analyzed", paste(genotype_labels[lev], collapse = ", ")),
         tab("Number of groups (genotypes)", length(lev)),
         tab("Number of samples (total)", ncol(pct)),
         tab("Number of samples per genotype", paste0(paste(table(group), collapse = ", "), "; ", sample_note)),
         tab("Number of cell types tested", n_ct),
         "",
         "Test details",
         tab("Software", paste0("propeller (speckle v", packageVersion("speckle"),
                                ") with limma v", packageVersion("limma"), "; ", R.version.string)),
         tab("Input", "Proportions per sample (data.type = \"proportions\")"),
         tab("Transformation", ifelse(transform == "asin", "Arcsine square root (transform = \"asin\")",
                                      "Logit (transform = \"logit\")")),
         tab("Model", "Linear model per cell type, one coefficient per genotype (~ 0 + genotype)"),
         tab("Variance estimation", paste0("Empirical Bayes moderation, robust = ", robust, ", trend = ", trend)),
         tab("Overall test", paste0("Moderated F-test across the ", length(lev), " genotypes (propeller.anova)")),
         tab("Pairwise tests", paste0("Moderated t-tests for all ", length(pairs), " genotype pairs (propeller.ttest), two-sided")),
         tab("Multiple testing correction", "Benjamini-Hochberg FDR"),
         tab("Significance threshold", "FDR < 0.05"),
         tab("Degrees of freedom", "Moderated (residual DF + prior DF), taken from the limma fit (df.total)"),
         "")

# Overall test
out <- c(out,
         paste0("Overall test (moderated F-test, ", length(lev), " genotypes)"),
         tab("Cell type", paste(paste("Mean %", genotype_labels[lev]), collapse = "\t"),
             "F (DFn, DFd)", "P value", "P value summary",
             paste0("FDR (across ", n_ct, " cell types)"), "FDR summary", "Significant (FDR < 0.05)?"))
for (ct in rownames(res_anova)) {
  r <- res_anova[ct, ]
  out <- c(out, tab(ct,
                    paste(f3(100 * unlist(r[paste0("PropMean.", lev)])), collapse = "\t"),
                    paste0("F (", r$DFn, ", ", f2(r$DFd), ") = ", f3(r$Fstatistic)),
                    fp(r$P.Value), stars(r$P.Value), fp(r$FDR), stars(r$FDR), yesno(r$FDR < 0.05)))
}
out <- c(out, "")

# Pairwise tests, grouped by comparison and sorted by P value within each
out <- c(out,
         "Pairwise comparisons (moderated t-tests, two-sided)",
         tab("Comparison", "Cell type", "Mean % group 1", "Mean % group 2", "Mean diff. (% points)",
             "Proportion ratio (group 1 / group 2)", "t", "DF", "P value", "P value summary",
             paste0("FDR (across ", n_ct, " cell types, within comparison)"),
             paste0("FDR (across all ", n_test, " tests)"),
             "Significant (FDR within comparison < 0.05)?",
             "Significant (FDR across all tests < 0.05)?"))
for (cmp in comparison_order) {
  s <- res_pairs[res_pairs$comparison == cmp, ]
  out <- c(out, tab(s$comparison, s$cell_type, f3(s$mean_pct_1), f3(s$mean_pct_2),
                    f3(s$mean_pct_1 - s$mean_pct_2), f3(s$PropRatio), f3(s$Tstatistic),
                    f2(s$DF), fp(s$P.Value), stars(s$P.Value), fp(s$FDR),
                    fp(s$FDR_all_comparisons), yesno(s$FDR < 0.05),
                    yesno(s$FDR_all_comparisons < 0.05)))
}
out <- c(out, "")

# Data summary
out <- c(out,
         "Data summary (% of nuclei)",
         tab("Cell type", paste(as.vector(t(outer(genotype_labels[lev], c("mean", "SD", "SEM", "n"),
                                                  paste))), collapse = "\t")))
for (ct in rownames(pct)) {
  cells <- unlist(lapply(lev, function(g) {
    v <- pct[ct, group == g]
    c(f3(mean(v)), f3(sd(v)), f3(sd(v) / sqrt(length(v))), length(v))
  }))
  out <- c(out, tab(ct, paste(cells, collapse = "\t")))
}
out <- c(out, "")

# Individual values
out <- c(out,
         "Individual values (% of nuclei per sample)",
         tab("Cell type", paste(colnames(pct), collapse = "\t")))
for (ct in rownames(pct)) {
  out <- c(out, tab(ct, paste(f3(pct[ct, ]), collapse = "\t")))
}

writeLines(out, "propeller_statistics.txt")
cat("\nStatistics summary written to propeller_statistics.txt\n")

# Plot: % of cells per genotype, one panel per cell type, with group means
pdf("cell_type_proportions.pdf", width = 10, height = 7)
par(mfrow = c(3, 4), mar = c(5.5, 4, 2, 1))
for (ct in rownames(pct)) {
  stripchart(pct[ct, ] ~ group, vertical = TRUE, method = "jitter", pch = 16,
             col = genotype_colours[levels(group)],
             group.names = genotype_labels[levels(group)],
             main = ct, ylab = "% of cells", las = 2, cex.main = 0.9)
  means <- tapply(pct[ct, ], group, mean)
  segments(seq_along(means) - 0.25, means, seq_along(means) + 0.25, means, lwd = 2)
}
dev.off()

sessionInfo()
