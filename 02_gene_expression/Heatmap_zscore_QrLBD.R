###############################################################################
#  Heatmap_zscore_QrLBD.R
#  --------------------------------------------------------------------------
#  Builds the publication-quality z-score expression heatmap for the
#  primary QrLBD genes across the six Q. rubra tissues from PRJNA273270.
#
#  Pipeline (matches the user's DESeq2 prefiltering convention from
#  4_DESeq.r:  dds <- dds[ rowSums(counts(dds)) > 20, ])
#  --------
#     1. Load the FULL StringTie / prepDE.py3 gene_count_matrix.csv
#        (~33k genes x 11 SRA libraries).
#     2. Pre-filter — drop any gene with rowSums(counts) <= 20 across the
#        11 libraries. This is performed on the WHOLE matrix, not on the
#        LBD subset, so the library-size normalisation in step 3 has a
#        stable basis.
#     3. Library-size normalisation — per-library size factor =
#        column-total / median(column-total); divide each cell by its
#        column factor (median-of-ratios style; equivalent to DESeq2's
#        estimateSizeFactors() and avoids the DESeq2 dependency).
#     4. Subset to the 42 primary QrLBD genes via
#        QrLBD_transcript_gene_mapping.csv (status == "primary"). QrLBDs
#        that did not survive the rowSums > 20 filter (because they are
#        essentially silent across the 6 tissues at the depth of this
#        dataset) are reported up front and dropped from the heatmap.
#     5. Average across biological replicates per tissue (Mother + Father
#        libraries within the same tissue are averaged).
#     6. Per-gene z-score across the 6 tissue averages, capped at +-2 to
#        prevent extreme outliers from saturating the colour ramp.
#     7. Re-order rows by phylogeny — the chapter-era order is remapped
#        to the new 42-gene names; the two old entries that became
#        sub-isoforms in the new catalogue (old QrLBD25 -> QrLBD19.1,
#        old QrLBD27 -> QrLBD19.2) are dropped, leaving up to 42 ordered
#        rows. Genes filtered out in step 2/4 are simply skipped.
#     8. Render with pheatmap, RdYlBu reversed palette, breaks [-2, 2],
#        cluster_rows = FALSE / cluster_cols = FALSE (so the manual
#        phylogeny order is preserved).
#
#  Run from the Quercus_rubra/ folder:
#       Rscript Heatmap_zscore_QrLBD.R
#  or in RStudio:  setwd() to this folder, then source().
#
#  Required packages: pheatmap, RColorBrewer
###############################################################################

# ---------------------------------------------------------------------------
# 0. Setup
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(pheatmap)
  library(RColorBrewer)
})

WORK_DIR  <- "."
COUNTS    <- "../data/QrLBD_counts_42_primary.csv"             # raw StringTie / prepDE output
MAPPING   <- "QrLBD_transcript_gene_mapping.csv"    # gene_id <-> QrLBD name
OUT_PNG   <- "QrLBD_zscore_heatmap.png"
OUT_PDF   <- "QrLBD_zscore_heatmap.pdf"
OUT_MAT   <- "QrLBD_zscore_matrix.csv"

PREFILTER_MIN_TOTAL <- 20    # matches DESeq2 prefilter (4_DESeq.r line 84)

cat("==========================================================\n")
cat("  QrLBD z-score expression heatmap\n")
cat("  prefilter:  rowSums(counts) >", PREFILTER_MIN_TOTAL,
    "  (matches DESeq2 convention)\n")
cat("==========================================================\n")
cat("  working dir : ", normalizePath(WORK_DIR), "\n")
cat("  counts in   : ", COUNTS, "\n")
cat("  mapping     : ", MAPPING, "\n")
cat("  output PNG  : ", OUT_PNG, "\n\n")

stopifnot(file.exists(COUNTS))
stopifnot(file.exists(MAPPING))

# ---------------------------------------------------------------------------
# 1. Sample metadata (PRJNA273270; Konar et al., 2017)
# ---------------------------------------------------------------------------
sample_meta <- list(
  SRR1778917 = c("Emerging leaf/bud",  "Mother"),
  SRR1778925 = c("Emerging leaf/bud",  "Father"),
  SRR1778923 = c("Immature twig",      "Mother"),
  SRR1778931 = c("Immature twig",      "Father"),
  SRR1778921 = c("Twig (late growth)", "Mother"),
  SRR1778929 = c("Twig (late growth)", "Father"),
  SRR1778919 = c("Mature leaf",        "Mother"),
  SRR1778927 = c("Mature leaf",        "Father"),
  SRR1778932 = c("Dormant twig",       "Father"),
  SRR1778922 = c("Developing acorn",   "Mother"),
  SRR1778930 = c("Developing acorn",   "Father")
)
tissue_order <- c("Emerging leaf/bud", "Immature twig", "Twig (late growth)",
                  "Mature leaf", "Dormant twig", "Developing acorn")

# ---------------------------------------------------------------------------
# 2. Load FULL count matrix  (33k genes x 11 samples)
# ---------------------------------------------------------------------------
counts_full <- as.matrix(read.csv(COUNTS, row.names = 1, check.names = FALSE))
mode(counts_full) <- "numeric"
cat("  full count matrix:", nrow(counts_full), "genes x ",
    ncol(counts_full), "samples\n")

# Sanity: every count column should have a sample-meta entry
miss_meta <- setdiff(colnames(counts_full), names(sample_meta))
if (length(miss_meta)) {
  warning("Samples without metadata (will be dropped): ",
          paste(miss_meta, collapse = ", "))
  counts_full <- counts_full[, !colnames(counts_full) %in% miss_meta,
                             drop = FALSE]
}

# ---------------------------------------------------------------------------
# 3. Pre-filter on the FULL matrix  (DESeq2-style)
# ---------------------------------------------------------------------------
keep <- rowSums(counts_full) > PREFILTER_MIN_TOTAL
cat("  prefilter:  rowSums >", PREFILTER_MIN_TOTAL, "  ->  kept",
    sum(keep), "/", length(keep), " genes\n")
counts_filt <- counts_full[keep, , drop = FALSE]

# ---------------------------------------------------------------------------
# 4. Library-size normalisation on the filtered matrix
#    (median-of-ratios style; equivalent in spirit to DESeq2's
#     estimateSizeFactors() but with no Bioconductor dependency)
# ---------------------------------------------------------------------------
col_tot   <- colSums(counts_filt)
size_fact <- col_tot / median(col_tot)
norm_full <- sweep(counts_filt, 2, size_fact, FUN = "/")
cat("  size factors (per library):\n")
print(round(size_fact, 3))

# ---------------------------------------------------------------------------
# 5. Subset to the 42 primary QrLBDs  (those that survive the prefilter)
# ---------------------------------------------------------------------------
mapping_all <- read.csv(MAPPING, stringsAsFactors = FALSE)
mapping     <- mapping_all[mapping_all$status == "primary", ]

qr_in_full   <- intersect(mapping$gene_id, rownames(counts_full))
qr_in_filt   <- intersect(mapping$gene_id, rownames(counts_filt))
qr_dropped   <- setdiff(qr_in_full, qr_in_filt)

cat("\n  42 primary QrLBDs:\n")
cat("     present in count matrix     :", length(qr_in_full), "/ 42\n")
cat("     surviving rowSums >",  PREFILTER_MIN_TOTAL, "       :",
    length(qr_in_filt), "/ 42\n")
cat("     dropped at prefilter step   :", length(qr_dropped), "\n")
if (length(qr_dropped)) {
  drop_names <- mapping$QrLBDs_new[match(qr_dropped, mapping$gene_id)]
  drop_sums  <- rowSums(counts_full[qr_dropped, , drop = FALSE])
  cat("     (genes silent or near-silent across the 6 tissues)\n")
  for (i in seq_along(qr_dropped)) {
    cat(sprintf("        %-9s rowSum = %d\n",
                drop_names[i], as.integer(drop_sums[i])))
  }
}

lbd_norm <- norm_full[qr_in_filt, , drop = FALSE]
gene_to_name <- setNames(mapping$QrLBDs_new, mapping$gene_id)
rownames(lbd_norm) <- gene_to_name[rownames(lbd_norm)]

# ---------------------------------------------------------------------------
# 6. Average across replicates by tissue
# ---------------------------------------------------------------------------
tissue_per_sample <- vapply(colnames(lbd_norm),
                            function(s) sample_meta[[s]][1],
                            character(1))

unique_tissues <- intersect(tissue_order, unique(tissue_per_sample))
tissue_avg <- sapply(unique_tissues, function(t) {
  cols <- which(tissue_per_sample == t)
  if (length(cols) == 1) lbd_norm[, cols] else rowMeans(lbd_norm[, cols, drop = FALSE])
})
tissue_avg <- as.matrix(tissue_avg)
colnames(tissue_avg) <- unique_tissues
cat("\n  tissue averages: ", nrow(tissue_avg), " genes x ",
    ncol(tissue_avg), " tissues\n")

# ---------------------------------------------------------------------------
# 7. Per-gene z-score, capped at +-2
# ---------------------------------------------------------------------------
z_score <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}
qr_z <- t(apply(tissue_avg, 1, z_score))
qr_z <- pmax(pmin(qr_z, 2), -2)
colnames(qr_z) <- colnames(tissue_avg)

# ---------------------------------------------------------------------------
# 8. Phylogeny order — NEW iTOL tree leaf order (42 primary genes)
#    Top-to-bottom from the user's updated phylogeny with subclass colour
#    strip:  Class Ie -> Ic -> Ia -> Id -> Ib -> Ia -> IIb -> IIa
# ---------------------------------------------------------------------------
phylo_qr_new <- c(
  # Class Ie cluster
  "QrLBD26","QrLBD8","QrLBD39","QrLBD10","QrLBD20","QrLBD6",
  # Class Ic cluster
  "QrLBD34","QrLBD16","QrLBD15","QrLBD17","QrLBD13","QrLBD14","QrLBD18","QrLBD12",
  # Class Ia (small clade) + Class Ic continuation
  "QrLBD27","QrLBD4","QrLBD33","QrLBD5",
  # Class Id cluster
  "QrLBD25","QrLBD29","QrLBD28","QrLBD30","QrLBD36","QrLBD37","QrLBD35",
  # Class Ib cluster
  "QrLBD31","QrLBD40","QrLBD9","QrLBD1","QrLBD21","QrLBD23","QrLBD22",
  # Class Ia (second clade)
  "QrLBD19","QrLBD11","QrLBD32","QrLBD3",
  # Class IIb cluster
  "QrLBD42","QrLBD41","QrLBD24",
  # Class IIa cluster
  "QrLBD2","QrLBD38","QrLBD7"
)
stopifnot(length(phylo_qr_new) == 42)

present <- phylo_qr_new[phylo_qr_new %in% rownames(qr_z)]
qr_z_ordered <- qr_z[present, , drop = FALSE]

write.csv(qr_z_ordered, OUT_MAT, row.names = TRUE)
cat("  wrote z-score matrix:", OUT_MAT,
    "(", nrow(qr_z_ordered), "x", ncol(qr_z_ordered), ")\n\n")

# ---------------------------------------------------------------------------
# 9. Heatmap
# ---------------------------------------------------------------------------
my_colors <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100)
my_breaks <- seq(-2, 2, length.out = 101)

# pheatmap height should scale with number of rows
plot_h <- max(5, 0.18 * nrow(qr_z_ordered) + 2)

pheatmap(qr_z_ordered,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         scale = "none",
         main = sprintf("QrLBD Expression Across Tissues (z-score, %d genes after prefilter)",
                        nrow(qr_z_ordered)),
         fontsize_row = 7,
         fontsize_col = 10,
         cellwidth = 45,
         cellheight = 11,
         color = my_colors,
         breaks = my_breaks,
         border_color = "grey80",
         angle_col = "45",
         filename = OUT_PNG,
         width = 7, height = plot_h)

pheatmap(qr_z_ordered,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         scale = "none",
         main = sprintf("QrLBD Expression Across Tissues (z-score, %d genes after prefilter)",
                        nrow(qr_z_ordered)),
         fontsize_row = 7,
         fontsize_col = 10,
         cellwidth = 45,
         cellheight = 11,
         color = my_colors,
         breaks = my_breaks,
         border_color = "grey80",
         angle_col = "45",
         filename = OUT_PDF,
         width = 7, height = plot_h)

cat("  wrote heatmap:", OUT_PNG, "\n")
cat("  wrote heatmap:", OUT_PDF, "\n")
cat("==========================================================\n")
cat("  Done.\n")
cat("==========================================================\n")
