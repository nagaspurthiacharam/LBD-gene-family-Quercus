###############################################################################
#  Heatmap_zscore_QsLBD.R
#  --------------------------------------------------------------------------
#  Builds the publication-quality z-score expression heatmap for the 37
#  primary QsLBD genes across the five Q. suber multi-tissue libraries
#  from PRJNA392919 (phellem, xylem, leaf, inner bark, pollen).
#
#  Data this script needs (all sit next to it in Counts_H_S/):
#       gene_count_matrix.csv          - StringTie / prepDE.py3 output
#                                        (43k genes x 5 SRA libraries)
#       my_LBD_raw_counts_MT_named.csv - existing 37-row mapping with
#                                        columns: ID, Name, <5 tissue cols>
#                                        (the 'Name' column gives QsLBD1..37)
#
#  Pipeline (matches the DESeq2 prefilter convention from 4_DESeq.r:
#  dds <- dds[ rowSums(counts(dds)) > 20, ])
#  --------
#     1. Load FULL gene_count_matrix.csv (~43k genes x 5 libraries).
#     2. Pre-filter on the full matrix: rowSums(counts) > 20.
#     3. Library-size normalisation - per-library size factor =
#        column-total / median(column-total) (median-of-ratios style;
#        equivalent in spirit to DESeq2 estimateSizeFactors() with no
#        Bioconductor dependency).
#     4. Subset to the 37 primary QsLBDs via the LOC IDs in
#        my_LBD_raw_counts_MT_named.csv. Each column is already a unique
#        tissue (no biological replicates), so no averaging is needed.
#     5. Per-gene z-score across the 5 tissues, capped at +-2 to keep the
#        diverging colour ramp readable.
#     6. Re-order rows by the mid-point-rooted QsLBD phylogeny (read from
#        QsLBD_phylogeny_order.txt under oak/Motif/ if available, otherwise
#        falls back to the embedded vector below).
#     7. Render with pheatmap, RdYlBu reversed palette, breaks [-2, 2],
#        cluster_rows = FALSE / cluster_cols = FALSE.
#
#  Run from the Counts_H_S/ folder:
#       Rscript Heatmap_zscore_QsLBD.R
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
COUNTS    <- "../data/QsLBD_counts_37_primary.csv"
NAMED     <- "my_LBD_raw_counts_MT_named.csv"
PHYLO     <- "../../../Motif/QsLBD_phylogeny_order.txt"   # optional override
OUT_PNG   <- "QsLBD_zscore_heatmap.png"
OUT_PDF   <- "QsLBD_zscore_heatmap.pdf"
OUT_MAT   <- "QsLBD_zscore_matrix.csv"

PREFILTER_MIN_TOTAL <- 20    # matches DESeq2 prefilter convention

cat("==========================================================\n")
cat("  QsLBD z-score expression heatmap\n")
cat("  prefilter:  rowSums(counts) >", PREFILTER_MIN_TOTAL,
    "  (DESeq2-style)\n")
cat("==========================================================\n")
cat("  working dir : ", normalizePath(WORK_DIR), "\n")
cat("  counts in   : ", COUNTS, "\n")
cat("  LBD mapping : ", NAMED,  "\n")
cat("  output PNG  : ", OUT_PNG, "\n\n")

stopifnot(file.exists(COUNTS))
stopifnot(file.exists(NAMED))

# ---------------------------------------------------------------------------
# 1. Load named LBD mapping (37 rows)
# ---------------------------------------------------------------------------
named <- read.csv(NAMED, stringsAsFactors = FALSE, check.names = FALSE)
cat("  LBD-named CSV: ", nrow(named), "rows,",
    ncol(named), "columns\n")

# expected columns: ID, Name, then 5 tissue columns
stopifnot(all(c("ID", "Name") %in% colnames(named)))
tissue_cols <- setdiff(colnames(named), c("ID", "Name"))
cat("  tissue columns (named CSV order):",
    paste(tissue_cols, collapse = ", "), "\n")

# ---------------------------------------------------------------------------
# 2. Load FULL count matrix
# ---------------------------------------------------------------------------
counts_full <- as.matrix(read.csv(COUNTS, row.names = 1, check.names = FALSE))
mode(counts_full) <- "numeric"
cat("\n  full count matrix:", nrow(counts_full), "genes x ",
    ncol(counts_full), "samples\n")
cat("  raw sample names :", paste(colnames(counts_full), collapse = ", "), "\n")

# Sample -> tissue map (PRJNA392919 / Post_preclean_data_logistics_QS.txt)
# Clean tissue names (no underscores) for readable heatmap column labels.
sample_to_tissue <- c(
  SRR5986737 = "Phellem",
  SRR5986738 = "Xylem",
  SRR5986739 = "Leaf",
  SRR5986740 = "Inner bark",
  SRR5986741 = "Pollen"
)

# Confirm every column has a tissue mapping; warn on any extra
miss_meta <- setdiff(colnames(counts_full), names(sample_to_tissue))
if (length(miss_meta)) {
  warning("Sample column(s) without metadata (will be dropped): ",
          paste(miss_meta, collapse = ", "))
  counts_full <- counts_full[, !colnames(counts_full) %in% miss_meta,
                             drop = FALSE]
}

# Rename the columns to tissue names so the heatmap labels are biological
colnames(counts_full) <- sample_to_tissue[colnames(counts_full)]
cat("  tissue columns now:", paste(colnames(counts_full), collapse = ", "), "\n")

# ---------------------------------------------------------------------------
# 3. Pre-filter on the FULL matrix (DESeq2-style rowSums > 20)
# ---------------------------------------------------------------------------
keep <- rowSums(counts_full) > PREFILTER_MIN_TOTAL
cat("\n  prefilter:  rowSums >", PREFILTER_MIN_TOTAL, "  ->  kept",
    sum(keep), "/", length(keep), " genes\n")
counts_filt <- counts_full[keep, , drop = FALSE]

# ---------------------------------------------------------------------------
# 4. Library-size normalisation on the filtered matrix
# ---------------------------------------------------------------------------
col_tot   <- colSums(counts_filt)
size_fact <- col_tot / median(col_tot)
norm_full <- sweep(counts_filt, 2, size_fact, FUN = "/")
cat("\n  size factors (per library):\n")
print(round(size_fact, 3))

# ---------------------------------------------------------------------------
# 5. Subset to the 37 primary QsLBDs (those that survive prefilter)
# ---------------------------------------------------------------------------
qs_in_full <- intersect(named$ID, rownames(counts_full))
qs_in_filt <- intersect(named$ID, rownames(counts_filt))
qs_dropped <- setdiff(qs_in_full, qs_in_filt)

cat("\n  37 primary QsLBDs:\n")
cat("     present in count matrix     :", length(qs_in_full), "/ 37\n")
cat("     surviving rowSums >",  PREFILTER_MIN_TOTAL, "       :",
    length(qs_in_filt), "/ 37\n")
cat("     dropped at prefilter step   :", length(qs_dropped), "\n")
if (length(qs_dropped)) {
  drop_names <- named$Name[match(qs_dropped, named$ID)]
  drop_sums  <- rowSums(counts_full[qs_dropped, , drop = FALSE])
  cat("     (genes silent or near-silent across the 5 tissues)\n")
  for (i in seq_along(qs_dropped)) {
    cat(sprintf("        %-9s rowSum = %d\n",
                drop_names[i], as.integer(drop_sums[i])))
  }
}

lbd_norm <- norm_full[qs_in_filt, , drop = FALSE]
loc_to_name <- setNames(named$Name, named$ID)
rownames(lbd_norm) <- loc_to_name[rownames(lbd_norm)]

# ---------------------------------------------------------------------------
# 6. Per-gene z-score, capped at +-2  (no replicate averaging needed)
# ---------------------------------------------------------------------------
z_score <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}
qs_z <- t(apply(lbd_norm, 1, z_score))
qs_z <- pmax(pmin(qs_z, 2), -2)
colnames(qs_z) <- colnames(lbd_norm)

cat("\n  z-scored matrix:", nrow(qs_z), "genes x ",
    ncol(qs_z), "tissues\n")

# ---------------------------------------------------------------------------
# 7. Phylogeny order (mid-point-rooted QsLBD tree leaf order, top-to-bottom)
# ---------------------------------------------------------------------------
phylo_qs <- if (file.exists(PHYLO)) {
  cat("\n  reading phylogeny order from:", PHYLO, "\n")
  readLines(PHYLO, warn = FALSE)
} else {
  c(
    "QsLBD27","QsLBD9","QsLBD8","QsLBD13","QsLBD2","QsLBD12","QsLBD24",
    "QsLBD29","QsLBD23","QsLBD20","QsLBD15","QsLBD1","QsLBD22","QsLBD21",
    "QsLBD25","QsLBD3","QsLBD34","QsLBD10","QsLBD30","QsLBD37","QsLBD35",
    "QsLBD11","QsLBD4","QsLBD36","QsLBD31","QsLBD33","QsLBD32","QsLBD28",
    "QsLBD5","QsLBD7","QsLBD6","QsLBD19","QsLBD18","QsLBD17","QsLBD16",
    "QsLBD26","QsLBD14"
  )
}
phylo_qs <- trimws(phylo_qs)
phylo_qs <- phylo_qs[nchar(phylo_qs) > 0]
stopifnot(length(phylo_qs) == 37)

present <- phylo_qs[phylo_qs %in% rownames(qs_z)]
qs_z_ordered <- qs_z[present, , drop = FALSE]

write.csv(qs_z_ordered, OUT_MAT, row.names = TRUE)
cat("  wrote z-score matrix:", OUT_MAT,
    "(", nrow(qs_z_ordered), "x", ncol(qs_z_ordered), ")\n\n")

# ---------------------------------------------------------------------------
# 8. Heatmap
# ---------------------------------------------------------------------------
my_colors <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100)
my_breaks <- seq(-2, 2, length.out = 101)

plot_h <- max(5, 0.18 * nrow(qs_z_ordered) + 2)

pheatmap(qs_z_ordered,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         scale = "none",
         main = sprintf("QsLBD Multi-tissue Expression (z-score, %d genes after prefilter)",
                        nrow(qs_z_ordered)),
         fontsize_row = 7,
         fontsize_col = 10,
         cellwidth = 50,
         cellheight = 11,
         color = my_colors,
         breaks = my_breaks,
         border_color = "grey80",
         angle_col = "45",
         filename = OUT_PNG,
         width = 6, height = plot_h)

pheatmap(qs_z_ordered,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         scale = "none",
         main = sprintf("QsLBD Multi-tissue Expression (z-score, %d genes after prefilter)",
                        nrow(qs_z_ordered)),
         fontsize_row = 7,
         fontsize_col = 10,
         cellwidth = 50,
         cellheight = 11,
         color = my_colors,
         breaks = my_breaks,
         border_color = "grey80",
         angle_col = "45",
         filename = OUT_PDF,
         width = 6, height = plot_h)

cat("  wrote heatmap:", OUT_PNG, "\n")
cat("  wrote heatmap:", OUT_PDF, "\n")
cat("==========================================================\n")
cat("  Done.\n")
cat("==========================================================\n")
