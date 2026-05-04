###############################################################################
#  extract_QsLBD_counts.R
#  --------------------------------------------------------------------------
#  Convenience wrapper (Qs counterpart of extract_QrLBD_counts.R).
#
#  The 37-row primary count CSV used by the heatmap pipeline is already
#  shipped in this repository as data/QsLBD_counts_37_primary.csv.
#  This script exists for transparency: if a user wants to regenerate it
#  from the full whole-genome StringTie matrix (which we do NOT ship for
#  size/relevance reasons), set FULL_MATRIX below to a local copy of the
#  StringTie / prepDE.py3 gene_count_matrix.csv (43,163 rows for the
#  PRJNA392919 alignment) and source this script.
#
#  Run from 02_gene_expression/:
#       Rscript extract_QsLBD_counts.R
#  or in RStudio:  setwd() to this folder, then source().
###############################################################################

WORK_DIR    <- "."
SHIPPED_OUT <- "../data/QsLBD_counts_37_primary.csv"   # what we already ship
FULL_MATRIX <- "gene_count_matrix.csv"                 # supply locally if regenerating

# QsLBD name <-> NCBI LOC id mapping. The 37 primary QsLBDs in
# QsLBD1 .. QsLBD37 order, with their gene_count_matrix gene_id values.
qs_map <- read.table(text = "
QsLBD,gene_id
QsLBD1,gene-LOC112040411|LOC112040411
QsLBD2,gene-LOC112032524|LOC112032524
QsLBD3,gene-LOC111998864|LOC111998864
QsLBD4,gene-LOC112000770|LOC112000770
QsLBD5,gene-LOC112004483|LOC112004483
QsLBD6,gene-LOC112004479|LOC112004479
QsLBD7,gene-LOC112004478|LOC112004478
QsLBD8,gene-LOC112039894|LOC112039894
QsLBD9,gene-LOC112037478|LOC112037478
QsLBD10,gene-LOC111996809|LOC111996809
QsLBD11,gene-LOC112027411|LOC112027411
QsLBD12,gene-LOC112038859|LOC112038859
QsLBD13,gene-LOC112035794|LOC112035794
QsLBD14,gene-LOC111997202|LOC111997202
QsLBD15,gene-LOC112038573|LOC112038573
QsLBD16,gene-LOC111995028|LOC111995028
QsLBD17,gene-LOC112024999|LOC112024999
QsLBD18,gene-LOC112022402|LOC112022402
QsLBD19,gene-LOC112005907|LOC112005907
QsLBD20,gene-LOC112039236|LOC112039236
QsLBD21,gene-LOC112016164|LOC112016164
QsLBD22,gene-LOC112020652|LOC112020652
QsLBD23,gene-LOC112020807|LOC112020807
QsLBD24,gene-LOC112019729|LOC112019729
QsLBD25,gene-LOC111986710|LOC111986710
QsLBD26,gene-LOC112017929|LOC112017929
QsLBD27,gene-LOC112018569|LOC112018569
QsLBD28,gene-LOC112019139|LOC112019139
QsLBD29,gene-LOC112016015|LOC112016015
QsLBD30,gene-LOC112015218|LOC112015218
QsLBD31,gene-LOC112015558|LOC112015558
QsLBD32,gene-LOC112023207|LOC112023207
QsLBD33,gene-LOC112016061|LOC112016061
QsLBD34,gene-LOC112041157|LOC112041157
QsLBD35,gene-LOC112020410|LOC112020410
QsLBD36,gene-LOC112019486|LOC112019486
QsLBD37,gene-LOC112024301|LOC112024301
", header = TRUE, sep = ",", stringsAsFactors = FALSE)

cat("==========================================================\n")
cat("  QsLBD count subset extractor\n")
cat("==========================================================\n")

if (!file.exists(FULL_MATRIX)) {
  cat("\n  Full whole-genome count matrix not found at:\n")
  cat("     ", FULL_MATRIX, "\n\n")
  cat("  This is expected for a fresh checkout. The 37-row LBD subset is\n")
  cat("  already shipped in the repository at:\n")
  cat("     ", SHIPPED_OUT, "\n\n")
  cat("  To REGENERATE the subset from scratch, place a copy of the\n")
  cat("  StringTie / prepDE.py3 gene_count_matrix.csv (BioProject\n")
  cat("  PRJNA392919 alignment) next to this script and rerun.\n")
  quit(status = 0)
}

# Load full Qs whole-genome matrix
counts <- read.csv(FULL_MATRIX, row.names = 1, check.names = FALSE)
cat("  full Qs count matrix :", nrow(counts), "genes x ",
    ncol(counts), "samples\n")

# Build the 37-row subset
sample_cols <- colnames(counts)
out <- data.frame(QsLBD = qs_map$QsLBD,
                  gene_id = qs_map$gene_id,
                  stringsAsFactors = FALSE)
for (s in sample_cols) {
  out[[s]] <- counts[match(qs_map$gene_id, rownames(counts)), s]
}

write.csv(out, SHIPPED_OUT, row.names = FALSE)
cat("\n  wrote (overwrote):", SHIPPED_OUT,
    "  (", nrow(out), "rows x", 2 + length(sample_cols), "cols )\n")
cat("==========================================================\n")
