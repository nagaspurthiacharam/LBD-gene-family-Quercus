###############################################################################
#  extract_QrLBD_counts.R
#  --------------------------------------------------------------------------
#  Convenience wrapper.
#
#  The 42-row primary count CSV used by the heatmap pipeline is already
#  shipped in this repository as data/QrLBD_counts_42_primary.csv.
#  This script exists for transparency: if a user wants to regenerate it
#  from the full whole-genome StringTie matrix (which we do NOT ship for
#  size/relevance reasons), set FULL_MATRIX below to a local copy of the
#  StringTie / prepDE.py3 gene_count_matrix.csv (33,247 rows for the
#  PRJNA273270 alignment) and source this script.
#
#  Run from 02_gene_expression/:
#       Rscript extract_QrLBD_counts.R
#  or in RStudio:  setwd() to this folder, then source().
###############################################################################

WORK_DIR    <- "."
CATALOGUE   <- "../data/QrLBD_gene_catalogue.csv"
FULL_MATRIX <- "gene_count_matrix.csv"   # supply locally if regenerating
OUT_42      <- "../data/QrLBD_counts_42_primary.csv"

cat("==========================================================\n")
cat("  QrLBD count subset extractor\n")
cat("==========================================================\n")

stopifnot(file.exists(CATALOGUE))
if (!file.exists(FULL_MATRIX)) {
  cat("\n  Full whole-genome count matrix not found at:\n")
  cat("     ", FULL_MATRIX, "\n\n")
  cat("  This is expected for a fresh checkout. The 42-row LBD subset is\n")
  cat("  already shipped in the repository at:\n")
  cat("     ", OUT_42, "\n\n")
  cat("  To REGENERATE the subset from scratch, place a copy of the\n")
  cat("  StringTie / prepDE.py3 gene_count_matrix.csv (BioProject\n")
  cat("  PRJNA273270 alignment) next to this script and rerun.\n")
  quit(status = 0)
}

# Build mapping of gene_id -> (Name, rna_id) for the 42 primary loci
cat_all   <- read.csv(CATALOGUE, stringsAsFactors = FALSE, check.names = FALSE)
primaries <- cat_all[is.na(cat_all[["Sub-isoform of"]]) |
                     cat_all[["Sub-isoform of"]] == "", ]
cat("  primaries in catalogue:", nrow(primaries), "\n")

counts <- read.csv(FULL_MATRIX, row.names = 1, check.names = FALSE)
cat("  full count matrix    :", nrow(counts), "genes x ",
    ncol(counts), "samples\n")

counts_df <- data.frame(gene_id = rownames(counts), counts,
                        check.names = FALSE, stringsAsFactors = FALSE)
joined <- merge(primaries[, c("Name", "rna_id", "gene_id")],
                counts_df, by = "gene_id", all.x = TRUE)
sample_cols <- colnames(counts)
joined <- joined[, c("Name", "rna_id", "gene_id", sample_cols)]
colnames(joined)[1] <- "QrLBD"
joined <- joined[order(suppressWarnings(
  as.integer(sub("QrLBD", "", joined$QrLBD)))), ]

write.csv(joined, OUT_42, row.names = FALSE)
cat("\n  wrote (overwrote):", OUT_42, " (", nrow(joined), "rows )\n")
cat("==========================================================\n")
