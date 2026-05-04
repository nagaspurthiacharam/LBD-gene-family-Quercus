###############################################################################
#  extract_QrLBD_counts.R
#  --------------------------------------------------------------------------
#  Pull tissue-resolved RNA-seq counts for the 42 QrLBD primary loci ONLY
#  from the StringTie / prepDE.py3 gene_count_matrix.csv. Sub-isoforms are
#  ignored: they share gene_id with their parent and StringTie counts at
#  gene level, so including them would just duplicate the parent row.
#
#  Why this script: the previous workflow tried to parse genomic.gff with
#  regex against an rna_ids.txt that contained typos ("Qururub" instead of
#  "Qurub"), so every match returned NA and the merge produced an empty
#  table. The mapping work has already been done once and saved in
#  qrlbd_to_rgq29_mapping_NEW.csv with five canonical columns:
#       rna_id , gene_id , QrLBDs_old , QrLBDs_new , status
#  This script reads that mapping, drops sub-isoform rows, joins to the
#  count matrix on gene_id, and writes ONE output:
#       QrLBD_counts_42_primary.csv     (42 primary loci only)
#
#  Run from the Quercus_rubra/ folder:
#       Rscript extract_QrLBD_counts.R
#  or in RStudio: setwd() to this folder and source().
###############################################################################

# ----------------------------------------------------------------------------
# 0. Setup
# ----------------------------------------------------------------------------
WORK_DIR <- "."         # script reads/writes in current working directory
MAPPING  <- "qrlbd_to_rgq29_mapping_NEW.csv"
COUNTS   <- "gene_count_matrix.csv"
OUT_42   <- "QrLBD_counts_42_primary.csv"

cat("==========================================================\n")
cat("  QrLBD count-matrix extractor (42 primary loci ONLY)\n")
cat("==========================================================\n")
cat("  working dir : ", normalizePath(WORK_DIR), "\n")
cat("  mapping     : ", MAPPING, "\n")
cat("  counts      : ", COUNTS,  "\n\n")

stopifnot(file.exists(MAPPING))
stopifnot(file.exists(COUNTS))

# ----------------------------------------------------------------------------
# 1. Load mapping and KEEP ONLY THE 42 PRIMARY ROWS
# ----------------------------------------------------------------------------
mapping_all <- read.csv(MAPPING, stringsAsFactors = FALSE)
need_cols <- c("rna_id", "gene_id", "QrLBDs_new", "status")
missing_cols <- setdiff(need_cols, colnames(mapping_all))
if (length(missing_cols)) {
  stop("Mapping file is missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

mapping <- mapping_all[mapping_all$status == "primary", ]
cat("  mapping rows kept:", nrow(mapping), "(primary only)\n")

# ----------------------------------------------------------------------------
# 2. Load count matrix (rownames = gene-RGQ29_NNNNNN)
# ----------------------------------------------------------------------------
counts <- read.csv(COUNTS, row.names = 1, check.names = FALSE)
cat("  count matrix:", nrow(counts), "genes x ",
    ncol(counts), "samples\n")
cat("  samples:", paste(colnames(counts), collapse = ", "), "\n\n")

# ----------------------------------------------------------------------------
# 3. Sanity: how many of our 42 gene_ids actually appear in the count matrix?
# ----------------------------------------------------------------------------
target_genes <- unique(mapping$gene_id)
hit          <- target_genes[target_genes %in% rownames(counts)]
miss         <- setdiff(target_genes, rownames(counts))
cat("  gene_ids in mapping that are in count matrix:",
    length(hit), "/", length(target_genes), "\n")
if (length(miss) > 0) {
  cat("  *** missing gene_ids:\n")
  for (m in miss) cat("        ", m, "\n")
}

# ----------------------------------------------------------------------------
# 4. Join: pull counts row by row, keep all 42 mapping rows
# ----------------------------------------------------------------------------
counts_df <- data.frame(
  gene_id = rownames(counts),
  counts,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

joined <- merge(mapping, counts_df, by = "gene_id", all.x = TRUE)

# Reorder columns: QrLBD name, rna_id, gene_id, then sample counts
sample_cols <- colnames(counts)
joined <- joined[, c("QrLBDs_new", "rna_id", "gene_id", sample_cols)]
colnames(joined)[1] <- "QrLBD"

# Sort by QrLBD number (1, 2, 3, ... 42)
qrlbd_num <- suppressWarnings(
  as.integer(sub("QrLBD", "", joined$QrLBD))
)
joined <- joined[order(qrlbd_num), ]

cat("\n  joined table:", nrow(joined), "rows\n")
cat("  preview (first 3 rows):\n")
print(head(joined[, 1:min(7, ncol(joined))], 3))

# ----------------------------------------------------------------------------
# 5. Write output (42-row CSV)
# ----------------------------------------------------------------------------
write.csv(joined, OUT_42, row.names = FALSE)

cat("\n  wrote:", OUT_42, "(", nrow(joined), "rows )\n")
cat("==========================================================\n")
cat("  Done. Use ", OUT_42, " for the heatmap analysis.\n", sep = "")
cat("==========================================================\n")