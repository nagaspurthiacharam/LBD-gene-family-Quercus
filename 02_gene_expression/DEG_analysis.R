
library(DESeq2)
library(dplyr)

############################################################
# STEP 1: Parse GFF to map RGQ29 gene IDs ↔ Qurub transcript IDs
############################################################

# Read GFF, only mRNA lines
gff <- read.table("genomic.gff", sep="\t", header=FALSE, 
                  comment.char="#", quote="", stringsAsFactors=FALSE,
                  fill=TRUE)

# Keep only mRNA rows
mrna <- gff[gff$V3 == "mRNA", ]

# Extract IDs using regex
mrna$rna_id  <- sub(".*ID=(rna-Qurub\\.[^;]+).*", "\\1", mrna$V9)
mrna$gene_id <- sub(".*Parent=(gene-RGQ29_[0-9]+).*", "\\1", mrna$V9)

# Mapping table
mapping <- mrna[, c("rna_id", "gene_id")]
head(mapping)

############################################################
# STEP 2: Read your 47 QrLBD rna IDs
############################################################
qrlbd_ids <- readLines("rna_ids.txt")
qrlbd_ids <- trimws(qrlbd_ids)
qrlbd_ids <- qrlbd_ids[qrlbd_ids != ""]  # remove blank lines
length(qrlbd_ids)  # should be ~47

############################################################
# STEP 3: Match QrLBD rna IDs to RGQ29 gene IDs
############################################################
qrlbd_mapping <- mapping[mapping$rna_id %in% qrlbd_ids, ]
cat("Matched:", nrow(qrlbd_mapping), "out of", length(qrlbd_ids), "\n")

# Check which ones didn't match
missing <- setdiff(qrlbd_ids, qrlbd_mapping$rna_id)
if(length(missing) > 0) {
  cat("Missing IDs:\n"); print(missing)
}

# Save the mapping
write.csv(qrlbd_mapping, "qrlbd_to_rgq29_mapping.csv", row.names=FALSE)

############################################################
# STEP 4: Load count matrix and filter to QrLBD genes
############################################################
countdata <- read.csv("gene_count_matrix.csv", row.names=1, check.names=FALSE)
dim(countdata)
head(countdata)

# Get the unique RGQ29 gene IDs for your QrLBDs
lbd_gene_ids <- unique(qrlbd_mapping$gene_id)

# Filter count matrix
lbd_counts <- countdata[rownames(countdata) %in% lbd_gene_ids, ]
cat("LBD genes found in count matrix:", nrow(lbd_counts), "\n")

# Save filtered counts
write.csv(lbd_counts, "QrLBD_counts.csv")

head(lbd_counts)



############################################################
# STEP 4 (revised): Keep ALL 47 rows with QrLBD names
############################################################

# Read your full mapping with QrLBD names (the Excel file you showed)
# Save it as a tab-separated txt file first, or read directly from CSV
qrlbd_full <- read.csv("qrlbd_to_rgq29_mapping.csv", stringsAsFactors=FALSE)
# Should have columns: rna_id, gene_id, QrLBDs
head(qrlbd_full)
nrow(qrlbd_full)  # 47

# Load count matrix
countdata <- read.csv("gene_count_matrix.csv", row.names=1, check.names=FALSE)

# Merge: for each of the 47 QrLBDs, pull the counts using gene_id
# This will duplicate counts for isoforms that share the same gene_id
lbd_counts_all <- merge(qrlbd_full, 
                        data.frame(gene_id=rownames(countdata), countdata, 
                                   check.names=FALSE),
                        by="gene_id", all.x=TRUE)

# Reorder columns: QrLBD name first, then rna_id, gene_id, then counts
lbd_counts_all <- lbd_counts_all[, c("QrLBDs", "rna_id", "gene_id", 
                                     colnames(countdata))]

# Sort by QrLBD number (optional, nicer for viewing)
lbd_counts_all <- lbd_counts_all[order(as.numeric(gsub("[^0-9.]", "", 
                                                       lbd_counts_all$QrLBDs))), ]

nrow(lbd_counts_all)  # should be 47
head(lbd_counts_all)

# Save
write.csv(lbd_counts_all, "QrLBD_counts_all47.csv", row.names=FALSE)