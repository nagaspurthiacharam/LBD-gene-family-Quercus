# Data Directory

This folder contains all curated input files used by the analysis pipeline. Every file here is read by one or more scripts in `02_gene_expression/`, `03_collinearity/`, or `04_cis_promoter/`. No file is generated on the fly during a normal pipeline run.

| # | File | Format | Records | Used by |
|---|---|---|---|---|
| 1 | `QrLBD_gene_catalogue.csv` | CSV (17 cols) | 47 (42 primary + 5 sub-iso) | extract_QrLBD_counts.R, all heatmaps |
| 2 | `LBD_gene_IDs.txt` | text, 1 ID/line | 79 | plot_synteny_all.py |
| 3 | `QrLBD_phylogeny_order.txt` | text, 1 name/line | 42 | Heatmap_zscore_QrLBD.R, Cis_tile_heatmap.R |
| 4 | `QsLBD_phylogeny_order.txt` | text, 1 name/line | 37 | Heatmap_zscore_QsLBD.R, Cis_tile_heatmap.R |
| 5 | `QrLBD_proteins.fasta` | FASTA | 42 sequences | phylogeny, motif scanning, structural modelling |
| 6 | `QsLBD_proteins.fasta` | FASTA | 37 sequences | phylogeny, motif scanning, structural modelling |
| 7 | `QrLBD_promoters_2kb.fasta` | FASTA | 42 sequences | PlantCARE cis-element prediction |
| 8 | `QsLBD_promoters_2kb.fasta` | FASTA | 37 sequences | PlantCARE cis-element prediction |
| 9 | `QrLBD_counts_42_primary.csv` | CSV | 42 rows x 11 libs | Heatmap_zscore_QrLBD.R |
|10 | `QsLBD_counts_37_primary.csv` | CSV | 37 rows x 5 libs | Heatmap_zscore_QsLBD.R |

---

## 1. `QrLBD_gene_catalogue.csv` — single source of truth

The canonical catalogue of every QrLBD locus (42 primary + 5 sub-isoform). Used by the count-extraction and heatmap scripts to (i) discover which RGQ29 gene IDs to pull from the StringTie count matrix, (ii) keep only primary loci, and (iii) attach human-readable names. **Sub-isoforms have a non-empty `Sub-isoform of` field; primaries have it empty.**

| Column | Type | Description |
|---|---|---|
| `Name` | text | Canonical paper name (`QrLBD1` ... `QrLBD42`, plus `QrLBD1.1`, `QrLBD19.1`, `QrLBD19.2`, `QrLBD22.1`, `QrLBD40.1`) |
| `rna_id` | text | NCBI RefSeq mRNA accession of the form `rna-Qurub.NNGNNNNNN.N.v2.1` |
| `gene_id` | text | NCBI RefSeq locus tag of the form `gene-RGQ29_NNNNNN` (links rows to the StringTie gene_count_matrix.csv) |
| `GeneId (Phytozome)` | text | Phytozome v13 transcript ID, e.g. `Qurub.01G008700.1` |
| `Locus` | text | Phytozome v13 gene locus (transcript-suffix stripped) |
| `Class` | text | LBD subclass: `Class Ia`, `Ib`, `Ic`, `Id`, `Ie`, `IIa`, `IIb` |
| `Chromosome` | text | Linkage group (LG1-LG12) |
| `Start`, `End` | integer | Genomic coordinates on the linkage group (bp) |
| `Length (aa)` | integer | Protein length |
| `MW (kDa)` | float | Theoretical molecular weight (ProtParam) |
| `pI` | float | Theoretical isoelectric point |
| `Instability` | float | Instability index (Guruprasad et al., 1990) |
| `Aliphatic` | float | Aliphatic index |
| `GRAVY` | float | Grand average of hydropathicity |
| `Localisation` | text | BUSCA-predicted subcellular localisation |
| `Sub-isoform of` | text | empty for primaries; QrLBD parent name for sub-isoforms |

## 2. `LBD_gene_IDs.txt`

Plain text, one identifier per line, used by `plot_synteny_all.py` to highlight LBD-anchored ribbons (red) on the synteny plots. Contains 79 entries:

- 42 *Q. rubra* loci as Phytozome IDs (`Qurub.01G008700.1`, etc.)
- 37 *Q. suber* loci as NCBI RefSeq XP_ accessions (`XP_065628701.1`, etc.)

## 3. `QrLBD_phylogeny_order.txt` (42 lines)

Mid-point-rooted leaf order of the QrLBD subtree (top-to-bottom). Used to enforce a consistent biological row-ordering across the expression heatmap, the gene-structure / motif diagram, and the cis-element heatmap. One name per line, no header.

## 4. `QsLBD_phylogeny_order.txt` (37 lines)

Same structure as item 3 but for the QsLBD subtree.

## 5. `QrLBD_proteins.fasta` (42 sequences)

Protein sequences for the 42 primary QrLBD loci, headers in the form `>QrLBD1`, `>QrLBD2`, ... `>QrLBD42`. Used as input for the multi-species phylogenetic alignment, MEME motif discovery, SWISS-MODEL homology modelling, and the STRING orthologue mapping.

## 6. `QsLBD_proteins.fasta` (37 sequences)

Same as item 5 but for the 37 QsLBD primary loci, headers `>QsLBD1` ... `>QsLBD37`.

## 7. `QrLBD_promoters_2kb.fasta`

The 2-kb upstream sequence immediately preceding the predicted translation start codon (ATG) of each QrLBD primary locus, extracted via the TBtools-II Promoter Sequence Extractor module. Submitted as-is to the PlantCARE database for cis-acting regulatory element prediction.

## 8. `QsLBD_promoters_2kb.fasta`

Same as item 7 but for the QsLBD primary loci.

## 9. `QrLBD_counts_42_primary.csv` (42 rows x 11 libraries)

LBD-family-only subset of the StringTie / `prepDE.py3` count matrix for *Quercus rubra*. Each row is a primary QrLBD locus; columns are `QrLBD`, `rna_id`, `gene_id`, then 11 SRA-accession sample columns from BioProject **PRJNA273270** (six tissues x two parental genotypes; one tissue is single-genotype). This is the direct input to `Heatmap_zscore_QrLBD.R`. The whole-genome 33,247-gene matrix from which this subset was derived is regenerable by re-running the upstream HISAT2 / StringTie / prepDE.py3 pipeline on PRJNA273270 (scripts in `02_gene_expression/PE_paired_end/`) — it is not shipped here for size and relevance reasons.

## 10. `QsLBD_counts_37_primary.csv` (37 rows x 5 libraries)

Same format as item 9 but for *Quercus suber*. Columns are `QsLBD`, `gene_id`, then five SRA-accession sample columns from BioProject **PRJNA392919** (one library per tissue: phellem, xylem, leaf, inner bark, pollen). Direct input to `Heatmap_zscore_QsLBD.R`.

---

## How the data files connect

```
QrLBD_gene_catalogue.csv
        |
        +--- gene_id ----> StringTie gene_count_matrix.csv  (extract_QrLBD_counts.R)
        |
        +--- Name -------> QrLBD_phylogeny_order.txt        (Heatmap_zscore_QrLBD.R,
        |                                                    Cis_tile_heatmap.R)
        |
        +--- Sub-isoform of (empty)  -> primary set used everywhere downstream

QrLBD_proteins.fasta  --> MEGA 7 (phylogeny, Figure 1)
                      --> MEME Suite (motifs, Figure 3)
                      --> SWISS-MODEL (structures, Figure 9)
                      --> STRING (PPI networks, Figure 8)

QrLBD_promoters_2kb.fasta --> PlantCARE
                            --> Cis_curated_QrLBD42.xlsx
                            --> Cis_tile_heatmap.R (Figure 6, Supp Figs S2/S3)
```

## Data provenance and authority

- **Q. rubra** sequences and coordinates are derived from the Phytozome v13 *Quercus rubra* v2.1 release (https://phytozome.jgi.doe.gov; Kapoor et al., 2023).
- **Q. suber** sequences and coordinates are derived from NCBI RefSeq assembly **GCF_002906115.3** (BioProject PRJNA392919; Ramos et al., 2018).
- **Subclass assignment** (Class Ia-Ie, IIa, IIb) follows clustering with the previously characterised *A. thaliana* LBD set as described in the manuscript Section 2.4.
- **Sub-isoform tagging** of QrLBD1.1 / QrLBD19.1 / QrLBD19.2 / QrLBD22.1 / QrLBD40.1 was made at the ProtParam quality-control step described in manuscript Section 2.2; these rows share genomic coordinates with their principal-isoform counterparts.

If you need to reproduce any figure from scratch you should not need any input outside this folder plus the public RNA-seq reads (PRJNA273270 and PRJNA392919).
