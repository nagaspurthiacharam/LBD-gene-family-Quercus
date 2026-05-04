#!/bin/bash
#============================================================
#  Genome-Wide Collinearity (Synteny) Analysis
#  LBD Gene Family — Quercus rubra Comparative Genomics
#
#  Species comparisons:
#    1. Q. rubra  vs  Q. suber         (congeneric)
#    2. Q. rubra  vs  P. trichocarpa   (woody eudicot)
#    3. Q. rubra  vs  A. thaliana      (herbaceous eudicot)
#    4. Q. rubra  vs  O. sativa        (monocot outgroup)
#
#  Tools: NCBI BLAST+ v2.12, MCScanX, Python 3 (matplotlib)
#  References:
#    Camacho et al. (2009) BMC Bioinformatics 10:421
#    Wang et al. (2012) Nucleic Acids Research 40(e49)
#============================================================

#------------------------------------------------------------
# STEP 1: Download and prepare genome resources
#------------------------------------------------------------

# Q. rubra v2.1 — Phytozome (12 linkage groups)
# Q. suber — NCBI RefSeq GCF_002906115.1 (scaffold-level)
# P. trichocarpa v4.1 — NCBI RefSeq GCF_000002775.5 (19 chromosomes)
# A. thaliana TAIR10.1 — NCBI GenBank GCA_000001735.2 (5 chromosomes)
# O. sativa indica 93-11 — NCBI GenBank GCA_000004655.2 (12 chromosomes)

# For each species, two files are needed:
#   protein.faa   — whole-proteome protein sequences
#   genomic.gff   — gene annotation in GFF3 format

#------------------------------------------------------------
# STEP 2: Fix GFF3 annotations for MCScanX compatibility
#------------------------------------------------------------

# NCBI GFF3 files use internal gene IDs that don't match
# protein FASTA headers. This script rebuilds the GFF3 so
# that gene entries carry the protein accession (e.g., XP_, EAY_)
# as their ID, enabling MCScanX to link BLAST hits to gene positions.

python3 fix_ncbi_gff3.py    # For RefSeq genomes (Q. suber, P. trichocarpa, A. thaliana)
python3 fix_rice_gff.py     # For GenBank rice genome (O. sativa indica)

# Phytozome GFF3 (Q. rubra) already uses matching IDs — no fix needed.

# Clean protein FASTA headers (remove descriptions, keep accession only)
# Output: protein_clean.faa for each species

#------------------------------------------------------------
# STEP 3: Convert GFF3 to MCScanX 4-column format
#------------------------------------------------------------

# MCScanX requires a simplified GFF with four tab-separated columns:
#   chromosome    geneID    start    end
#
# IMPORTANT: For Q. rubra (Phytozome), extract from "mRNA" features
# to obtain transcript-level IDs (e.g., Qurub.01G008700.1) that match
# the protein FASTA. Using "gene" features gives IDs without the .1
# suffix, causing zero inter-species matches.
#
# For NCBI-derived species, extract from "gene" features since the
# fixed GFF3 already maps protein IDs at the gene level.

# Q. rubra — use mRNA (critical: gives .1 suffix IDs)
awk -F'\t' '$3=="mRNA"{ split($9,a,";"); split(a[1],b,"="); gsub(/\./,"_",$1); print $1"\t"b[2]"\t"$4"\t"$5 }' \
    Qr_fixed.gff3 > Qr.gff

# Other species — use gene
awk -F'\t' '$3=="gene"{ split($9,a,";"); split(a[1],b,"="); gsub(/\./,"_",$1); print $1"\t"b[2]"\t"$4"\t"$5 }' \
    Sp2_fixed.gff3 > Sp2.gff

# Merge into combined GFF for each comparison
cat Qr.gff Sp2.gff > QrSp2.gff

#------------------------------------------------------------
# STEP 4: BLASTP all-vs-all
#------------------------------------------------------------

# Combine protein sequences from both species
cat Qr_protein.fa Sp2_protein.faa > QrSp2.fa

# Build BLAST database
makeblastdb -in QrSp2.fa -dbtype prot -out QrSp2_db

# Run BLASTP
#   -evalue 1e-10    : significance threshold
#   -outfmt 6        : tabular output (required by MCScanX)
#   -num_alignments 5: top 5 hits per query
#   -num_threads 4   : parallel threads
blastp -query QrSp2.fa \
       -db QrSp2_db \
       -out QrSp2.blast \
       -evalue 1e-10 \
       -num_threads 4 \
       -outfmt 6 \
       -num_alignments 5

#------------------------------------------------------------
# STEP 5: Run MCScanX
#------------------------------------------------------------

# MCScanX requires the .gff and .blast files to share the same
# prefix and be in the same directory.
#
# Parameters (defaults):
#   Match score: 50
#   Gap penalty: -1
#   E-value: 1e-05
#   Match size: 5 (minimum genes per collinear block)
#   Max gaps: 25

MCScanX QrSp2

# Output files:
#   QrSp2.collinearity  — collinear blocks with gene pairs
#   QrSp2.tandem         — tandem duplicate gene pairs
#   QrSp2.html/          — HTML visualization files

#------------------------------------------------------------
# STEP 6: Generate control (.ctl) files
#------------------------------------------------------------

# The .ctl file specifies visualization parameters:
#   Line 1: plot width (pixels)
#   Line 2: plot height (pixels)
#   Line 3: species 1 chromosomes (comma-separated)
#   Line 4: species 2 chromosomes (comma-separated)

python3 -c "
import re
chroms = set()
with open('QrSp2.gff') as f:
    for l in f:
        p = l.split('\t')
        chroms.add(p[0]) if len(p) >= 4 else None
sp1 = sorted([c for c in chroms if c.startswith('LG')],
             key=lambda n: [int(x) if x.isdigit() else x for x in re.split(r'(\d+)', n)])
sp2 = sorted([c for c in chroms if not c.startswith('LG')],
             key=lambda n: [int(x) if x.isdigit() else x for x in re.split(r'(\d+)', n)])
with open('QrSp2.ctl', 'w') as o:
    o.write('800\n800\n' + ','.join(sp1) + '\n' + ','.join(sp2) + '\n')
print(f'{len(sp1)} sp1 chromosomes + {len(sp2)} sp2 chromosomes')
"

#------------------------------------------------------------
# STEP 7: Generate dual synteny plots
#------------------------------------------------------------

# Python/matplotlib script that:
#   - Parses .collinearity for inter-species blocks
#   - Removes scaffolds with zero collinear blocks
#   - Draws bezier curves connecting syntenic gene pairs
#   - Highlights LBD gene pairs in red
#   - Outputs PNG (300 dpi) and PDF

python3 plot_synteny.py

#------------------------------------------------------------
# SUMMARY OF RESULTS
#------------------------------------------------------------

# Comparison         | Blocks | Gene Pairs | LBD Pairs
# -------------------+--------+------------+----------
# Qr vs Qs (oak)     | 1,057  |   11,012   |    14
# Qr vs Pt (poplar)  |   850  |   19,413   |    41
# Qr vs At (arabid.) |   711  |    8,916   |    15
# Qr vs Os (rice)    |   319  |    2,833   |     9
#
# The progressive decline in collinear blocks from congeneric
# (Q. suber) to monocot outgroup (O. sativa) is consistent with
# increasing evolutionary distance and genome rearrangement.
