#!/bin/bash
#============================================================
#  LBD Gene Family Identification using HMMER
#
#  Identifies LBD (Lateral Organ Boundary Domain) gene family
#  members in plant proteomes using the LOB domain profile.
#
#  Tools: HMMER 3.3, NCBI BLAST+, NCBI CDD
#  References:
#    Eddy (2011) PLoS Comput Biol 7:e1002195
#    Camacho et al. (2009) BMC Bioinformatics 10:421
#============================================================

#------------------------------------------------------------
# STEP 1: Download the LOB domain HMM profile
#------------------------------------------------------------

# The LOB domain (PF03195) defines the LBD gene family.
# Download the profile from Pfam/InterPro:
wget -O LOB.hmm "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/PF03195?annotation=hmm"
hmmpress LOB.hmm

#------------------------------------------------------------
# STEP 2: Prepare target proteomes
#------------------------------------------------------------

# Download whole-proteome FASTA files for each species.
# Sources: Phytozome (Q. rubra), NCBI (Q. suber, others)
#
# Required files:
#   Qr_protein.fa   — Q. rubra proteome
#   Qs_protein.faa  — Q. suber proteome

#------------------------------------------------------------
# STEP 3: Run HMMER search
#------------------------------------------------------------

# hmmsearch scans each proteome for proteins matching
# the LOB domain profile. The E-value threshold of 1e-5
# balances sensitivity and specificity.

# Q. rubra
hmmsearch --tblout Qr_LBD.tbl \
          --domtblout Qr_LBD_dom.tbl \
          -E 1e-5 \
          LOB.hmm \
          Qr_protein.fa > Qr_LBD_full.out

# Q. suber
hmmsearch --tblout Qs_LBD.tbl \
          --domtblout Qs_LBD_dom.tbl \
          -E 1e-5 \
          LOB.hmm \
          Qs_protein.faa > Qs_LBD_full.out

# Output format (--tblout):
#   Column 1: target protein ID
#   Column 3: query profile name (LOB)
#   Column 5: full-sequence E-value
#   Column 6: full-sequence score

#------------------------------------------------------------
# STEP 4: Extract candidate gene IDs
#------------------------------------------------------------

# Parse the table output to get unique protein IDs
# Skip comment lines (starting with #)

grep -v "^#" Qr_LBD.tbl | awk '{print $1}' | sort -u > Qr_LBD_candidates.txt
grep -v "^#" Qs_LBD.tbl | awk '{print $1}' | sort -u > Qs_LBD_candidates.txt

echo "Q. rubra candidates: $(wc -l < Qr_LBD_candidates.txt)"
echo "Q. suber candidates: $(wc -l < Qs_LBD_candidates.txt)"

#------------------------------------------------------------
# STEP 5: Extract candidate protein sequences
#------------------------------------------------------------

# Pull full-length sequences for all candidates
# using seqtk or a simple Python/awk extraction

python3 -c "
import sys
ids = set(open('Qr_LBD_candidates.txt').read().split())
write = False
with open('Qr_protein.fa') as f:
    for line in f:
        if line.startswith('>'):
            pid = line[1:].split()[0]
            write = pid in ids
        if write:
            sys.stdout.write(line)
" > QrLBD.faa

python3 -c "
import sys
ids = set(open('Qs_LBD_candidates.txt').read().split())
write = False
with open('Qs_protein.faa') as f:
    for line in f:
        if line.startswith('>'):
            pid = line[1:].split()[0]
            write = pid in ids
        if write:
            sys.stdout.write(line)
" > QsLBD.faa

echo "QrLBD sequences: $(grep -c '^>' QrLBD.faa)"
echo "QsLBD sequences: $(grep -c '^>' QsLBD.faa)"

#------------------------------------------------------------
# STEP 6: Cross-verify with BLASTP against known LBDs
#------------------------------------------------------------

# Use all 43 known A. thaliana LBD proteins as queries
# against the oak proteomes to catch any HMMER misses.

makeblastdb -in Qr_protein.fa -dbtype prot -out Qr_db
blastp -query AtLBD_43_proteins.fa \
       -db Qr_db \
       -out Qr_blastp_vs_AtLBD.txt \
       -evalue 1e-10 \
       -outfmt 6 \
       -num_alignments 5

# Extract unique oak hits from BLAST
awk '{print $2}' Qr_blastp_vs_AtLBD.txt | sort -u > Qr_blast_hits.txt

# Compare HMMER vs BLAST candidates
comm -23 Qr_blast_hits.txt Qr_LBD_candidates.txt > Qr_blast_only.txt
echo "Additional candidates found by BLAST only: $(wc -l < Qr_blast_only.txt)"

# Repeat for Q. suber
makeblastdb -in Qs_protein.faa -dbtype prot -out Qs_db
blastp -query AtLBD_43_proteins.fa \
       -db Qs_db \
       -out Qs_blastp_vs_AtLBD.txt \
       -evalue 1e-10 \
       -outfmt 6 \
       -num_alignments 5

#------------------------------------------------------------
# STEP 7: Validate with NCBI CDD
#------------------------------------------------------------

# Submit candidate sequences to NCBI Conserved Domain Database
# (https://www.ncbi.nlm.nih.gov/cdd/) for batch search.
# Retain only proteins confirmed to contain the LOB domain.
#
# This step is performed via the NCBI web interface:
#   1. Upload QrLBD.faa and QsLBD.faa
#   2. Run CD-Search with default parameters
#   3. Download results and filter for LOB domain hits
#   4. Remove any false positives lacking the LOB domain

#------------------------------------------------------------
# STEP 8: Verify with Pfam and SMART
#------------------------------------------------------------

# Additional cross-verification using:
#   Pfam:  https://pfam.xfam.org/search
#   SMART: https://smart.embl-heidelberg.de/
#
# Proteins confirmed by at least two databases are retained
# as final LBD family members.

#------------------------------------------------------------
# SUMMARY
#------------------------------------------------------------

# Final counts after all validation steps:
#   Q. rubra: 44 LBD genes (QrLBD1-QrLBD44)
#   Q. suber: 37 LBD genes (QsLBD1-QsLBD37)
#
# All confirmed members possess intact LOB domains
# containing the C-block (CX2CX6CX3C), GAS block,
# and leucine-zipper motif (Class I only).
