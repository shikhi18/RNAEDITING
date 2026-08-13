# Codon Optimality and Structural Position Analysis

Scripts used to test (1) whether A-to-I RNA-editing sites shift codons toward
faster/more-common or slower/rarer synonymous codons, and (2) where
recurring nonsynonymous editing sites fall relative to known protein
domains, topology, and disordered regions, for *Microbotryum intermedium*
(MI), *M. lychnidis-dioicae* (MVLG), and *M. superbum* (MvSup).

This accompanies the corresponding Methods write-up ("Methods: Codon
Optimality and Structural Position Analysis of A-to-I Editing Sites").

## Requirements

```
pip install openpyxl requests
```

Python 3.8+. No other third-party dependencies (FASTA/GFF3/BED/VCF parsing
is done with the standard library in `codon_utils.py`).

## Pipeline overview

| Script | Purpose | Methods section |
|---|---|---|
| `codon_utils.py` | Shared library: genetic code, FASTA/GFF3/BED parsing, CDS splicing, SnapGene `.dna` parser, RSCU calculation, SnpEff ANN-field parsing | — |
| `01_wobble_position_check.py` | Confirms synonymous edits fall at the 3rd (wobble) codon position | A.3 |
| `02_build_rscu_table.py` | Builds a genome-wide codon usage (RSCU) table + per-transcript CDS sequences | A.4 |
| `03_extract_codon_changes.py` | Extracts exact ref/alt codons at synonymous sites, classifies direction of RSCU change, cross-validates against SnpEff's HGVS.p (MI, MVLG) | A.5, A.6 |
| `03b_extract_codon_changes_mvsup.py` | Same as above but for M. superbum, using genomic-position mapping instead of the (locally-numbered, unreliable) HGVS.c position in that species' SnpEff tables | A.7 |
| `04_identify_recurring_nonsynonymous_sites.py` | Finds nonsynonymous sites detected independently in 2+ conditions | B.3 |
| `05_fetch_domain_annotations.py` | Retrieves Pfam/PROSITE/PANTHER/Gene3D/SuperFamily/Phobius/TMHMM/disorder annotations for a transcript from Ensembl (+ UniProt cross-reference) | B.4 |

## Typical run (M. lychnidis-dioicae example)

```bash
# 1. Confirm wobble-position bias (no reference needed, just the VCFs)
python 01_wobble_position_check.py \
    "p1A1"=raw_snpsp1A1-filtered.ann.vcf \
    "p1A2"=raw_snpsp1A2-filtered.ann.vcf \
    "Mated"=raw_snpsmated-filtered.ann.vcf

# 2. Build the genome-wide codon usage table (splices CDS from genome + GFF3)
python 02_build_rscu_table.py --mode genome \
    --fasta MVLG.genome.fa --annotation MVLG.gff3 --format gff3 \
    --out mvlg_rscu.pkl

# 3. Extract + validate exact codon changes, classify fast/slow direction
python 03_extract_codon_changes.py \
    --rscu mvlg_rscu.pkl \
    --vcf "p1A1"=raw_snpsp1A1-filtered.ann.vcf \
    --vcf "p1A2"=raw_snpsp1A2-filtered.ann.vcf \
    --vcf "Mated"=raw_snpsmated-filtered.ann.vcf \
    --out mvlg_codon_changes.pkl

# 4. Find recurring nonsynonymous sites for structural follow-up
python 04_identify_recurring_nonsynonymous_sites.py \
    "p1A1"=raw_snpsp1A1-filtered.ann.vcf \
    "p1A2"=raw_snpsp1A2-filtered.ann.vcf \
    "Mated"=raw_snpsmated-filtered.ann.vcf

# 5. Pull domain/topology context for the strongest candidate genes
python 05_fetch_domain_annotations.py MVLG_00116T0 MVLG_02820T0 MVLG_06756T0
```

For M. intermedium, run step 2 with `--mode cdna --fasta MI.cdna.all.fa`
(a pre-spliced reference is used directly, since it's already indexed by
the same transcript IDs SnpEff reports).

For M. superbum, run step 2 with `--mode genome --format bed
--save-tx-info mvsup_tx_info.pkl`, then use
`03b_extract_codon_changes_mvsup.py` instead of `03_...py` — see the
docstring at the top of that script for why (M. superbum's SnpEff
Codon_Change field uses a non-standard, per-CDS-segment-local position
numbering that does not match the full spliced transcript, discovered via
the HGVS.p cross-validation step failing for the majority of sites until
this was worked around).


