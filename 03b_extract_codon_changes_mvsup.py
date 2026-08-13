#!/usr/bin/env python3
"""
03b_extract_codon_changes_mvsup.py

M. superbum-specific version of step 3-4 (see Methods, section A.7).

M. superbum has no public reference, so its pre-computed SnpEff annotation
tables were supplied as Excel sheets (columns: Allele, Annotation, Impact,
Gene_Name, Gene_ID, Feature_Type, Transcript_ID, Codon_Change, AA_Change,
CHROM, POS, REF, ALT), where Gene_ID looks like "g8129.t1.CDS6" -- i.e. one
row per individual CDS *segment* of a multi-exon gene, not per transcript.

IMPORTANT: the "Codon_Change" (HGVS.c-style, e.g. "c.2229T>C") column in
these tables was discovered to be numbered *locally within each CDS
segment* (e.g. relative to the start of "CDS6"), not across the full
spliced transcript. Using that position directly to index into the full
transcript sequence gives the wrong codon for any segment other than the
first one in a gene. This was only caught by the HGVS.p cross-validation
in 03_extract_codon_changes.py failing for the majority of sites -- see
Methods A.7 for the full account.

The fix used here: ignore the local Codon_Change position entirely and
instead map the genomic CHROM/POS + REF/ALT columns (also present in the
same tables) onto the transcript using the AUGUSTUS CDS-segment coordinates
directly, via codon_utils.genomic_pos_to_transcript_pos(). Even after this
fix, a meaningful fraction of M. superbum sites still fail HGVS.p cross-
validation (most likely because the gene-model files available do not
exactly match whichever iteration of the AUGUSTUS annotation was used to
originally run SnpEff on this data) -- those sites are dropped, and only
the surviving, independently validated subset should be reported.

Requires openpyxl (`pip install openpyxl`).

Usage:
    python 03b_extract_codon_changes_mvsup.py \\
        --rscu mvsup_rscu.pkl \\
        --tx-info mvsup_tx_info.pkl \\
        --xlsx "6D"=POSTDNA6D-6D_STRAIN_EDITED_SITES.xlsx:"POST6D-AGTC-6DSNPEFF" \\
        --xlsx "6P"=ALL_6P_MVDP-POSTDNASEQ-BOTH-REPLICATES.xlsx:GENIC \\
        --out mvsup_codon_changes.pkl

    (--tx-info is produced as a side effect of 02_build_rscu_table.py when
    run with --mode genome; see that script's build_transcript_index() /
    save it separately if needed, or regenerate here from the same BED +
    genome FASTA with codon_utils.parse_bed_cds / build_transcript_index.)
"""

import argparse
import pickle
import re

import openpyxl

from codon_utils import CODON_TABLE, parse_hgvs_p, genomic_pos_to_transcript_pos, complement_base

GENE_ID_RE = re.compile(r'(g\d+\.t\d+)')


def extract_and_validate(xlsx_path, sheet_name, label, sequences, rscu, tx_info):
    wb = openpyxl.load_workbook(xlsx_path, read_only=True, data_only=True)
    ws = wb[sheet_name]

    results = []
    header = None
    n_checked = n_mismatch = n_not_found = n_base_mismatch = 0

    for i, row in enumerate(ws.iter_rows(values_only=True)):
        if i == 0:
            header = row
            continue
        d = dict(zip(header, row))
        effect = d.get('Annotation')
        gene_id = d.get('Gene_ID')
        aa_change = d.get('AA_Change')
        chrom, pos, gref, galt = d.get('CHROM'), d.get('POS'), d.get('REF'), d.get('ALT')
        if effect not in ('synonymous_variant', 'missense_variant') or not gene_id or pos is None:
            continue
        m = GENE_ID_RE.match(str(gene_id))
        if not m:
            continue
        tid = m.group(1)

        mapped = genomic_pos_to_transcript_pos(tx_info, tid, int(pos))
        if not mapped:
            n_not_found += 1
            continue
        abs_pos, strand = mapped

        seq = sequences.get(tid)
        if not seq or abs_pos < 1 or abs_pos > len(seq):
            n_not_found += 1
            continue

        gref, galt = str(gref).upper(), str(galt).upper()
        if len(gref) != 1 or len(galt) != 1 or gref not in 'ACGT' or galt not in 'ACGT':
            continue
        tref, talt = (gref, galt) if strand == '+' else (complement_base(gref), complement_base(galt))

        if seq[abs_pos - 1] != tref:
            n_base_mismatch += 1
            continue

        codon_start = ((abs_pos - 1) // 3) * 3
        ref_codon = seq[codon_start:codon_start + 3]
        pos_in_codon = (abs_pos - 1) - codon_start
        alt_codon = ref_codon[:pos_in_codon] + talt + ref_codon[pos_in_codon + 1:]
        if len(ref_codon) != 3 or ref_codon not in CODON_TABLE or alt_codon not in CODON_TABLE:
            continue

        derived_ref_aa = CODON_TABLE[ref_codon]
        derived_alt_aa = CODON_TABLE[alt_codon]

        parsed_p = parse_hgvs_p(aa_change) if aa_change else None
        if parsed_p:
            n_checked += 1
            p_ref, _, p_alt = parsed_p
            if not (p_ref == derived_ref_aa and p_alt == derived_alt_aa):
                n_mismatch += 1
                continue

        results.append({
            'label': label, 'tid': tid, 'chrom': chrom, 'pos': pos, 'effect': effect,
            'ref_codon': ref_codon, 'alt_codon': alt_codon,
            'ref_rscu': rscu.get(ref_codon), 'alt_rscu': rscu.get(alt_codon),
            'aa_change': aa_change,
        })

    print(f"{label}: kept={len(results)}  not_found={n_not_found}  base_mismatch={n_base_mismatch}  "
          f"hgvsp_checked={n_checked}  hgvsp_mismatch={n_mismatch}")
    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--rscu', required=True, help='pickle from 02_build_rscu_table.py (mode=genome, M. superbum)')
    ap.add_argument('--tx-info', required=True,
                     help='pickle containing the tx_info dict from codon_utils.build_transcript_index()')
    ap.add_argument('--xlsx', action='append', required=True, metavar='LABEL=PATH:SHEET')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    with open(args.rscu, 'rb') as f:
        rscu_data = pickle.load(f)
    sequences, rscu = rscu_data['sequences'], rscu_data['rscu']

    with open(args.tx_info, 'rb') as f:
        tx_info = pickle.load(f)

    all_results = {}
    for item in args.xlsx:
        label, rest = item.split('=', 1)
        path, sheet = rest.rsplit(':', 1)
        all_results[label] = extract_and_validate(path, sheet, label, sequences, rscu, tx_info)

    with open(args.out, 'wb') as f:
        pickle.dump(all_results, f)
    print(f"Saved: {args.out}")


if __name__ == '__main__':
    main()
