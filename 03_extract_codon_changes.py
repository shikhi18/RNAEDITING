#!/usr/bin/env python3
"""
03_extract_codon_changes.py

Steps 3-4 of the codon-optimality analysis (see Methods, sections A.5-A.6).

For each synonymous (and, optionally, missense) A-to-I site, extracts the
exact ref/alt codon by indexing directly into the transcript's spliced CDS
sequence at its HGVS.c coding-sequence position, looks up RSCU for both
codons, and classifies the direction of change ("faster" = toward the more
common codon, "slower" = toward the rarer one).

Every extracted codon is cross-validated against SnpEff's own HGVS.p amino-
acid call for the same site: if translating ref_codon/alt_codon doesn't
reproduce the same amino-acid change SnpEff already reported, the site is
dropped rather than silently kept. This is the check that caught the
M. superbum position-numbering issue described in the Methods (section A.7)
-- if you see a high mismatch rate here, the transcript sequences and the
VCF's HGVS.c numbering are not on the same coordinate system and should not
be trusted without investigating further (see 03b_extract_codon_changes_mvsup.py
for how that was resolved for M. superbum).

This script assumes the *transcript* HGVS.c position in the VCF is numbered
across the full spliced transcript (true for standard SnpEff runs, and
confirmed true for the M. intermedium and M. lychnidis-dioicae VCFs used in
this study). It was NOT true for the M. superbum SnpEff tables -- use
03b_extract_codon_changes_mvsup.py for that species instead.

Usage:
    python 03_extract_codon_changes.py \\
        --rscu mi_rscu.pkl \\
        --vcf "MI A1"=raw_snps-filtered.ann.vcf \\
        --vcf "MI A2"=raw_snps979899-filtered.ann.vcf \\
        --vcf "MI Mated"=raw_snpsM1M2M3-filtered.ann.vcf \\
        --out mi_codon_changes.pkl
"""

import argparse
import pickle

from codon_utils import iter_ann_records, parse_hgvs_p, CODON_TABLE


def extract_and_validate(vcf_path, label, sequences, rscu, effects=("synonymous_variant", "missense_variant")):
    results = []
    n_checked = 0
    n_hgvsp_mismatch = 0
    n_seq_missing = 0

    for rec in iter_ann_records(vcf_path, effects=effects):
        tid = rec['feature_id']
        cds_pos, tref, talt = rec['cds_pos'], rec['tref'], rec['talt']
        if cds_pos is None:
            continue
        seq = sequences.get(tid)
        if not seq:
            n_seq_missing += 1
            continue
        if cds_pos < 1 or cds_pos > len(seq):
            continue
        if seq[cds_pos - 1] != tref:
            continue  # reference-base mismatch: sequence/position system disagree, skip
        codon_start = ((cds_pos - 1) // 3) * 3
        ref_codon = seq[codon_start:codon_start + 3]
        pos_in_codon = (cds_pos - 1) - codon_start
        alt_codon = ref_codon[:pos_in_codon] + talt + ref_codon[pos_in_codon + 1:]
        if len(ref_codon) != 3 or ref_codon not in CODON_TABLE or alt_codon not in CODON_TABLE:
            continue

        derived_ref_aa = CODON_TABLE[ref_codon]
        derived_alt_aa = CODON_TABLE[alt_codon]

        parsed_p = parse_hgvs_p(rec['hgvs_p'])
        if parsed_p:
            n_checked += 1
            p_ref, _, p_alt = parsed_p
            if not (p_ref == derived_ref_aa and p_alt == derived_alt_aa):
                n_hgvsp_mismatch += 1
                continue  # failed independent validation against SnpEff's own call -- drop

        results.append({
            'label': label, 'tid': tid, 'chrom': rec['chrom'], 'pos': rec['pos'],
            'effect': rec['effect'], 'ref_codon': ref_codon, 'alt_codon': alt_codon,
            'ref_rscu': rscu.get(ref_codon), 'alt_rscu': rscu.get(alt_codon),
            'hgvs_p': rec['hgvs_p'],
        })

    print(f"{label}: kept={len(results)}  seq_missing={n_seq_missing}  "
          f"hgvsp_checked={n_checked}  hgvsp_mismatch={n_hgvsp_mismatch}")
    return results


def summarize_direction(results, effect='synonymous_variant'):
    faster = slower = same = 0
    for r in results:
        if r['effect'] != effect:
            continue
        rr, ar = r['ref_rscu'], r['alt_rscu']
        if rr is None or ar is None:
            continue
        if ar > rr + 1e-9:
            faster += 1
        elif ar < rr - 1e-9:
            slower += 1
        else:
            same += 1
    total = faster + slower + same
    if total:
        print(f"  {effect}: n={total}  toward more common codon={faster} ({100*faster/total:.1f}%)  "
              f"toward less common codon={slower} ({100*slower/total:.1f}%)  unchanged={same}")
    return faster, slower, same


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--rscu', required=True, help='pickle produced by 02_build_rscu_table.py')
    ap.add_argument('--vcf', action='append', required=True, metavar='LABEL=PATH',
                     help='repeatable; one SnpEff ANN VCF per condition')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    with open(args.rscu, 'rb') as f:
        rscu_data = pickle.load(f)
    sequences = rscu_data['sequences']
    rscu = rscu_data['rscu']

    all_results = {}
    for item in args.vcf:
        label, path = item.split('=', 1)
        results = extract_and_validate(path, label, sequences, rscu)
        summarize_direction(results, 'synonymous_variant')
        all_results[label] = results

    with open(args.out, 'wb') as f:
        pickle.dump(all_results, f)
    print(f"Saved: {args.out}")


if __name__ == '__main__':
    main()
