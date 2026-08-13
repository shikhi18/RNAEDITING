#!/usr/bin/env python3
"""
01_wobble_position_check.py

Step 1 of the codon-optimality analysis (see Methods, section A.3).

For every PASS, A-to-I (A>G / T>C) SnpEff-annotated site, determines which
position within its codon (1st, 2nd, or 3rd / wobble) the edited base falls
on, and tabulates the distribution separately for synonymous vs.
missense (nonsynonymous) sites. This uses only the HGVS.c coding-sequence
position already present in the VCF's SnpEff ANN field -- no reference
FASTA is required for this step.

Usage:
    python 01_wobble_position_check.py <label1>=<vcf1> [<label2>=<vcf2> ...]

Example:
    python 01_wobble_position_check.py \\
        "MI A1"=raw_snps-filtered.ann.vcf \\
        "MI A2"=raw_snps979899-filtered.ann.vcf \\
        "MI Mated"=raw_snpsM1M2M3-filtered.ann.vcf
"""

import sys
from collections import defaultdict

from codon_utils import iter_ann_records, codon_position


def analyze_vcf(vcf_path, label):
    syn_pos = defaultdict(int)
    nonsyn_pos = defaultdict(int)
    for rec in iter_ann_records(vcf_path, effects=("synonymous_variant", "missense_variant")):
        if rec['cds_pos'] is None:
            continue
        pic = codon_position(rec['cds_pos'])
        if rec['effect'] == 'synonymous_variant':
            syn_pos[pic] += 1
        else:
            nonsyn_pos[pic] += 1

    syn_total = sum(syn_pos.values())
    nonsyn_total = sum(nonsyn_pos.values())
    print(f"{label}:")
    if syn_total:
        print(f"  synonymous   n={syn_total:6d}  pos1={syn_pos[1]:5d} ({100*syn_pos[1]/syn_total:5.1f}%)  "
              f"pos2={syn_pos[2]:5d} ({100*syn_pos[2]/syn_total:5.1f}%)  "
              f"pos3={syn_pos[3]:5d} ({100*syn_pos[3]/syn_total:5.1f}%)")
    if nonsyn_total:
        print(f"  nonsynonymous n={nonsyn_total:6d}  pos1={nonsyn_pos[1]:5d} ({100*nonsyn_pos[1]/nonsyn_total:5.1f}%)  "
              f"pos2={nonsyn_pos[2]:5d} ({100*nonsyn_pos[2]/nonsyn_total:5.1f}%)  "
              f"pos3={nonsyn_pos[3]:5d} ({100*nonsyn_pos[3]/nonsyn_total:5.1f}%)")
    return syn_pos, nonsyn_pos


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    for arg in sys.argv[1:]:
        label, path = arg.split('=', 1)
        analyze_vcf(path, label)
