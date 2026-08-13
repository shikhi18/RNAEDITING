#!/usr/bin/env python3
"""
04_identify_recurring_nonsynonymous_sites.py

Step 1 of the structural-position analysis (see Methods, section B.3).

Scans SnpEff ANN-annotated VCFs across multiple biological conditions
(e.g. two haploid mating types and a mated/dikaryotic stage) for
missense_variant A-to-I sites, and finds the ones observed independently
in two or more conditions -- these recurring sites are the most robust
individual editing events to examine further (e.g. for structural context,
see 05_fetch_domain_annotations.py).

Usage:
    python 04_identify_recurring_nonsynonymous_sites.py \\
        "A1"=raw_snps-filtered.ann.vcf \\
        "A2"=raw_snps979899-filtered.ann.vcf \\
        "Mated"=raw_snpsM1M2M3-filtered.ann.vcf
"""

import sys
from collections import defaultdict

from codon_utils import iter_ann_records


def find_recurring_sites(vcf_paths_by_label, min_conditions=2):
    """vcf_paths_by_label: {label: vcf_path}. Returns:
      recurring: {(transcript_id, chrom, pos): [label, ...]}
      per_site_detail: {(transcript_id, chrom, pos): [(label, hgvs_c, hgvs_p), ...]}
    restricted to sites seen in >= min_conditions distinct labels.
    """
    by_site = defaultdict(list)
    detail = defaultdict(list)
    for label, path in vcf_paths_by_label.items():
        for rec in iter_ann_records(path, effects=("missense_variant",)):
            key = (rec['feature_id'], rec['chrom'], rec['pos'])
            by_site[key].append(label)
            detail[key].append((label, rec['hgvs_c'], rec['hgvs_p']))

    recurring = {k: v for k, v in by_site.items() if len(set(v)) >= min_conditions}
    recurring_detail = {k: detail[k] for k in recurring}
    return recurring, recurring_detail


def rank_genes_by_recurring_site_count(recurring):
    counts = defaultdict(int)
    for (tid, chrom, pos) in recurring:
        counts[tid] += 1
    return sorted(counts.items(), key=lambda x: -x[1])


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    vcf_paths_by_label = {}
    for arg in sys.argv[1:]:
        label, path = arg.split('=', 1)
        vcf_paths_by_label[label] = path

    recurring, detail = find_recurring_sites(vcf_paths_by_label)
    print(f"Sites recurring in >=2 conditions: {len(recurring)}\n")

    ranked = rank_genes_by_recurring_site_count(recurring)
    print("Top genes by number of recurring nonsynonymous sites:")
    for tid, n in ranked[:20]:
        print(f"  {tid}: {n} recurring site(s)")

    print("\nDetail for top gene:")
    top_tid = ranked[0][0]
    for (tid, chrom, pos), records in sorted(detail.items()):
        if tid != top_tid:
            continue
        labels_seen = sorted(set(l for l, _, _ in records))
        _, hgvs_c, hgvs_p = records[0]
        print(f"  {chrom}:{pos}  {hgvs_c}  {hgvs_p}  (seen in: {', '.join(labels_seen)})")
