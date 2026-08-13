#!/usr/bin/env python3
"""
02_build_rscu_table.py

Step 2 of the codon-optimality analysis (see Methods, section A.4).

Builds a genome-wide codon usage / RSCU table and a set of per-transcript
CDS sequences for one species, using one of two reference styles:

  --mode cdna   : a pre-spliced cDNA/CDS FASTA (one sequence per transcript,
                  already in frame from base 1). Used for M. intermedium.

  --mode genome : a genome FASTA + a GFF3/GTF/BED file listing CDS segments
                  per transcript, which this script splices itself
                  (handling multi-exon transcripts and strand). Used for
                  M. lychnidis-dioicae (GFF3 input) and, with --feature-type
                  and --id-attr adjusted, M. superbum (BED input).

Output: a pickle containing {'codon_counts': {...}, 'rscu': {...},
'sequences': {transcript_id: seq}}, plus a summary printed to stdout
(number of transcripts used, total codons, 3rd-position GC%).

Usage examples:
    # M. intermedium (pre-spliced Ensembl cDNA FASTA)
    python 02_build_rscu_table.py --mode cdna \\
        --fasta MI.cdna.all.fa --out mi_rscu.pkl

    # M. lychnidis-dioicae (genome FASTA + GFF3 CDS features)
    python 02_build_rscu_table.py --mode genome \\
        --fasta MVLG.genome.fa --annotation MVLG.gff3 --format gff3 \\
        --out mvlg_rscu.pkl

    # M. superbum (genome FASTA + AUGUSTUS BED CDS features)
    python 02_build_rscu_table.py --mode genome \\
        --fasta MvSup.genome.fa --annotation MvSup.augustus.bed --format bed \\
        --out mvsup_rscu.pkl
"""

import argparse
import pickle
import re
from collections import defaultdict

from codon_utils import load_fasta, build_transcript_index, splice_transcript, build_rscu_table, gc3_content


def parse_gff3_cds(path):
    """Return {transcript_id: [(start0, end0, strand, contig), ...]} from a
    GFF3 file's CDS features, matched via 'Parent=transcript:<id>'."""
    cds_by_tx = defaultdict(list)
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'CDS':
                continue
            contig, start, end, strand, attrs = fields[0], int(fields[3]), int(fields[4]), fields[6], fields[8]
            m = re.search(r'Parent=transcript:([^;]+)', attrs)
            if not m:
                continue
            tid = m.group(1)
            # GFF3 coordinates are 1-based inclusive; convert to 0-based half-open
            cds_by_tx[tid].append((start - 1, end, strand, contig))
    return cds_by_tx


def parse_bed_cds(path, feature_col=7, id_attr_regex=r'Parent=([^;]+)'):
    """Return {transcript_id: [(start0, end0, strand, contig), ...]} from an
    AUGUSTUS-style extended BED file's CDS rows. `feature_col` is the
    0-based column index holding the feature type (default 7, matching
    <chrom> <start> <end> <name> <score> <strand> <source> <feature> ...).
    Contig names with a trailing '_len=...' style suffix are trimmed to
    the part before the first underscore.
    """
    cds_by_tx = defaultdict(list)
    id_re = re.compile(id_attr_regex)
    with open(path) as f:
        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) <= feature_col or fields[feature_col] != 'CDS':
                continue
            contig_full, start, end, strand = fields[0], int(fields[1]), int(fields[2]), fields[5]
            contig = contig_full.split('_')[0]
            attrs = fields[9] if len(fields) > 9 else ''
            m = id_re.search(attrs)
            if not m:
                continue
            tid = m.group(1)
            # BED coordinates are already 0-based half-open
            cds_by_tx[tid].append((start, end, strand, contig))
    return cds_by_tx


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--mode', choices=['cdna', 'genome'], required=True)
    ap.add_argument('--fasta', required=True, help='cDNA FASTA (mode=cdna) or genome FASTA (mode=genome)')
    ap.add_argument('--annotation', help='GFF3 or BED file with CDS features (mode=genome only)')
    ap.add_argument('--format', choices=['gff3', 'bed'], help='annotation file format (mode=genome only)')
    ap.add_argument('--out', required=True, help='output pickle path')
    ap.add_argument('--save-tx-info', metavar='PATH',
                     help='(mode=genome only) also save the transcript CDS-segment index here; '
                          'required as input to 03b_extract_codon_changes_mvsup.py for M. superbum')
    args = ap.parse_args()

    tx_info = None
    if args.mode == 'cdna':
        sequences = load_fasta(args.fasta)
    else:
        if not args.annotation or not args.format:
            ap.error('--annotation and --format are required when --mode genome')
        genome = load_fasta(args.fasta)
        cds_by_tx = parse_gff3_cds(args.annotation) if args.format == 'gff3' else parse_bed_cds(args.annotation)
        tx_info = build_transcript_index(cds_by_tx)
        sequences = {}
        for tid in tx_info:
            seq = splice_transcript(genome, tx_info, tid)
            if seq:
                sequences[tid] = seq
        if args.save_tx_info:
            with open(args.save_tx_info, 'wb') as f:
                pickle.dump(tx_info, f)
            print(f"Saved transcript index: {args.save_tx_info}")

    codon_counts, rscu = build_rscu_table(sequences.values())
    total_codons = sum(codon_counts.values())
    gc3 = gc3_content(codon_counts)

    print(f"Transcripts used   : {len(sequences)}")
    print(f"Total codons       : {total_codons}")
    print(f"3rd-position GC%   : {100*gc3:.1f}%" if gc3 is not None else "3rd-position GC%   : n/a")

    with open(args.out, 'wb') as f:
        pickle.dump({'codon_counts': codon_counts, 'rscu': rscu, 'sequences': sequences}, f)
    print(f"Saved: {args.out}")


if __name__ == '__main__':
    main()
