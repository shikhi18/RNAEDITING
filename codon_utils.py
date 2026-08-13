"""
codon_utils.py

Shared utilities for the codon-optimality and structural-position analysis
of A-to-I RNA editing sites (M. intermedium, M. lychnidis-dioicae, M. superbum).

Contains:
  - the standard genetic code
  - FASTA / GFF3 / GTF / BED parsing helpers
  - CDS splicing from genome + CDS-segment coordinates
  - a SnapGene .dna binary sequence parser (used for M. superbum, which has
    no public reference genome; per-transcript coding sequences were
    pre-exported to SnapGene .dna files)
  - RSCU (Relative Synonymous Codon Usage) table construction
  - SnpEff ANN-field parsing helpers
"""

import re
import struct
from collections import defaultdict

# ---------------------------------------------------------------------------
# Standard genetic code
# ---------------------------------------------------------------------------

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

AA3TO1 = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Gln': 'Q', 'Glu': 'E', 'Gly': 'G',
    'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S',
    'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', 'Ter': '*',
}

_COMPLEMENT = str.maketrans('ACGTacgt', 'TGCAtgca')
_COMPLEMENT_BASE = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}


def revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def complement_base(base):
    return _COMPLEMENT_BASE[base.upper()]


# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------

def load_fasta(path):
    """Return {sequence_id: sequence} for a (multi-)FASTA file. sequence_id
    is the first whitespace-delimited token of the header line."""
    seqs = {}
    cur_id, cur_chunks = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if cur_id is not None:
                    seqs[cur_id] = ''.join(cur_chunks)
                cur_id = line[1:].split()[0]
                cur_chunks = []
            else:
                cur_chunks.append(line)
        if cur_id is not None:
            seqs[cur_id] = ''.join(cur_chunks)
    return seqs


# ---------------------------------------------------------------------------
# SnapGene .dna binary sequence parser
# ---------------------------------------------------------------------------

def parse_snapgene_dna(path):
    """Extract the raw nucleotide sequence from a SnapGene .dna file.

    SnapGene .dna files are a sequence of blocks: [1-byte type][4-byte
    big-endian length][length bytes of data]. Block type 0 holds the DNA
    sequence: 1 flag byte (topology) followed by the ASCII sequence.
    Returns the sequence in upper case, or None if no sequence block found.
    """
    with open(path, 'rb') as f:
        data = f.read()
    pos = 0
    while pos + 5 <= len(data):
        block_type = data[pos]
        block_len = struct.unpack('>I', data[pos + 1:pos + 5])[0]
        block_data = data[pos + 5:pos + 5 + block_len]
        if block_type == 0 and len(block_data) >= 1:
            return block_data[1:].decode('ascii', errors='replace').upper()
        pos += 5 + block_len
    return None


# ---------------------------------------------------------------------------
# CDS splicing from a genome FASTA + CDS-segment coordinates (BED/GFF3-style)
# ---------------------------------------------------------------------------

def build_transcript_index(cds_segments_by_tx):
    """Given {transcript_id: [(seg_start_0based, seg_end_0based_exclusive,
    strand, contig), ...]}, return {transcript_id: (strand, contig,
    [(seg_start, seg_end, cumulative_offset_before_this_segment), ...])}
    with segments ordered in transcript (5'->3') direction and cumulative
    offsets in transcript-relative, 0-based coordinates.
    """
    tx_info = {}
    for tid, segs in cds_segments_by_tx.items():
        strand = segs[0][2]
        contig = segs[0][3]
        ordered = sorted(segs, key=lambda x: x[0]) if strand == '+' else sorted(segs, key=lambda x: -x[0])
        cum = 0
        seglist = []
        for seg_start, seg_end, _, _ in ordered:
            seglist.append((seg_start, seg_end, cum))
            cum += (seg_end - seg_start)
        tx_info[tid] = (strand, contig, seglist)
    return tx_info


def splice_transcript(genome, tx_info, tid):
    """Build the spliced CDS sequence (5'->3', in reading frame from the
    first base of the CDS) for one transcript, given the genome dict and
    the tx_info produced by build_transcript_index()."""
    info = tx_info.get(tid)
    if not info:
        return None
    strand, contig, seglist = info
    if contig not in genome:
        return None
    parts = []
    for seg_start, seg_end, _ in sorted(seglist, key=lambda x: x[0]):
        parts.append(genome[contig][seg_start:seg_end])
    seq = ''.join(parts)
    if strand == '-':
        seq = revcomp(seq)
    return seq.upper()


def genomic_pos_to_transcript_pos(tx_info, tid, pos_1based):
    """Map a 1-based genomic position to a 1-based position within the
    fully spliced transcript, using the transcript's CDS segment map.
    Returns (abs_transcript_pos_1based, strand) or None if the position
    doesn't fall within any CDS segment of this transcript."""
    info = tx_info.get(tid)
    if not info:
        return None
    strand, contig, seglist = info
    pos0 = pos_1based - 1
    for seg_start, seg_end, cum in seglist:
        if seg_start <= pos0 < seg_end:
            offset = (pos0 - seg_start) if strand == '+' else ((seg_end - 1) - pos0)
            return cum + offset + 1, strand
    return None


# ---------------------------------------------------------------------------
# RSCU (Relative Synonymous Codon Usage)
# ---------------------------------------------------------------------------

def build_rscu_table(cds_sequences):
    """Given an iterable of CDS sequences (each assumed in-frame from base
    1), return (codon_counts, rscu) where rscu[codon] = observed count /
    expected count under equal synonymous usage for that amino acid.
    RSCU > 1 = used more often than expected ("optimal"/"faster");
    RSCU < 1 = used less often ("rare"/"slower")."""
    codon_counts = defaultdict(int)
    for seq in cds_sequences:
        seq = seq.upper()
        trim = len(seq) - (len(seq) % 3)
        for i in range(0, trim, 3):
            codon = seq[i:i + 3]
            if all(b in 'ACGT' for b in codon):
                codon_counts[codon] += 1

    aa_to_codons = defaultdict(list)
    for codon, aa in CODON_TABLE.items():
        aa_to_codons[aa].append(codon)

    rscu = {}
    for aa, codons in aa_to_codons.items():
        counts = [codon_counts.get(c, 0) for c in codons]
        total = sum(counts)
        n = len(codons)
        if total == 0:
            continue
        for c, cnt in zip(codons, counts):
            rscu[c] = cnt / (total / n)
    return dict(codon_counts), rscu


def gc3_content(codon_counts):
    """Fraction of codons with G or C at the third (wobble) position."""
    gc3 = sum(v for c, v in codon_counts.items() if c[2] in 'GC')
    total = sum(codon_counts.values())
    return gc3 / total if total else None


# ---------------------------------------------------------------------------
# SnpEff ANN-field parsing (VCF INFO string)
# ---------------------------------------------------------------------------

ANN_FIELDS = [
    "allele", "effect", "impact", "gene_name", "gene_id", "feature_type",
    "feature_id", "biotype", "rank", "hgvs_c", "hgvs_p", "cdna_pos",
    "cds_pos", "aa_pos", "distance", "errors",
]

_HGVS_C_RE = re.compile(r'c\.(-?\d+)([ACGT])>([ACGT])')
_HGVS_P_RE = re.compile(r'p\.([A-Za-z]{3})(\d+)([A-Za-z]{3}|\*)')
_ANN_RE = re.compile(r'ANN=([^;]+)')


def iter_ann_records(vcf_path, effects=("synonymous_variant", "missense_variant"),
                      pass_only=True, editing_only=True):
    """Yield dicts for each PASS, A-to-I (A>G or T>C), single-annotation
    record in a SnpEff ANN-annotated VCF whose effect is in `effects`.
    Only the first matching annotation per VCF line is yielded (a variant
    overlapping multiple transcripts produces multiple ANN entries; take
    the first genic one, matching the approach used throughout this
    analysis).

    Yielded dict keys: chrom, pos, ref, alt, effect, feature_id (transcript
    ID, "transcript:" prefix stripped), hgvs_c, hgvs_p, cds_pos (int or
    None), tref, talt (transcript-relative ref/alt bases parsed from
    hgvs_c).
    """
    with open(vcf_path, errors='ignore') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 8:
                continue
            chrom, pos, _id, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
            filt = fields[6]
            if pass_only and filt != 'PASS':
                continue
            if editing_only:
                if len(ref) != 1 or len(alt) != 1:
                    continue
                if not ((ref == 'A' and alt == 'G') or (ref == 'T' and alt == 'C')):
                    continue
            m = _ANN_RE.search(fields[7])
            if not m:
                continue
            for ann in m.group(1).split(','):
                parts = ann.split('|')
                if len(parts) < 11:
                    continue
                rec = dict(zip(ANN_FIELDS, parts))
                if rec['effect'] not in effects:
                    continue
                feature_id = rec['feature_id'].replace('transcript:', '')
                cm = _HGVS_C_RE.search(rec['hgvs_c'])
                cds_pos, tref, talt = (None, None, None)
                if cm:
                    cds_pos, tref, talt = int(cm.group(1)), cm.group(2), cm.group(3)
                yield {
                    'chrom': chrom, 'pos': int(pos), 'ref': ref, 'alt': alt,
                    'effect': rec['effect'], 'feature_id': feature_id,
                    'hgvs_c': rec['hgvs_c'], 'hgvs_p': rec['hgvs_p'],
                    'cds_pos': cds_pos, 'tref': tref, 'talt': talt,
                }
                break  # only the first matching annotation per record


def parse_hgvs_p(hgvs_p):
    """Parse a 'p.Ser743Ser'-style string to (ref_aa_1letter, aa_pos,
    alt_aa_1letter) or None if it doesn't match."""
    m = _HGVS_P_RE.match(str(hgvs_p))
    if not m:
        return None
    ref3, aa_pos, alt3 = m.group(1), int(m.group(2)), m.group(3)
    ref1 = AA3TO1.get(ref3, '?')
    alt1 = '*' if alt3 == '*' else AA3TO1.get(alt3, '?')
    return ref1, aa_pos, alt1


def codon_position(cds_pos):
    """1-based position within its codon (1, 2, or 3) for a 1-based
    coding-sequence position."""
    return ((cds_pos - 1) % 3) + 1
