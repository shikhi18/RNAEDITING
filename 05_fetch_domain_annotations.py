#!/usr/bin/env python3
"""
05_fetch_domain_annotations.py

Step 2 of the structural-position analysis (see Methods, section B.4).

For a given Ensembl translation/transcript ID, fetches predicted protein
domain, topology, and disorder features from the Ensembl REST API
(overlap/translation endpoint), and optionally the corresponding UniProt
TrEMBL flat-file record for cross-referencing.

Note: this analysis originally used a hosted fetch tool rather than direct
HTTP calls; this script reproduces the same two API calls with the
`requests` library so it's runnable standalone. An attempt to also
download real AlphaFold 3D structures (alphafold.ebi.ac.uk) to compute true
solvent accessibility was not successful in that environment (see Methods,
section B.4) -- this script only retrieves linear domain/topology/disorder
annotations, which is what the structural-position findings in the paper
are based on.

Requires `requests` (`pip install requests`). Be considerate of Ensembl's
and UniProt's REST rate limits if querying many IDs (add a short sleep
between calls; both endpoints return empty responses rather than a clear
error when rate-limited, which is why iter_features() below explicitly
checks for and reports empty responses).

Usage:
    python 05_fetch_domain_annotations.py MVLG_00116T0 MVLG_02820T0 MVLG_06756T0
"""

import sys
import time
import json

import requests

ENSEMBL_OVERLAP_URL = "https://rest.ensembl.org/overlap/translation/{id}?content-type=application/json"
UNIPROT_FLATFILE_URL = "https://rest.uniprot.org/uniprotkb/{accession}.txt"


def fetch_ensembl_domains(translation_id, retries=3, wait_s=20):
    """Return the parsed JSON list of protein_feature records for a
    translation ID, or None if the request failed / was empty after
    retrying (Ensembl's REST API returns an empty body rather than an
    HTTP error when rate-limited)."""
    url = ENSEMBL_OVERLAP_URL.format(id=translation_id)
    for attempt in range(retries):
        resp = requests.get(url, timeout=30)
        if resp.status_code == 200 and resp.text.strip():
            try:
                return resp.json()
            except json.JSONDecodeError:
                pass
        if attempt < retries - 1:
            time.sleep(wait_s)
    return None


def fetch_uniprot_record(accession):
    """Return the raw UniProt flat-file text for a TrEMBL/Swiss-Prot
    accession, or None on failure."""
    url = UNIPROT_FLATFILE_URL.format(accession=accession)
    resp = requests.get(url, timeout=30)
    if resp.status_code == 200 and resp.text.strip():
        return resp.text
    return None


def summarize_domains(features):
    """Print a compact summary of the informative (non-disorder,
    non-low-complexity) domain/topology features for one translation."""
    if not features:
        print("  (no features returned)")
        return
    informative_types = {
        'Pfam', 'PROSITE_profiles', 'Prosite_profiles', 'PROSITE_patterns', 'Prosite_patterns',
        'PANTHER', 'Gene3D', 'SuperFamily', 'CDD', 'NCBIfam', 'TIGRfam', 'Phobius', 'TMHMM',
    }
    for feat in features:
        ftype = feat.get('type')
        if ftype not in informative_types:
            continue
        desc = feat.get('description') or feat.get('hseqname')
        start, end = feat.get('start'), feat.get('end')
        interpro = feat.get('interpro')
        print(f"  [{ftype}] {feat.get('id') or feat.get('hseqname')}  {start}-{end}  "
              f"{desc}" + (f"  (InterPro {interpro})" if interpro else ""))
    alphafold = [f for f in features if f.get('type') == 'alphafold']
    if alphafold:
        print(f"  AlphaFold model available: {alphafold[0].get('id')}")
    else:
        print("  No AlphaFold model listed for this translation.")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    for translation_id in sys.argv[1:]:
        print(f"=== {translation_id} ===")
        features = fetch_ensembl_domains(translation_id)
        summarize_domains(features)
        print()
        time.sleep(15)  # be polite between requests -- Ensembl's REST API rate-limits repeated calls
