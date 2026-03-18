#!/usr/bin/env python3
"""
Benchmark validation: check ASV recovery from mock community data.

Expected: 22 reference sequences in mock_sequences_V4.fasta
Checks:
  - ASV count is >= MIN_ASVS (allowing for genuine noise ASVs above the reference)
  - ASV count is <= MAX_ASVS (not exploding with spurious sequences)
  - All 22 reference sequences have at least one near-identical ASV (>= IDENTITY)
"""

import sys
import os
import argparse

MIN_ASVS   = 21   # allow 1 ASV to be missed
MAX_ASVS   = 23   # allow 1 extra ASV above the 22 references
IDENTITY   = 0.97 # minimum identity to consider a reference sequence "recovered"


def read_fasta(path):
    seqs = {}
    header = None
    buf = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header:
                    seqs[header] = "".join(buf)
                header = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.upper())
    if header:
        seqs[header] = "".join(buf)
    return seqs


def simple_identity(a, b):
    """Ungapped identity over the shorter sequence (rough but dependency-free)."""
    if not a or not b:
        return 0.0
    # align by trimming to the shorter length from both ends
    n = min(len(a), len(b))
    matches = sum(x == y for x, y in zip(a[:n], b[:n]))
    return matches / n


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--asvs",      required=True, help="Path to asvs_nonchimeras.fasta")
    parser.add_argument("--reference", required=True, help="Path to mock_sequences_V4.fasta")
    args = parser.parse_args()

    asvs = read_fasta(args.asvs)
    refs = read_fasta(args.reference)

    n_asvs = len(asvs)
    n_refs = len(refs)
    errors = []

    print(f"Reference sequences : {n_refs}")
    print(f"Recovered ASVs      : {n_asvs}")

    # Check 1: ASV count bounds
    if n_asvs < MIN_ASVS:
        errors.append(f"FAIL: Only {n_asvs} ASVs recovered, expected >= {MIN_ASVS}")
    else:
        print(f"PASS: ASV count {n_asvs} >= {MIN_ASVS}")

    if n_asvs > MAX_ASVS:
        errors.append(f"FAIL: {n_asvs} ASVs recovered, expected <= {MAX_ASVS} (too many spurious ASVs)")
    else:
        print(f"PASS: ASV count {n_asvs} <= {MAX_ASVS}")

    # Check 2: each reference sequence is recovered by at least one ASV
    missing = []
    for ref_name, ref_seq in refs.items():
        best = max(simple_identity(ref_seq, asv_seq) for asv_seq in asvs.values())
        if best < IDENTITY:
            missing.append((ref_name, best))
        else:
            print(f"PASS: {ref_name} recovered (best identity {best:.3f})")

    if missing:
        for name, best in missing:
            errors.append(f"FAIL: Reference '{name}' not recovered (best identity {best:.3f} < {IDENTITY})")

    if errors:
        print("\n=== BENCHMARK VALIDATION FAILED ===")
        for e in errors:
            print(e)
        sys.exit(1)
    else:
        print(f"\n=== BENCHMARK VALIDATION PASSED ({n_asvs} ASVs, {n_refs - len(missing)}/{n_refs} references recovered) ===")


if __name__ == "__main__":
    main()
