#!/usr/bin/env python3
"""Define shared LD regions and per-study leads from harmonized sumstats.

Usage: tools/prep/regions.py PREP_DIR [FLANK_BP]

Every variant with P < 5e-8 (from z) in any study, extended by FLANK_BP (default 500000) on each
side, and overlapping intervals merged, gives one region, named chr<c>_<start>_<end> in GRCh38.
Writes PREP_DIR/regions.tsv (one row per region) and PREP_DIR/leads.tsv (one row per study and
region: the study's top variant in the region and how many variants reach P < 5e-8 there).
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

SIG = -np.log10(5e-8)
MHC_HG38 = (6, 25_000_000, 34_000_000)
COLS = ["vkey", "chr", "bp", "rsid", "ea", "oa", "z", "log10p", "eaf", "n", "chr_hg19", "bp_hg19"]


def merge(intervals):
    out = []
    for c, s, e in sorted(intervals):
        if out and out[-1][0] == c and s <= out[-1][2]:
            out[-1][2] = max(out[-1][2], e)
        else:
            out.append([c, s, e])
    return out


def main():
    prep = sys.argv[1]
    flank = int(sys.argv[2]) if len(sys.argv) > 2 else 500_000
    studies = {os.path.basename(p)[:-8]: p for p in sorted(glob.glob(os.path.join(prep, "sumstats/*.parquet")))}
    sig = {s: pd.read_parquet(p, columns=COLS, filters=[("log10p", ">", SIG)]) for s, p in studies.items()}
    regions = merge([(c, max(1, b - flank), b + flank) for d in sig.values() for c, b in zip(d["chr"], d["bp"])])
    regions = pd.DataFrame(regions, columns=["chr", "start", "end"])
    regions.insert(0, "region", "chr" + regions["chr"].astype(str) + "_" + regions["start"].astype(str)
                   + "_" + regions["end"].astype(str))
    c, s, e = MHC_HG38
    regions["mhc"] = (regions["chr"] == c) & (regions["start"] < e) & (regions["end"] > s)

    leads, spans = [], []
    for study, path in studies.items():
        for r in regions.itertuples():
            d = pd.read_parquet(path, columns=COLS, filters=[("chr", "==", r.chr), ("bp", ">=", r.start),
                                                             ("bp", "<=", r.end)])
            spans.append(d[["chr_hg19", "bp_hg19"]].dropna().assign(region=r.region))
            top = d.loc[d["log10p"].idxmax()]
            leads.append({"study": study, "region": r.region, "n_variants": len(d),
                          "n_sig": int((d["log10p"] > SIG).sum()), **top[COLS].to_dict()})
    leads = pd.DataFrame(leads)
    h19 = pd.concat(spans)
    h19 = h19[h19["chr_hg19"] == h19["region"].map(regions.set_index("region")["chr"])]
    h19 = h19.groupby("region")["bp_hg19"].agg(start_hg19="min", end_hg19="max")
    regions = regions.join(h19, on="region")
    best = leads.sort_values("log10p", ascending=False).drop_duplicates("region").set_index("region")
    regions["max_log10p"] = regions["region"].map(best["log10p"])
    regions["top_study"] = regions["region"].map(best["study"])
    regions["top_rsid"] = regions["region"].map(best["rsid"])
    regions.to_csv(os.path.join(prep, "regions.tsv"), sep="\t", index=False)
    leads.to_csv(os.path.join(prep, "leads.tsv"), sep="\t", index=False)
    widths = (regions["end"] - regions["start"]) / 1e6
    print(f"{len(regions)} regions, width median {widths.median():.2f} Mb, max {widths.max():.2f} Mb, "
          f"total {widths.sum():.0f} Mb; per study with P < 5e-8: "
          + ", ".join(f"{s} {int((g['n_sig'] > 0).sum())}" for s, g in leads.groupby("study")))


if __name__ == "__main__":
    main()
