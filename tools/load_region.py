#!/usr/bin/env python3
"""Load one region of a shared-prep study: z aligned to a panel's a1, the matching LD submatrix, and a table.

Every analysis choice is a required argument, so nothing is decided silently:
  ambiguous    "drop" or "keep" strand-ambiguous (A/T, C/G) SNPs
  min_n_frac   keep variants with n >= min_n_frac * max(n in the region); 0 keeps all
  min_maf      keep variants with min(eaf, 1 - eaf) >= min_maf in the sumstats; 0 keeps all
Always removed, and counted in `counts`: variants with a duplicated vkey in the sumstats file, variants
without finite z, and variants absent from the panel's matrix. For other policies, load with
ambiguous="keep", min_n_frac=0, min_maf=0 and filter the returned table yourself (its `row` column
indexes the panel matrix).

Python:
  import sys; sys.path.insert(0, "tools"); from load_region import load_region
  d = load_region("prep/eadb2026", "GCST90704647", "chr19_...", "ukb_hg19", ambiguous="drop", min_n_frac=0.9, min_maf=0)
  d["z"], d["R"] (float64, m x m), d["table"] (pandas), d["counts"]
CLI (used by tools/load_region.R):
  tools/load_region.py PREP SUMSTATS REGION PANEL --ambiguous drop --min-n-frac 0.9 --min-maf 0 --out PREFIX
  writes PREFIX.tsv (the table, z in column z_aligned), PREFIX.ld.bin (float32 m x m) and PREFIX.counts.json
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

SS_COLS = ["vkey", "ea", "oa", "beta", "se", "z", "log10p", "eaf", "n", "neff", "ambiguous", "dup_vkey"]


def load_region(prep, sumstats, region, panel, *, ambiguous, min_n_frac, min_maf):
    if ambiguous not in ("drop", "keep"):
        raise ValueError('ambiguous must be "drop" or "keep"')
    ld_dir = os.path.join(prep, "ld", region)
    vars_ = pd.read_csv(os.path.join(ld_dir, f"{panel}.vars.tsv"), sep="\t", dtype={"a1": str, "a2": str})
    vars_["row"] = np.arange(len(vars_))
    chrom, start, end = region.split("_")
    ss = pd.read_parquet(os.path.join(prep, "sumstats", f"{sumstats}.parquet"), columns=SS_COLS,
                         filters=[("chr", "==", int(chrom[3:])), ("bp", ">=", int(start)), ("bp", "<=", int(end))])
    counts = {"sumstats_in_region": len(ss)}
    ss = ss[~ss["dup_vkey"]]
    counts["after_dropping_duplicated_vkeys"] = len(ss)
    ss = ss[np.isfinite(ss["z"])]
    counts["after_dropping_nonfinite_z"] = len(ss)
    t = vars_.merge(ss, on="vkey", how="inner")
    counts["in_panel"] = len(t)
    if ambiguous == "drop":
        t = t[~t["ambiguous"]]
    counts["after_ambiguous"] = len(t)
    if min_n_frac > 0:
        t = t[t["n"] >= min_n_frac * t["n"].max()]
    counts["after_min_n_frac"] = len(t)
    if min_maf > 0:
        t = t[np.minimum(t["eaf"], 1 - t["eaf"]) >= min_maf]
    counts["after_min_maf"] = len(t)
    t = t.sort_values("row").reset_index(drop=True)
    same, swapped = (t["ea"] == t["a1"]).to_numpy(), (t["ea"] == t["a2"]).to_numpy()
    assert np.all(same | swapped), "vkey matched but alleles differ; the prep files are inconsistent"
    t["sign"] = np.where(same, 1, -1)
    t["z_aligned"] = t["z"] * t["sign"]
    m_all = len(vars_)
    full = np.memmap(os.path.join(ld_dir, f"{panel}.bin"), np.float32, "r", shape=(m_all, m_all))
    rows = t["row"].to_numpy()
    R = np.asarray(full[rows][:, rows], dtype=np.float64)
    return {"z": t["z_aligned"].to_numpy(), "R": R, "table": t, "counts": counts,
            "choices": {"ambiguous": ambiguous, "min_n_frac": min_n_frac, "min_maf": min_maf}}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for a in ("prep", "sumstats", "region", "panel"):
        ap.add_argument(a)
    ap.add_argument("--ambiguous", required=True, choices=["drop", "keep"])
    ap.add_argument("--min-n-frac", required=True, type=float)
    ap.add_argument("--min-maf", required=True, type=float)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    d = load_region(a.prep, a.sumstats, a.region, a.panel, ambiguous=a.ambiguous, min_n_frac=a.min_n_frac,
                    min_maf=a.min_maf)
    d["table"].to_csv(a.out + ".tsv", sep="\t", index=False)
    d["R"].astype(np.float32).tofile(a.out + ".ld.bin")
    with open(a.out + ".counts.json", "w") as f:
        json.dump({"counts": d["counts"], "choices": d["choices"], "m": len(d["z"])}, f, indent=1)


if __name__ == "__main__":
    main()
