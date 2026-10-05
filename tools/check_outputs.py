#!/usr/bin/env python3
"""Check that a run folder holds the standard result tables, so runs can be compared by joining them.

Usage: tools/check_outputs.py RUN_DIR [--prep PREP_DIR]

Required files (tab-separated, with a header):
  configs.tsv  config, description
      One row per combination of choices you ran; `config` is your short label, `description` says
      in words what it is (panel, filters, method, prior, ...). Record the reasons in the run log.
  pip.tsv      config, study, region, method, vkey, rsid, pip, cs
      One row per variant fine-mapped; `cs` is the credible set it belongs to (1, 2, ...) or 0.
  cs.tsv       config, study, region, method, cs, size, coverage, min_abs_r, top_vkey, top_pip, vkeys
      One row per credible set; `vkeys` is a comma-separated list of its members.
`study` is the sumstats file (e.g. GCST90704647), `region` is a name from PREP_DIR/regions.tsv when you
use shared-prep regions (otherwise chr<c>_<start>_<end> in GRCh38), and `vkey` is chr:bp:A:B (GRCh38,
alleles sorted). With --prep, regions and vkeys are also checked against the shared prep.
Exits 1 and lists the problems if anything is off.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

SCHEMA = {
    "configs.tsv": ["config", "description"],
    "pip.tsv": ["config", "study", "region", "method", "vkey", "rsid", "pip", "cs"],
    "cs.tsv": ["config", "study", "region", "method", "cs", "size", "coverage", "min_abs_r", "top_vkey", "top_pip",
               "vkeys"],
}
VKEY = r"^\d+:\d+:[ACGTN]+:[ACGTN]+$"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--prep")
    a = ap.parse_args()
    problems, t = [], {}
    for name, cols in SCHEMA.items():
        path = os.path.join(a.run_dir, name)
        if not os.path.exists(path):
            problems.append(f"{name} is missing")
            continue
        df = pd.read_csv(path, sep="\t", dtype={"config": str, "study": str, "region": str, "vkey": str})
        missing = [c for c in cols if c not in df.columns]
        if missing:
            problems.append(f"{name} lacks columns {missing}")
        t[name] = df
    if "pip.tsv" in t and not problems:
        p = t["pip.tsv"]
        if not p["pip"].between(0, 1).all():
            problems.append("pip.tsv: pip outside [0, 1]")
        if (p["cs"] < 0).any() or not np.issubdtype(p["cs"].dtype, np.integer):
            problems.append("pip.tsv: cs must be a non-negative integer")
        bad = ~p["vkey"].str.match(VKEY)
        if bad.any():
            problems.append(f"pip.tsv: {int(bad.sum())} vkeys not of the form chr:bp:A:B, e.g. {p['vkey'][bad].iloc[0]}")
        dup = p.duplicated(["config", "study", "region", "method", "vkey"])
        if dup.any():
            problems.append(f"pip.tsv: {int(dup.sum())} duplicated (config, study, region, method, vkey) rows")
    if "cs.tsv" in t and "pip.tsv" in t and not problems:
        c, p = t["cs.tsv"], t["pip.tsv"]
        members = p[p["cs"] > 0].groupby(["config", "study", "region", "method", "cs"]).size()
        sizes = c.set_index(["config", "study", "region", "method", "cs"])["size"]
        if not members.sort_index().equals(sizes.sort_index().rename(None)):
            problems.append("cs.tsv sizes do not match the cs labels in pip.tsv")
        if (c["vkeys"].str.count(",") + 1 != c["size"]).any():
            problems.append("cs.tsv: vkeys list length differs from size")
    if "configs.tsv" in t:
        known = set(t["configs.tsv"]["config"])
        for name in ("pip.tsv", "cs.tsv"):
            if name in t and "config" in t[name]:
                extra = set(t[name]["config"]) - known
                if extra:
                    problems.append(f"{name}: configs not described in configs.tsv: {sorted(extra)[:5]}")
    if a.prep and "pip.tsv" in t and not problems:
        regions = set(pd.read_csv(os.path.join(a.prep, "regions.tsv"), sep="\t")["region"])
        used = set(t["pip.tsv"]["region"])
        if used - regions:
            print(f"note: regions not in {a.prep}/regions.tsv (fine if you defined your own): {sorted(used - regions)[:5]}")
    if problems:
        print("\n".join(f"problem: {x}" for x in problems))
        sys.exit(1)
    p = t["pip.tsv"]
    print(f"ok: {p['config'].nunique()} configs, {p['region'].nunique()} regions, {len(p):,} PIP rows, "
          f"{len(t['cs.tsv']):,} credible sets")


if __name__ == "__main__":
    main()
