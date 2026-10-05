#!/usr/bin/env python3
"""Compare the universes of a cohort: spread, agreement, attribution to pinned axes, and score against truth.

Usage: compare.py COHORT_DIR [--keys shared|union] [--table pip.tsv] [--key study,region,vkey] [--value pip] [--truth truth.tsv] [--top 20]

Every subfolder of COHORT_DIR holding TABLE is a universe (spec.json, if present, gives its pinned choices).
The defaults read the fine-mapping PIP table; any table with a key and a numeric value works,
and when a universe table has several `config` rows the value is averaged over them (use --by-config to keep them apart).

--keys shared (default) compares only keys that every universe reports, so a difference in scope is not mistaken for disagreement;
--keys union keeps all keys and counts a key a universe does not report as 0.
Writes into COHORT_DIR:
  coverage.tsv     per universe: keys reported, keys shared by all, keys only this universe reports
  spread.tsv       per key: universes covering it, mean, sd, min, max of the value
  agreement.tsv    pairwise Pearson correlation of the value across universes
  attribution.tsv  per pinned axis: share of between-universe variance explained (eta^2); free-design cohorts have none
  score.tsv        with --truth (columns = the key columns, one row per true item): per universe the value mass on true items
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("cohort_dir")
    p.add_argument("--keys", choices=("shared", "union"), default="shared")
    p.add_argument("--table", default="pip.tsv")
    p.add_argument("--key", default="study,region,vkey")
    p.add_argument("--value", default="pip")
    p.add_argument("--truth")
    p.add_argument("--top", type=int, default=20)
    a = p.parse_args()
    key = a.key.split(",")

    vecs, choices = {}, {}
    for path in sorted(glob.glob(os.path.join(a.cohort_dir, "*", a.table))):
        uid = os.path.basename(os.path.dirname(path))
        df = pd.read_csv(path, sep="\t", usecols=key + [a.value])
        vecs[uid] = df.groupby(key)[a.value].mean()
        sp = os.path.join(a.cohort_dir, uid, "spec.json")
        choices[uid] = json.load(open(sp)).get("choices", {}) if os.path.exists(sp) else {}
    if len(vecs) < 2:
        raise SystemExit(f"need at least 2 universes with {a.table} under {a.cohort_dir}, found {len(vecs)}")

    M = pd.DataFrame(vecs)                      # rows = keys, columns = universes
    covered = M.notna().sum(axis=1)
    shared = covered == M.shape[1]
    pd.DataFrame({"keys": M.notna().sum(), "shared_by_all": int(shared.sum()),
                  "only_this": M.notna().mul(covered == 1, axis=0).sum()}).rename_axis("universe").to_csv(
        os.path.join(a.cohort_dir, "coverage.tsv"), sep="\t")
    if a.keys == "shared":
        if not shared.any():
            raise SystemExit("no key is reported by every universe; use --keys union or align the scope")
        M, covered = M[shared], covered[shared]
    F = M.fillna(0.0)                           # union mode: a key a universe never reports carries 0
    spread = pd.DataFrame({"n_universes": covered, "mean": F.mean(axis=1), "sd": F.std(axis=1, ddof=0),
                           "min": F.min(axis=1), "max": F.max(axis=1)}).sort_values("mean", ascending=False)
    spread.to_csv(os.path.join(a.cohort_dir, "spread.tsv"), sep="\t")
    F.corr().to_csv(os.path.join(a.cohort_dir, "agreement.tsv"), sep="\t")

    rows = []
    total = ((F.sub(F.mean(axis=1), axis=0)) ** 2).to_numpy().sum()
    axes = sorted({x for c in choices.values() for x in c})
    for ax in axes:
        lev = pd.Series({u: choices[u].get(ax) for u in F.columns})
        ss = 0.0
        for _, us in lev.groupby(lev):
            g = F[us.index]
            ss += len(us) * ((g.mean(axis=1) - F.mean(axis=1)) ** 2).sum()
        rows.append({"axis": ax, "levels": lev.nunique(), "eta2": ss / total if total else np.nan})
    pd.DataFrame(rows, columns=["axis", "levels", "eta2"]).to_csv(
        os.path.join(a.cohort_dir, "attribution.tsv"), sep="\t", index=False)

    if a.truth:
        t = pd.read_csv(a.truth, sep="\t", usecols=key).drop_duplicates()
        idx = pd.MultiIndex.from_frame(t)
        sc = []
        for u in F.columns:
            v = F[u].reindex(idx).fillna(0.0)
            sc.append({"universe": u, "n_truth": len(idx), "mass_on_truth": v.sum(),
                       "mean_value_on_truth": v.mean(), "n_truth_missed": int((v < 0.5).sum())})
        pd.DataFrame(sc).to_csv(os.path.join(a.cohort_dir, "score.tsv"), sep="\t", index=False)

    print(f"{len(F.columns)} universes, {len(F)} keys ({a.keys}); mean pairwise r = "
          f"{(F.corr().to_numpy().sum() - len(F.columns)) / (len(F.columns) ** 2 - len(F.columns)):.3f}")
    print("most contested keys (largest sd):")
    print(spread.sort_values("sd", ascending=False).head(a.top).round(3).to_string())
    if rows:
        print("variance explained by pinned axes (eta^2):")
        print(pd.DataFrame(rows).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
