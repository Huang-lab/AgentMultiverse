#!/usr/bin/env python3
"""Join PolyFun's precomputed functional priors (snpvar_bin, GRCh37 UK Biobank SNPs) to a study's vkeys.

Usage: tools/prep/priors.py PREP_DIR
Writes PREP_DIR/priors/polyfun_snpvar_bin.parquet with vkey, polyfun_snp, snpvar_bin, for every
sumstats variant (any file of the study) whose GRCh37 position and allele pair match a PolyFun SNP.
Values are PolyFun's as published: no flooring or normalization.
"""
import glob
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ld import ROOT, match_hg19  # noqa: E402

POLYFUN = sorted(glob.glob(os.path.join(ROOT, "tools/polyfun/snpvar_meta.chr*.parquet")))


def main():
    prep = sys.argv[1]
    studies = sorted(glob.glob(os.path.join(prep, "sumstats/*.parquet")))
    cols = ["vkey", "rsid", "chr", "bp", "chr_hg19", "bp_hg19", "hg19_minus"]
    pf = pd.concat(pd.read_parquet(p, columns=["CHR", "BP", "SNP", "A1", "A2", "snpvar_bin"]) for p in POLYFUN)
    out = []
    for chrom in range(1, 23):
        ss = pd.concat(pd.read_parquet(p, columns=cols, filters=[("chr_hg19", "==", chrom)]) for p in studies)
        ss = ss.drop_duplicates("vkey").reset_index(drop=True)
        ss["rsid"] = ss["rsid"].astype(object)
        panel = pf[pf["CHR"] == chrom].rename(columns={"BP": "pos37", "A1": "p1", "A2": "p2", "SNP": "polyfun_snp"})
        v, n_ss = match_hg19(ss, chrom, panel)
        out.append(v[["vkey", "polyfun_snp", "snpvar_bin"]])
        print(f"chr{chrom}: {len(v):,} of {n_ss:,} variants matched", flush=True)
    out = pd.concat(out, ignore_index=True)
    os.makedirs(os.path.join(prep, "priors"), exist_ok=True)
    out.to_parquet(os.path.join(prep, "priors/polyfun_snpvar_bin.parquet"), index=False)
    print(f"{len(out):,} variants with PolyFun snpvar_bin")


if __name__ == "__main__":
    main()
