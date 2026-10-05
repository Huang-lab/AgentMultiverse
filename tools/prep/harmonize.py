#!/usr/bin/env python3
"""Harmonize GWAS-SSF summary statistics (GRCh38) into one parquet per study, without dropping variants.

Usage: tools/prep/harmonize.py OUT_DIR SUMSTATS.tsv.gz [...]

Every input row is kept. Instead of filtering, rows get flags that an analysis can act on:
  vkey        chr:bp:A:B with the two alleles in sorted order, the same key the LD variant tables use
  ambiguous   A/T or C/G SNP
  dup_vkey    the vkey occurs more than once in this study
z is beta/se, and log10p is computed from z, so neither underflows at strong signals.
Positions are lifted to GRCh37 with UCSC liftOver; hg19_minus marks sites that map to the minus
strand, whose alleles must be complemented to match GRCh37 references.
"""
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.csv as pcsv
import pyarrow.parquet as pq
from scipy.stats import norm

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CHAIN = os.path.join(ROOT, "ref/liftover/hg38ToHg19.over.chain.gz")
LIFTOVER = os.path.join(ROOT, "bin/liftOver")

TYPES = {
    "chromosome": pa.string(), "base_pair_location": pa.int64(), "effect_allele": pa.string(),
    "other_allele": pa.string(), "beta": pa.float64(), "standard_error": pa.float64(),
    "effect_allele_frequency": pa.float64(), "p_value": pa.float64(), "variant_id": pa.string(),
    "rs_id": pa.string(), "N_cases": pa.float64(), "N_controls": pa.float64(), "Neff_cases": pa.float64(),
    "Neff_controls": pa.float64(), "Neff_total": pa.float64(), "HetISq": pa.float64(), "HetPVal": pa.float64(),
}
RENAME = {
    "base_pair_location": "bp", "effect_allele": "ea", "other_allele": "oa", "standard_error": "se",
    "effect_allele_frequency": "eaf", "p_value": "p", "rs_id": "rsid", "N_cases": "n_cases",
    "N_controls": "n_controls", "Neff_cases": "neff_cases", "Neff_controls": "neff_controls",
    "Neff_total": "neff", "HetISq": "het_isq", "HetPVal": "het_p",
}


def read_ssf(path):
    t = pcsv.read_csv(path, parse_options=pcsv.ParseOptions(delimiter="\t"),
                      convert_options=pcsv.ConvertOptions(column_types=TYPES, null_values=["NA", ""]))
    df = t.to_pandas(types_mapper={pa.string(): pd.ArrowDtype(pa.string())}.get).rename(columns=RENAME)
    df = df.drop(columns=["variant_id"])
    df.insert(0, "chr", df.pop("chromosome").astype(str).str.removeprefix("chr").astype(np.int8))
    return df


def add_columns(df):
    ea, oa = df["ea"].astype(str).str.upper(), df["oa"].astype(str).str.upper()
    df["ea"], df["oa"] = ea, oa
    lo, hi = np.where(ea < oa, ea, oa), np.where(ea < oa, oa, ea)
    df["vkey"] = df["chr"].astype(str) + ":" + df["bp"].astype(str) + ":" + lo + ":" + hi
    df["z"] = df["beta"] / df["se"]
    df["log10p"] = -(np.log(2) + norm.logsf(np.abs(df["z"]))) / np.log(10)
    df["n"] = df["n_cases"] + df["n_controls"]
    pair = ea + oa
    df["ambiguous"] = pair.isin(["AT", "TA", "CG", "GC"])
    df["dup_vkey"] = df["vkey"].duplicated(keep=False)
    return df


def liftover(sites):
    """sites: DataFrame chr, bp (GRCh38). Returns it with chr_hg19, bp_hg19, hg19_minus (NA if unmapped)."""
    sites = sites.drop_duplicates().reset_index(drop=True)
    with tempfile.TemporaryDirectory() as tmp:
        bed, out, unm = (os.path.join(tmp, f) for f in ("in.bed", "out.bed", "unmapped.bed"))
        pd.DataFrame({"c": "chr" + sites["chr"].astype(str), "s": sites["bp"] - 1, "e": sites["bp"],
                      "i": np.arange(len(sites)), "sc": 0, "st": "+"}).to_csv(bed, sep="\t", header=False, index=False)
        subprocess.run([LIFTOVER, bed, CHAIN, out, unm], check=True, capture_output=True)
        lifted = pd.read_csv(out, sep="\t", header=None, names=["c", "s", "e", "i", "sc", "st"])
    lifted = lifted[lifted["c"].str.fullmatch(r"chr\d+")]
    sites["chr_hg19"] = pd.array([pd.NA] * len(sites), dtype="Int8")
    sites["bp_hg19"] = pd.array([pd.NA] * len(sites), dtype="Int64")
    sites["hg19_minus"] = pd.array([pd.NA] * len(sites), dtype="boolean")
    idx = lifted["i"].to_numpy()
    sites.loc[idx, "chr_hg19"] = lifted["c"].str[3:].astype(int).to_numpy()
    sites.loc[idx, "bp_hg19"] = lifted["e"].to_numpy()
    sites.loc[idx, "hg19_minus"] = (lifted["st"] == "-").to_numpy()
    return sites


def main():
    out_dir, paths = sys.argv[1], sys.argv[2:]
    if not paths:
        sys.exit(__doc__)
    os.makedirs(out_dir, exist_ok=True)
    for path in paths:
        study = os.path.basename(path).split(".")[0]
        df = add_columns(read_ssf(path))
        df = df.merge(liftover(df[["chr", "bp"]]), on=["chr", "bp"], how="left")
        df = df.sort_values(["chr", "bp", "vkey"], kind="stable").reset_index(drop=True)
        order = ["vkey", "chr", "bp", "ea", "oa", "rsid", "beta", "se", "z", "p", "log10p", "eaf", "n", "n_cases",
                 "n_controls", "neff", "neff_cases", "neff_controls", "het_isq", "het_p", "ambiguous", "dup_vkey",
                 "chr_hg19", "bp_hg19", "hg19_minus"]
        table = pa.Table.from_pandas(df[order], preserve_index=False)
        pq.write_table(table, os.path.join(out_dir, f"{study}.parquet"), row_group_size=500_000, compression="zstd")
        print(f"{study}: {len(df):,} rows, {int(df['dup_vkey'].sum()):,} in duplicated vkeys, "
              f"{int(df['ambiguous'].sum()):,} ambiguous, {int(df['bp_hg19'].isna().sum()):,} not lifted, "
              f"{int((df['log10p'] > -np.log10(5e-8)).sum()):,} with P < 5e-8", flush=True)


if __name__ == "__main__":
    main()
