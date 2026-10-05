#!/usr/bin/env python3
"""Gene table from the GENCODE v50 basic GTFs: one row per gene, for GRCh38 and GRCh37.

Usage: tools/prep/genes.py OUT.tsv.gz
Columns: gene_id, gene_name, gene_type, chr, start, end, strand (GRCh38), then start_hg19, end_hg19, strand_hg19
(from the lift37 GTF, empty where GENCODE did not map the gene). Coordinates are 1-based and inclusive.
"""
import os
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
GTF38 = os.path.join(ROOT, "ref/gencode/gencode.v50.basic.annotation.gtf.gz")
GTF37 = os.path.join(ROOT, "ref/gencode/gencode.v50lift37.basic.annotation.gtf.gz")


def genes(path):
    g = pd.read_csv(path, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4, 6, 8],
                    names=["chr", "feature", "start", "end", "strand", "attr"])
    g = g[g["feature"] == "gene"].copy()
    for key in ("gene_id", "gene_name", "gene_type"):
        g[key] = g["attr"].str.extract(rf'{key} "([^"]+)"')[0]
    g["gene_id"] = g["gene_id"].str.replace(r"_\d+$", "", regex=True)  # lift37 appends a remap suffix
    g["chr"] = g["chr"].str.removeprefix("chr")
    return g[["gene_id", "gene_name", "gene_type", "chr", "start", "end", "strand"]]


def main():
    g38 = genes(GTF38)
    g37 = genes(GTF37).drop_duplicates("gene_id")
    g37 = g37[["gene_id", "chr", "start", "end", "strand"]].rename(
        columns={"chr": "chr_hg19", "start": "start_hg19", "end": "end_hg19", "strand": "strand_hg19"})
    out = g38.merge(g37, on="gene_id", how="left")
    out = out[out["chr_hg19"].isna() | (out["chr_hg19"] == out["chr"])].drop(columns="chr_hg19")
    for c in ("start_hg19", "end_hg19"):
        out[c] = out[c].astype("Int64")
    out.to_csv(sys.argv[1], sep="\t", index=False)
    print(f"{len(out):,} genes, {int(out['start_hg19'].notna().sum()):,} with GRCh37 coordinates")


if __name__ == "__main__":
    main()
