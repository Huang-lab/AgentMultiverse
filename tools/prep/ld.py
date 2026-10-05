#!/usr/bin/env python3
"""Precompute per-region LD matrices for every reference panel, restricted to variants in the sumstats.

Usage:
  tools/prep/ld.py 1kg PREP_DIR [REGION ...]   1000 Genomes GRCh38 EUR (525 unrelated), MAF >= 0.001
  tools/prep/ld.py ukb PREP_DIR [REGION ...]   UK Biobank GRCh37 (337K British, PolyFun), all variants in the window
  tools/prep/ld.py panukb PREP_DIR [REGION ...] Pan-UKBB EUR GRCh37 (about 420K), all variants in the region

Reads PREP_DIR/regions.tsv and PREP_DIR/sumstats/*.parquet and writes, per region and panel,
PREP_DIR/ld/<region>/<panel>.bin (float32 m x m, row-major, signed r) and <panel>.vars.tsv
(one row per matrix row: vkey, rsid, chr, bp, bp_hg19, a1, a2, a1_freq). Signs are relative to a1,
with alleles on the GRCh38 forward strand, so sign = +1 where a study's ea equals a1 and -1 otherwise.
Only panel variants whose chr, position and allele pair match a sumstats variant (in any study) are kept.
Existing outputs are skipped, so the command can be rerun to finish an interrupted build.
Appends one row per region and panel to PREP_DIR/ld/index.tsv.
"""
import fcntl
import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
from scipy import sparse

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PLINK2 = os.path.join(ROOT, "bin/plink2")
KG = os.path.join(ROOT, "ref/1kg/hg38/all_hg38")
KG_KEEP = os.path.join(ROOT, "ref/1kg/hg38/EUR.keep")
UKB_DIR = os.path.join(ROOT, "ref/ukbb_ld")
UKB_INDEX = os.path.join(ROOT, "tools/polyfun/ukb_regions.tsv.gz")
FETCH = os.path.join(ROOT, "tools/fetch.py")
PANUKB_BUCKET = "https://pan-ukb-us-east-1.s3.amazonaws.com"
PANUKB_PREFIX = "ld_release"
PANUKB_DIR = os.path.join(ROOT, "ref/panukb_ld")
PANUKB_BM = os.path.join(PANUKB_DIR, PANUKB_PREFIX, "UKBB.EUR.ldadj.bm")
PANUKB_VARS = os.path.join(PANUKB_DIR, "UKBB.EUR.variants.parquet")
PANUKB_READ = os.path.join(ROOT, "tools/prep/panukb_read.py")
S3_MIRROR = os.path.join(ROOT, "tools/prep/s3_mirror.py")
HAIL_PY = os.path.join(ROOT, "env_hail/bin/python")
COMPLEMENT = str.maketrans("ACGT", "TGCA")
VARS_COLS = ["vkey", "rsid", "chr", "bp", "bp_hg19", "a1", "a2", "a1_freq"]
INDEX_COLS = ["region", "panel", "m", "source", "source_start", "source_end", "n_sumstats", "note"]


def sorted_key(chrom, bp, a, b):
    lo, hi = np.where(a < b, a, b), np.where(a < b, b, a)
    return chrom.astype(str) + ":" + bp.astype(str) + ":" + lo + ":" + hi


def region_sumstats(prep, r):
    """Union of the region's variants across studies: vkey, rsid, bp, chr_hg19, bp_hg19, hg19_minus."""
    cols = ["vkey", "rsid", "chr", "bp", "chr_hg19", "bp_hg19", "hg19_minus"]
    parts = [pd.read_parquet(p, columns=cols, filters=[("chr", "==", r.chr), ("bp", ">=", r.start), ("bp", "<=", r.end)])
             for p in sorted(glob.glob(os.path.join(prep, "sumstats/*.parquet")))]
    d = pd.concat(parts).drop_duplicates("vkey").reset_index(drop=True)
    d["rsid"] = d["rsid"].astype(object)
    return d


def write(out, name, R, vars_):
    assert R.shape == (len(vars_), len(vars_)) and R.dtype == np.float32 and np.isfinite(R).all()
    tmp = os.path.join(out, f".{name}.bin.part")
    R.tofile(tmp)
    vars_[VARS_COLS].to_csv(os.path.join(out, f"{name}.vars.tsv"), sep="\t", index=False)
    os.replace(tmp, os.path.join(out, f"{name}.bin"))


def append_index(prep, row):
    """Both panels may be built at once, so appends are serialized with a lock."""
    path = os.path.join(prep, "ld/index.tsv")
    with open(path + ".lock", "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        new = not os.path.exists(path)
        pd.DataFrame([row], columns=INDEX_COLS).to_csv(path, sep="\t", index=False, header=new, mode="a")


# ---------------------------------------------------------------- 1000 Genomes GRCh38 EUR

def kg_chrom(work, chrom):
    """EUR-only, MAF >= 0.001, biallelic subset of one chromosome with chr:pos:ref:alt IDs (made once)."""
    prefix = os.path.join(work, f"chr{chrom}")
    if not os.path.exists(prefix + ".afreq"):
        subprocess.run([PLINK2, "--pfile", KG, "vzs", "--keep", KG_KEEP, "--chr", str(chrom), "--maf", "0.001",
                        "--max-alleles", "2", "--set-all-var-ids", "@:#:$r:$a", "--new-id-max-allele-len", "1000",
                        "missing", "--rm-dup", "exclude-all", "--make-pgen", "--freq", "--memory", "9000",
                        "--threads", "4", "--out", prefix], check=True, capture_output=True)
    return prefix


def kg_region(prep, work, r):
    prefix = kg_chrom(work, r.chr)
    with open(prefix + ".pvar") as f:
        n_meta = next(i for i, line in enumerate(f) if line.startswith("#CHROM"))
    pvar = pd.read_csv(prefix + ".pvar", sep="\t", skiprows=n_meta + 1, header=None, usecols=range(5),
                       dtype={2: str, 3: str, 4: str})
    pvar.columns = ["chr", "bp", "id", "ref", "alt"]
    pvar = pvar[(pvar["bp"] >= r.start) & (pvar["bp"] <= r.end) & (pvar["id"] != ".")]
    pvar["vkey"] = sorted_key(pvar["chr"], pvar["bp"], pvar["ref"].to_numpy(), pvar["alt"].to_numpy())
    ss = region_sumstats(prep, r)
    pvar = pvar[pvar["vkey"].isin(ss["vkey"])]
    with tempfile.TemporaryDirectory(dir=work) as tmp:
        ids = os.path.join(tmp, "ids")
        pvar["id"].to_csv(ids, index=False, header=False)
        subprocess.run([PLINK2, "--pfile", prefix, "--extract", ids, "--r-unphased", "square", "bin4", "ref-based",
                        "yes-really", "--threads", "4", "--memory", "4000", "--out", os.path.join(tmp, "ld")],
                       check=True, capture_output=True)
        order = pd.read_csv(os.path.join(tmp, "ld.unphased.vcor1.bin.vars"), header=None)[0]
        m = len(order)
        R = np.fromfile(os.path.join(tmp, "ld.unphased.vcor1.bin"), np.float32).reshape(m, m)
    freq = pd.read_csv(prefix + ".afreq", sep="\t", usecols=["ID", "ALT_FREQS"]).set_index("ID")["ALT_FREQS"]
    v = pvar.set_index("id").loc[order].reset_index()
    v = v.merge(ss[["vkey", "rsid", "bp_hg19"]], on="vkey", how="left")
    v["a1"], v["a2"], v["a1_freq"] = v["alt"], v["ref"], v["id"].map(freq).round(6).to_numpy()
    np.fill_diagonal(R, 1.0)
    return R, v, {"source": "1000G NYGC 30x EUR.keep", "source_start": r.start, "source_end": r.end,
                  "n_sumstats": len(ss), "note": ""}


# ---------------------------------------------------------------- GRCh37 panels

def match_hg19(ss, chrom, panel):
    """Match region sumstats to a GRCh37 panel with columns pos37, p1, p2 (panel alleles, GRCh37 strand).

    Sumstats alleles are complemented at minus-strand liftOver sites; minus-strand indels cannot be matched
    reliably and are left out. Keys that are not unique on either side are dropped.
    Returns (panel rows with vkey, rsid, bp, minus added; number of region sumstats variants).
    """
    n_ss = len(ss)
    ss = ss[(ss["chr_hg19"] == chrom).fillna(False).to_numpy()].reset_index(drop=True)
    alleles = ss["vkey"].str.split(":", expand=True)
    a, b = alleles[2].astype(str), alleles[3].astype(str)
    ss["minus"] = ss["hg19_minus"].fillna(False).astype(bool)
    keep = ~ss["minus"] | ((a.str.len() == 1) & (b.str.len() == 1))
    ss, a, b = ss[keep].reset_index(drop=True), a[keep].reset_index(drop=True), b[keep].reset_index(drop=True)
    a37 = np.where(ss["minus"], a.str.translate(COMPLEMENT), a)
    b37 = np.where(ss["minus"], b.str.translate(COMPLEMENT), b)
    ss["key37"] = sorted_key(ss["chr_hg19"].astype(int), ss["bp_hg19"].astype(int), a37, b37)
    panel = panel.copy()
    panel["key37"] = sorted_key(pd.Series(chrom, index=panel.index), panel["pos37"], panel["p1"].to_numpy(),
                                panel["p2"].to_numpy())
    ss = ss[~ss["key37"].duplicated(keep=False)]
    panel = panel[~panel["key37"].duplicated(keep=False)]
    v = panel.merge(ss[["key37", "vkey", "rsid", "bp", "minus"]], on="key37", suffixes=("_panel", ""))
    v = v[~v["vkey"].duplicated(keep=False)]
    return v, n_ss


def orient(v, chrom):
    """a1/a2 on the GRCh38 forward strand from GRCh37 panel alleles p1 (counted) and p2."""
    comp = lambda s: s.str.translate(COMPLEMENT)
    v["a1"] = np.where(v["minus"], comp(v["p1"]), v["p1"])
    v["a2"] = np.where(v["minus"], comp(v["p2"]), v["p2"])
    v["chr"], v["bp_hg19"] = chrom, v["pos37"]
    return v


# ---------------------------------------------------------------- UK Biobank GRCh37 (PolyFun)

def ukb_window(r):
    w = pd.read_csv(UKB_INDEX, sep="\t")
    w = w[w["CHR"] == r.chr].copy()
    w["overlap"] = np.minimum(w["END"], r.end_hg19) - np.maximum(w["START"], r.start_hg19)
    w["offcentre"] = abs((w["START"] + w["END"]) / 2 - (r.start_hg19 + r.end_hg19) / 2)
    w = w[w["overlap"] > 0].sort_values(["overlap", "offcentre"], ascending=[False, True])
    return None if w.empty else w.iloc[0]


def ukb_region(prep, work, r):
    w = ukb_window(r)
    ss = region_sumstats(prep, r)
    if w is None:
        return None, None, {"source": "none", "source_start": "", "source_end": "", "n_sumstats": len(ss),
                            "note": "no UK Biobank window overlaps this region"}
    name = os.path.basename(w["URL_PREFIX"])
    base = os.path.join(UKB_DIR, name)
    for ext in (".gz", ".npz"):
        subprocess.run([sys.executable, FETCH, w["URL_PREFIX"] + ext, base + ext], check=True)
    snps = pd.read_csv(base + ".gz", sep="\t", dtype={"allele1": str, "allele2": str})
    snps["row"] = np.arange(len(snps))

    v, n_ss = match_hg19(ss, r.chr, snps.rename(columns={"position": "pos37", "allele1": "p1", "allele2": "p2"}))
    v = v.rename(columns={"p1": "allele1", "p2": "allele2", "pos37": "position"}).sort_values("row").reset_index(drop=True)

    z = np.load(base + ".npz")
    L = sparse.coo_matrix((z["data"], (z["row"], z["col"])), shape=tuple(z["shape"])).tocsr()
    del z
    idx = v["row"].to_numpy()
    sub = L[idx][:, idx].toarray()
    del L
    R = (sub + sub.T).astype(np.float32)
    assert np.allclose(np.diag(R), 1.0, atol=1e-4), "UK Biobank diagonal is not 1"
    np.fill_diagonal(R, 1.0)
    v = orient(v.rename(columns={"allele1": "p1", "allele2": "p2", "position": "pos37"}), r.chr)
    v["a1_freq"] = np.nan
    covered = (max(w["START"], r.start_hg19), min(w["END"], r.end_hg19))
    note = "" if covered == (r.start_hg19, r.end_hg19) else \
        f"window covers GRCh37 {covered[0]}-{covered[1]} of {r.start_hg19}-{r.end_hg19}"
    return R, v, {"source": name, "source_start": int(w["START"]), "source_end": int(w["END"]),
                  "n_sumstats": n_ss, "note": note}


# ---------------------------------------------------------------- Pan-UKBB EUR GRCh37

def panukb_region(prep, work, r):
    """Download only the BlockMatrix blocks the region touches, read them with Hail, then delete them."""
    ss = region_sumstats(prep, r)
    panel = pd.read_parquet(PANUKB_VARS, filters=[("chr_hg19", "==", str(r.chr)), ("bp_hg19", ">=", int(r.start_hg19)),
                                                   ("bp_hg19", "<=", int(r.end_hg19))])
    panel = panel.drop(columns="rsid").rename(columns={"bp_hg19": "pos37", "alt": "p1", "ref": "p2"})
    v, n_ss = match_hg19(ss, r.chr, panel)
    v = orient(v.sort_values("idx").reset_index(drop=True), r.chr)
    v["a1_freq"] = v["alt_freq"].round(6)
    idx = v["idx"].to_numpy()

    with open(os.path.join(PANUKB_BM, "metadata.json")) as f:
        meta = json.load(f)
    size, nbr = meta["blockSize"], -(-meta["nRows"] // meta["blockSize"])
    parts = dict(zip(meta["maybeFiltered"], meta["partFiles"]))
    blocks = np.unique(idx // size)
    keys = [f"{PANUKB_PREFIX}/UKBB.EUR.ldadj.bm/parts/{parts[i + j * nbr]}"
            for i in blocks for j in blocks if i <= j and i + j * nbr in parts]
    missing = [(i, j) for i in blocks for j in blocks if i <= j and i + j * nbr not in parts]
    with tempfile.TemporaryDirectory(dir=work) as tmp:
        with open(os.path.join(tmp, "keys"), "w") as f:
            f.write("\n".join(keys))
        subprocess.run([sys.executable, S3_MIRROR, PANUKB_BUCKET, "x", PANUKB_DIR, "--keys", os.path.join(tmp, "keys"),
                        "--jobs", "8"], check=True)
        np.save(os.path.join(tmp, "idx.npy"), idx)
        out = os.path.join(tmp, "R.bin")
        subprocess.run([HAIL_PY, PANUKB_READ, PANUKB_BM, os.path.join(tmp, "idx.npy"), out], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        R = np.fromfile(out, np.float32).reshape(len(idx), len(idx))
    for k in keys:
        os.remove(os.path.join(PANUKB_DIR, k))
    note = f"{len(missing)} block pairs outside the stored band are 0" if missing else ""
    return R, v, {"source": "Pan-UKBB UKBB.EUR.ldadj.bm", "source_start": int(r.start_hg19),
                  "source_end": int(r.end_hg19), "n_sumstats": n_ss, "note": note}


PANELS = {"1kg": ("1kg_hg38_EUR", kg_region), "ukb": ("ukb_hg19", ukb_region),
          "panukb": ("panukb_hg19_EUR", panukb_region)}


def main():
    if len(sys.argv) < 3 or sys.argv[1] not in PANELS:
        sys.exit(__doc__)
    panel, prep, only = sys.argv[1], sys.argv[2], set(sys.argv[3:])
    name, build = PANELS[panel]
    regions = pd.read_csv(os.path.join(prep, "regions.tsv"), sep="\t")
    work = os.path.join(prep, ".work", panel)
    os.makedirs(work, exist_ok=True)
    for r in regions.itertuples():
        if only and r.region not in only:
            continue
        out = os.path.join(prep, "ld", r.region)
        if os.path.exists(os.path.join(out, f"{name}.bin")):
            continue
        os.makedirs(out, exist_ok=True)
        R, v, info = build(prep, work, r)
        if R is not None:
            write(out, name, R, v)
        append_index(prep, {"region": r.region, "panel": name, "m": 0 if R is None else len(v), **info})
        print(f"{r.region} {name}: m = {0 if R is None else len(v):,} of {info['n_sumstats']:,} sumstats variants"
              + (f" ({info['note']})" if info["note"] else ""), flush=True)


if __name__ == "__main__":
    main()
