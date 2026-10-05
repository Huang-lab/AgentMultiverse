#!/usr/bin/env python3
"""Export the Pan-UKBB EUR LD variant index to parquet (run with env_hail/bin/python from the project root).

Reads ref/panukb_ld/ld_release/UKBB.EUR.ldadj.variant.ht (mirrored with tools/prep/s3_mirror.py) and writes
ref/panukb_ld/UKBB.EUR.variants.parquet: idx (matrix row), chr_hg19, bp_hg19, ref, alt, alt_freq, rsid,
sorted by idx and checked to cover every matrix row.
"""
import os
import sys

import pandas as pd

os.environ.setdefault("JAVA_HOME", os.path.join(os.path.dirname(os.path.dirname(sys.executable)), "lib/jvm"))
import hail as hl  # noqa: E402

SRC = "ref/panukb_ld/ld_release/UKBB.EUR.ldadj.variant.ht"
TSV = "ref/panukb_ld/UKBB.EUR.variants.tsv.bgz"
OUT = "ref/panukb_ld/UKBB.EUR.variants.parquet"


def main():
    hl.init(quiet=True, master="local[4]", spark_conf={"spark.driver.memory": "8g"}, log="/dev/null")
    ht = hl.read_table(SRC).key_by()
    ht = ht.select(idx=ht.idx, chr_hg19=ht.locus.contig, bp_hg19=ht.locus.position, ref=ht.alleles[0],
                   alt=ht.alleles[1], alt_freq=ht.AF, rsid=ht.rsid)
    ht.export(TSV)
    hl.stop()
    d = pd.read_csv(TSV, sep="\t", compression="gzip", dtype={"chr_hg19": str, "ref": str, "alt": str})
    d = d.sort_values("idx").reset_index(drop=True)
    assert (d["idx"].to_numpy() == range(len(d))).all(), "variant index does not cover every matrix row"
    d.to_parquet(OUT, index=False)
    for f in (TSV, os.path.join(os.path.dirname(TSV), "." + os.path.basename(TSV) + ".crc")):
        if os.path.exists(f):
            os.remove(f)
    print(f"{len(d):,} variants -> {OUT}")


if __name__ == "__main__":
    main()
