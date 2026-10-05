#!/usr/bin/env python3
"""Read a symmetric submatrix of the Pan-UKBB LD BlockMatrix (run with env_hail/bin/python).

Usage: panukb_read.py BM_DIR IDX.npy OUT.bin
  IDX.npy holds sorted matrix indices. Only the blocks those indices touch need to be present
  under BM_DIR/parts; the release stores the upper triangle, which is mirrored here.
  Writes float32 m x m, row-major, with the diagonal set to 1.
"""
import os
import sys

import numpy as np

os.environ.setdefault("JAVA_HOME", os.path.join(os.path.dirname(os.path.dirname(sys.executable)), "lib/jvm"))
import hail as hl  # noqa: E402


def main():
    bm_dir, idx_file, out = sys.argv[1:4]
    idx = np.load(idx_file).astype(np.int64)
    assert np.all(np.diff(idx) > 0), "indices must be sorted and unique"
    m = len(idx)
    driver_gb = max(8, int(np.ceil(m * m * 8 * 2.5 / 1e9)) + 2)
    hl.init(quiet=True, master="local[4]", spark_conf={"spark.driver.memory": f"{driver_gb}g"}, log="/dev/null")
    U = hl.linalg.BlockMatrix.read(bm_dir).filter(idx.tolist(), idx.tolist()).to_numpy()
    hl.stop()
    R = U.astype(np.float32)  # upper triangle only: absent blocks and the lower half of diagonal blocks are 0
    del U
    step = 2048
    for s in range(step, m, step):  # mirror the upper triangle into the lower one, a column band at a time
        R[s:s + step, :s] = R[:s, s:s + step].T
    for s in range(0, m, step):
        blk = R[s:s + step, s:s + step]
        blk += np.triu(blk, 1).T
    np.fill_diagonal(R, 1.0)
    tmp = out + ".part"
    R.tofile(tmp)
    os.replace(tmp, out)


if __name__ == "__main__":
    main()
