#!/usr/bin/env python3
"""Shared result cache for slow fine-mapping runs (Python side; the R side is tools/fmcache.R).

FINEMAP through the cache:
  tools/fmcache.py finemap --z X.z --ld X.ld --n-samples N --out PREFIX [--replicate K] [-- FINEMAP ARGS ...]
    Runs `bin/finemap --sss --in-files <master> FINEMAP ARGS` on copies of X.z and X.ld, and writes
    PREFIX.snp, PREFIX.config, PREFIX.cred*, PREFIX.log_sss. The key is the md5 of both files' bytes,
    N, the FINEMAP arguments, the FINEMAP version and K. FINEMAP 1.4.2 cannot be seeded, so a hit
    returns the exact earlier run; use --replicate 2, 3, ... for independent runs.

From Python:
  from fmcache import cached_call
  res = cached_call("mymodule.fn", fn, {"z": z, "R": R, "n": n}, version="1.0")
    Calls fn(**args) once per key; arrays are keyed by the md5 of their float64 bytes (column-major).

Both print whether the result was a hit; record that with tools/runlog.py.
"""
import argparse
import glob
import hashlib
import json
import os
import pickle
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, ".fmcache")
FINEMAP = os.path.join(ROOT, "bin/finemap")


def md5_file(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


def spec(x):
    if isinstance(x, dict):
        return {k: spec(x[k]) for k in sorted(x)}
    if isinstance(x, (list, tuple)):
        return [spec(v) for v in x]
    if isinstance(x, np.ndarray) and x.size > 1:
        a = np.asarray(x, dtype="<f8")
        h = hashlib.md5()
        flat = a.reshape(-1, order="F")
        for s in range(0, flat.size, 1 << 24):
            h.update(np.ascontiguousarray(flat[s:s + (1 << 24)]).tobytes())
        return {"md5": h.hexdigest(), "length": int(a.size), "dim": list(a.shape) if a.ndim > 1 else None}
    if isinstance(x, np.generic):
        return x.item()
    if callable(x):
        raise TypeError("fmcache: functions cannot be part of a cache key")
    return x


def alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def run_cached(name, key_spec, compute):
    """compute(workdir) writes the result into workdir; returns (dir holding the result, hit, seconds)."""
    text = json.dumps(key_spec, sort_keys=True, separators=(",", ":"))
    key = hashlib.md5(text.encode()).hexdigest()
    d = os.path.join(CACHE, name.replace("/", "_").replace(":", "_"), key[:2], key)
    lock = d + ".lock"
    os.makedirs(os.path.dirname(d), exist_ok=True)
    while True:
        if os.path.exists(os.path.join(d, "meta.json")):
            with open(os.path.join(d, "meta.json")) as f:
                secs = json.load(f)["seconds"]
            print(f"fmcache: hit {name} {key} (computed once in {secs:.0f} s)", file=sys.stderr)
            return d, True, secs
        try:
            os.mkdir(lock)
            break
        except FileExistsError:
            try:
                pid = int(open(os.path.join(lock, "pid")).read())
            except (OSError, ValueError):
                pid = None
            if pid is not None and not alive(pid):
                shutil.rmtree(lock, ignore_errors=True)
            else:
                time.sleep(2)
    try:
        with open(os.path.join(lock, "pid"), "w") as f:
            f.write(str(os.getpid()))
        work = tempfile.mkdtemp(dir=os.path.dirname(d), prefix=".work-")
        t0 = time.time()
        compute(work)
        secs = time.time() - t0
        with open(os.path.join(work, "meta.json"), "w") as f:
            json.dump({"key": key, "spec": key_spec, "seconds": secs, "created": time.strftime("%Y-%m-%dT%H:%M:%S")}, f)
        os.replace(work, d)
    finally:
        shutil.rmtree(lock, ignore_errors=True)
    print(f"fmcache: computed {name} {key} in {secs:.0f} s", file=sys.stderr)
    return d, False, secs


def cached_call(name, fn, args, version="", replicate=1):
    key_spec = {"fun": name, "version": version, "args": spec(args), "replicate": int(replicate)}

    def compute(work):
        with open(os.path.join(work, "result.pkl"), "wb") as f:
            pickle.dump(fn(**args), f)

    d, hit, secs = run_cached(name, key_spec, compute)
    with open(os.path.join(d, "result.pkl"), "rb") as f:
        return pickle.load(f)


def finemap_version():
    out = subprocess.run([FINEMAP, "--help"], capture_output=True, text=True).stdout
    line = next((l for l in out.splitlines() if "FINEMAP v" in l), "")
    return line.strip()


def cmd_finemap(a, extra):
    key_spec = {"fun": "finemap --sss", "version": finemap_version(), "z_md5": md5_file(a.z), "ld_md5": md5_file(a.ld),
                "n_samples": a.n_samples, "args": extra, "replicate": a.replicate}

    def compute(work):
        os.symlink(os.path.abspath(a.z), os.path.join(work, "data.z"))
        os.symlink(os.path.abspath(a.ld), os.path.join(work, "data.ld"))
        with open(os.path.join(work, "master"), "w") as f:
            f.write("z;ld;snp;config;cred;log;n_samples\n")
            f.write(f"data.z;data.ld;data.snp;data.config;data.cred;data.log;{a.n_samples}\n")
        subprocess.run([FINEMAP, "--sss", "--in-files", "master", *extra], cwd=work, check=True,
                       stdout=open(os.path.join(work, "stdout.txt"), "w"), stderr=subprocess.STDOUT)
        for f in ("data.z", "data.ld"):
            os.remove(os.path.join(work, f))

    d, hit, secs = run_cached("finemap", key_spec, compute)
    for f in glob.glob(os.path.join(d, "data.*")):
        shutil.copyfile(f, a.out + f[len(os.path.join(d, "data")):])
    shutil.copyfile(os.path.join(d, "stdout.txt"), a.out + ".stdout")
    print(json.dumps({"hit": hit, "key": os.path.basename(d), "seconds": round(secs, 1)}))


def main():
    argv = sys.argv[1:]
    extra = argv[argv.index("--") + 1:] if "--" in argv else []
    argv = argv[:argv.index("--")] if "--" in argv else argv
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    fm = sub.add_parser("finemap")
    fm.add_argument("--z", required=True), fm.add_argument("--ld", required=True)
    fm.add_argument("--n-samples", type=int, required=True), fm.add_argument("--out", required=True)
    fm.add_argument("--replicate", type=int, default=1)
    a = ap.parse_args(argv)
    cmd_finemap(a, extra)


if __name__ == "__main__":
    main()
