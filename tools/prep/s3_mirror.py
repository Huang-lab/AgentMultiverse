#!/usr/bin/env python3
"""Mirror files from a public S3 bucket over plain HTTPS, keeping the key layout.

Usage: tools/prep/s3_mirror.py BUCKET_URL PREFIX DEST_DIR [--keys FILE] [--jobs N]
  Without --keys, every object under PREFIX is mirrored; with --keys, only the keys listed in FILE.
  Existing files are skipped; downloads go to DEST_DIR/.partial/ and are moved into place when complete,
  so the mirrored tree (which Hail reads directly) never holds half-written or temporary files.
"""
import argparse
import os
import re
import subprocess
import sys
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor

def list_keys(bucket, prefix):
    keys, token = [], None
    while True:
        q = {"list-type": "2", "prefix": prefix, **({"continuation-token": token} if token else {})}
        xml = urllib.request.urlopen(f"{bucket}/?{urllib.parse.urlencode(q)}").read().decode()
        keys += re.findall(r"<Key>([^<]+)</Key>", xml)
        token = re.search(r"<NextContinuationToken>([^<]+)</NextContinuationToken>", xml)
        if not token:
            return keys
        token = token.group(1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bucket"), ap.add_argument("prefix"), ap.add_argument("dest")
    ap.add_argument("--keys"), ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    keys = open(a.keys).read().split() if a.keys else list_keys(a.bucket, a.prefix)

    def get(key):
        dest = os.path.join(a.dest, key)
        if os.path.exists(dest):
            return 0
        part = os.path.join(a.dest, ".partial", key.replace("/", "__"))
        os.makedirs(os.path.dirname(part), exist_ok=True)
        rc = subprocess.run(["curl", "-fsSL", "--retry", "5", "-C", "-", "-o", part, f"{a.bucket}/{key}"]).returncode
        if rc == 0:
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            os.replace(part, dest)
        return rc

    with ThreadPoolExecutor(a.jobs) as ex:
        failed = sum(rc != 0 for rc in ex.map(get, keys))
    print(f"{len(keys) - failed} of {len(keys)} objects present under {a.dest}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
