#!/usr/bin/env python3
"""Download a file once, safely, when several analyses may ask for it at the same time.

Usage:
  tools/fetch.py URL DEST

If DEST exists it is left alone. Otherwise the download goes to a partial file in a hidden
.fetch/ folder next to DEST (resumed if an earlier attempt was interrupted) and is moved to DEST
only when its size matches the server's, so DEST is never half-written. A second request for the
same DEST waits for the first one and then reuses its file. Progress is printed every minute.
"""
import fcntl
import os
import re
import subprocess
import sys
import threading
import time

REPORT_EVERY_SECONDS = 60.0


def remote_size(url):
    """Content-Length of the final response after redirects, or None if the server does not say."""
    head = subprocess.run(["curl", "-fsSIL", "--retry", "5", url], capture_output=True, text=True)
    sizes = re.findall(r"(?im)^content-length:\s*(\d+)", head.stdout)
    return int(sizes[-1]) if head.returncode == 0 and sizes else None


def report_progress(part, total, t0, done):
    while not done.wait(REPORT_EVERY_SECONDS):
        size = os.path.getsize(part) if os.path.exists(part) else 0
        of = f" of {total / 1e6:,.0f}" if total else ""
        print(f"fetch: {os.path.basename(part)[:-5]} {size / 1e6:,.0f}{of} MB after {time.time() - t0:.0f} s",
              file=sys.stderr, flush=True)


def main():
    if len(sys.argv) == 2 and sys.argv[1] in ("-h", "--help"):
        print(__doc__.strip())
        return
    if len(sys.argv) != 3:
        sys.exit("usage: tools/fetch.py URL DEST")
    url, dest = sys.argv[1], os.path.abspath(sys.argv[2])
    if os.path.exists(dest):
        print(f"fetch: {sys.argv[2]} already present")
        return
    work = os.path.join(os.path.dirname(dest), ".fetch")
    os.makedirs(work, exist_ok=True)
    name = os.path.basename(dest)
    lock = os.open(os.path.join(work, name + ".lock"), os.O_RDWR | os.O_CREAT, 0o666)
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        print(f"fetch: another process is downloading {name}; waiting for it", file=sys.stderr, flush=True)
        fcntl.flock(lock, fcntl.LOCK_EX)
    if os.path.exists(dest):
        print(f"fetch: {sys.argv[2]} already present")
        return

    part = os.path.join(work, name + ".part")
    total = remote_size(url)
    if total is not None and os.path.exists(part) and os.path.getsize(part) > total:
        os.remove(part)
    t0 = time.time()
    if total is None or not os.path.exists(part) or os.path.getsize(part) < total:
        done = threading.Event()
        threading.Thread(target=report_progress, args=(part, total, t0, done), daemon=True).start()
        rc = subprocess.run(["curl", "-fsSL", "--retry", "5", "--retry-delay", "10", "-C", "-", "-o", part, url]).returncode
        done.set()
        if rc != 0:
            kept = "; the partial download is kept and the next call resumes it" if os.path.exists(part) else ""
            sys.exit(f"fetch: curl failed with exit code {rc}{kept}")
    size = os.path.getsize(part)
    if total is not None and size != total:
        sys.exit(f"fetch: got {size} bytes but the server reports {total}; the next call resumes the download")
    os.replace(part, dest)
    secs = time.time() - t0
    print(f"fetch: {sys.argv[2]} ({size / 1e6:,.0f} MB in {secs:.0f} s)")


if __name__ == "__main__":
    main()
