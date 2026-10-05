#!/usr/bin/env python3
"""Shared memory gate for analyses that run side by side on this machine.

Usage:
  tools/memgate.py GB -- COMMAND [ARGS ...]   wait until GB of the shared budget is free, run COMMAND, release it
  tools/memgate.py status                      print the budget and what is held or queued

Cooperative: each analysis routes its memory-heavy steps through the gate, so their combined
reservations stay within the budget. GB is your estimate of the step's peak RAM.
Requests are served first come, first served; a request larger than the whole budget waits
until nothing else holds memory and then runs alone. Reservations of processes that died are
dropped automatically. The exit code is COMMAND's.
Gate single heavy steps rather than whole pipelines, so other analyses are not kept waiting.
For a shell pipeline: tools/memgate.py 4 -- bash -c '...'

The budget (GB) is read from .memgate/budget in the project root at every check (default 8),
so it can be changed while jobs are queued. Each finished job is appended to .memgate/history.tsv.
"""
import fcntl
import json
import os
import signal
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STATE = os.path.join(ROOT, ".memgate")
DEFAULT_BUDGET_GB = 8.0
POLL_SECONDS = 2.0
REPORT_EVERY_SECONDS = 60.0
USAGE = "usage: tools/memgate.py GB -- COMMAND [ARGS ...] | tools/memgate.py status"


def budget():
    try:
        with open(os.path.join(STATE, "budget")) as f:
            return float(f.read().strip())
    except (OSError, ValueError):
        return DEFAULT_BUDGET_GB


def alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


class Ledger:
    """Exclusive, locked view of the reservations; dead holders are dropped on entry, changes saved on exit."""

    def __enter__(self):
        os.makedirs(STATE, exist_ok=True)
        self.fd = os.open(os.path.join(STATE, "lock"), os.O_RDWR | os.O_CREAT, 0o666)
        fcntl.flock(self.fd, fcntl.LOCK_EX)
        self.path = os.path.join(STATE, "ledger.json")
        try:
            with open(self.path) as f:
                entries = json.load(f)
        except (OSError, ValueError):
            entries = []
        self.entries = [e for e in entries if alive(e["pid"])]
        return self

    def __exit__(self, *exc):
        tmp = self.path + ".tmp"
        with open(tmp, "w") as f:
            json.dump(self.entries, f)
        os.replace(tmp, self.path)
        fcntl.flock(self.fd, fcntl.LOCK_UN)
        os.close(self.fd)

    def running(self):
        return [e for e in self.entries if e["started"]]

    def waiting(self):
        return sorted((e for e in self.entries if not e["started"]), key=lambda e: (e["requested"], e["pid"]))


def acquire(gb):
    pid = os.getpid()
    requested = time.time()
    last_report = requested
    while True:
        with Ledger() as ledger:
            if not any(e["pid"] == pid for e in ledger.entries):
                ledger.entries.append({"pid": pid, "gb": gb, "requested": requested, "started": None})
            running, waiting = ledger.running(), ledger.waiting()
            held, cap = sum(e["gb"] for e in running), budget()
            ahead = [e["pid"] for e in waiting].index(pid)
            if ahead == 0 and (held + gb <= cap or not running):
                for e in ledger.entries:
                    if e["pid"] == pid:
                        e["started"] = time.time()
                break
        now = time.time()
        if now - last_report >= REPORT_EVERY_SECONDS:
            print(f"memgate: waiting {now - requested:.0f} s for {gb:g} GB; {held:g} of {cap:g} GB held by "
                  f"{len(running)} job(s); {ahead} request(s) ahead", file=sys.stderr, flush=True)
            last_report = now
        time.sleep(POLL_SECONDS)
    return requested, time.time()


def release(gb, requested, started, rc):
    pid = os.getpid()
    with Ledger() as ledger:
        ledger.entries = [e for e in ledger.entries if e["pid"] != pid]
    history = os.path.join(STATE, "history.tsv")
    new = not os.path.exists(history)
    ended = time.time()
    stamp = lambda t: time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(t))
    with open(history, "a") as f:
        if new:
            f.write("requested\tstarted\tended\twait_s\trun_s\tgb\texit\tcwd\n")
        f.write(f"{stamp(requested)}\t{stamp(started)}\t{stamp(ended)}\t{started - requested:.0f}\t"
                f"{ended - started:.0f}\t{gb:g}\t{rc}\t{os.getcwd()}\n")


def run(gb, cmd):
    requested, started = acquire(gb)
    waited = started - requested
    if waited >= REPORT_EVERY_SECONDS:
        print(f"memgate: got {gb:g} GB after {waited:.0f} s", file=sys.stderr, flush=True)
    rc = 127
    try:
        child = subprocess.Popen(cmd)
        forward = lambda signum, frame: child.send_signal(signum)
        for s in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
            signal.signal(s, forward)
        rc = child.wait()
    except OSError as e:
        print(f"memgate: cannot run {cmd[0]}: {e}", file=sys.stderr)
    finally:
        release(gb, requested, started, rc)
    sys.exit(rc if rc >= 0 else 128 - rc)


def status():
    with Ledger() as ledger:
        running, waiting, cap = ledger.running(), ledger.waiting(), budget()
    held = sum(e["gb"] for e in running)
    queued = sum(e["gb"] for e in waiting)
    print(f"budget {cap:g} GB; held {held:g} GB by {len(running)} job(s); "
          f"{len(waiting)} waiting for {queued:g} GB")


def main():
    args = sys.argv[1:]
    if not args or args[0] in ("-h", "--help"):
        print(__doc__.strip())
        return
    if args == ["status"]:
        status()
        return
    try:
        gb = float(args[0])
    except ValueError:
        sys.exit(USAGE)
    cmd = args[2:] if len(args) > 1 and args[1] == "--" else args[1:]
    if gb <= 0 or not cmd:
        sys.exit(USAGE)
    run(gb, cmd)


if __name__ == "__main__":
    main()
