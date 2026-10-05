#!/usr/bin/env python3
"""Append-only record of what an analysis run did, decided and produced.

Each run folder gets three files, all append-only:
  log.jsonl     one JSON object per entry (machine-readable)
  LOG.md        the same entries, human-readable
  manifest.tsv  every file registered with `output`, with size and md5

Usage:
  runlog.py RUN_DIR init "title" [--why TEXT] [--data JSON]
  runlog.py RUN_DIR action   "what was done"      [--cmd TEXT] [--file PATH ...]
  runlog.py RUN_DIR decision "what was chosen"    --why "alternatives and reason" [--axis NAME --choice VALUE]
  runlog.py RUN_DIR output   "what the file holds" --file PATH [--file PATH ...]
  runlog.py RUN_DIR issue    "problem, deviation or failure, and how it was handled"
  runlog.py RUN_DIR result   "key finding with numbers" [--data JSON]
  runlog.py RUN_DIR show
  runlog.py ROOT collect     merge every log.jsonl under ROOT into ROOT/all_logs.tsv
"""
import argparse
import datetime
import hashlib
import json
import os
import sys

KINDS = ("action", "decision", "output", "issue", "result")
MD5_MAX_BYTES = 2 * 1024**3


def now():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def human(n):
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024


def append(path, text):
    with open(path, "a", encoding="utf-8") as f:
        f.write(text)


def describe_file(run_dir, path):
    if not os.path.exists(path):
        sys.exit(f"runlog: file not found: {path}")
    size = os.path.getsize(path)
    rel = os.path.relpath(path, run_dir)
    if rel.startswith(".."):
        rel = os.path.relpath(path)
    return {"path": rel, "bytes": size, "md5": md5sum(path) if size <= MD5_MAX_BYTES else ""}


def log_entry(args):
    run_dir = args.run_dir
    if args.kind == "init":
        os.makedirs(run_dir, exist_ok=True)
    elif not os.path.isdir(run_dir):
        sys.exit(f"runlog: {run_dir} does not exist; start with `runlog.py {run_dir} init \"title\"`")
    if args.kind == "decision" and not args.why:
        sys.exit("runlog: a decision needs --why (the alternatives considered and the reason)")
    if bool(args.axis) != bool(args.choice):
        sys.exit("runlog: --axis and --choice go together (the dimension decided, and the value chosen)")
    if args.axis and args.kind != "decision":
        sys.exit("runlog: --axis/--choice only apply to a decision")
    if args.kind == "output" and not args.file:
        sys.exit("runlog: an output needs at least one --file")

    files = [describe_file(run_dir, p) for p in (args.file or [])]
    data = json.loads(args.data) if args.data else None
    entry = {"ts": now(), "run": os.path.basename(os.path.abspath(run_dir)), "kind": args.kind,
             "msg": args.message, "why": args.why, "cmd": args.cmd, "files": files, "data": data,
             "axis": args.axis, "choice": args.choice}
    append(os.path.join(run_dir, "log.jsonl"), json.dumps(entry) + "\n")

    stamp = entry["ts"][:19].replace("T", " ")
    if args.kind == "init":
        md = f"# {args.message}\n\nRun folder: `{os.path.relpath(run_dir)}`, started {stamp}.\n"
        if args.why:
            md += f"\n{args.why}\n"
        if data:
            md += "\n```json\n" + json.dumps(data, indent=2) + "\n```\n"
        md += "\n## Log\n\n"
    else:
        md = f"- {stamp} **{args.kind.upper()}** {args.message}"
        if args.why:
            md += f" _Why:_ {args.why}"
        if args.axis:
            md += f" [{args.axis} = {args.choice}]"
        if args.cmd:
            md += f"\n  `{args.cmd}`"
        for f in files:
            md += f"\n  `{f['path']}` ({human(f['bytes'])}{', md5 ' + f['md5'][:10] if f['md5'] else ''})"
        if data:
            md += f"\n  `{json.dumps(data)}`"
        md += "\n"
    append(os.path.join(run_dir, "LOG.md"), md)

    if args.kind == "output":
        manifest = os.path.join(run_dir, "manifest.tsv")
        if not os.path.exists(manifest):
            append(manifest, "path\tbytes\tmd5\tlogged_at\tdescription\n")
        for f in files:
            append(manifest, f"{f['path']}\t{f['bytes']}\t{f['md5']}\t{entry['ts']}\t{args.message}\n")


def collect(root):
    rows = []
    for dirpath, _, filenames in os.walk(root):
        if "log.jsonl" in filenames:
            with open(os.path.join(dirpath, "log.jsonl"), encoding="utf-8") as f:
                rows += [json.loads(line) for line in f if line.strip()]
    rows.sort(key=lambda e: (e["run"], e["ts"]))
    clean = lambda s: (s or "").replace("\t", " ").replace("\n", " ")
    out = os.path.join(root, "all_logs.tsv")
    with open(out, "w", encoding="utf-8") as f:
        f.write("run\tts\tkind\taxis\tchoice\tmsg\twhy\tcmd\tfiles\n")
        for e in rows:
            files = ",".join(x["path"] for x in e.get("files") or [])
            f.write("\t".join([e["run"], e["ts"], e["kind"], clean(e.get("axis")), clean(e.get("choice")), clean(e["msg"]), clean(e.get("why")),
                               clean(e.get("cmd")), files]) + "\n")
    counts = {}
    for e in rows:
        counts.setdefault(e["run"], {}).setdefault(e["kind"], 0)
        counts[e["run"]][e["kind"]] += 1
    for run, c in sorted(counts.items()):
        print(run, " ".join(f"{k}={v}" for k, v in sorted(c.items())))
    print(f"wrote {out} ({len(rows)} entries)")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("run_dir")
    p.add_argument("kind", choices=("init",) + KINDS + ("show", "collect"))
    p.add_argument("message", nargs="?")
    p.add_argument("--why")
    p.add_argument("--cmd")
    p.add_argument("--axis", help="decision only: the dimension decided, e.g. the kind of choice, in your own words")
    p.add_argument("--choice", help="decision only: the value chosen on that dimension")
    p.add_argument("--file", action="append")
    p.add_argument("--data", help="JSON object with structured details (numbers, parameters)")
    args = p.parse_args()

    if args.kind == "show":
        with open(os.path.join(args.run_dir, "LOG.md"), encoding="utf-8") as f:
            print(f.read(), end="")
    elif args.kind == "collect":
        collect(args.run_dir)
    else:
        if not args.message:
            p.error(f"{args.kind} needs a message")
        log_entry(args)


if __name__ == "__main__":
    main()
