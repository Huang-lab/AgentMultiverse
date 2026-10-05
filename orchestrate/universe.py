#!/usr/bin/env python3
"""Expand a cohort spec into universe folders, each with a frozen guide, a prompt and a spec.

A cohort is a set of blinded agent runs ("universes") on one task.
  design "free"       every universe gets the same prompt and chooses everything itself
  design "factorial"  every combination of the pinned choice axes is a universe; all other choices stay free
Each universe is run by one agent started from the printed prompt; this script launches nothing.

Usage: universe.py SPEC.json [--root results] [--guide CLAUDE.md]

SPEC.json
  cohort      folder name under ROOT
  task        what the agent must do (domain text; never earlier findings)
  design      "free" | "factorial"
  replicates  independent universes per combination (default 1)
  scope       list of items every universe must cover exactly (regions, files, samples...); keeps universes comparable
  pins        {axis: {level: "instruction given to the agent", ...}, ...}   (factorial only)
  note        free text kept in cohort.json (not shown to agents)

Writes ROOT/COHORT/{cohort.json,cohort.tsv} and, per universe, ROOT/COHORT/UNIVERSE/{GUIDE.md,PROMPT.txt,spec.json}.
"""
import argparse
import hashlib
import itertools
import json
import os
import re
import shutil
import sys


def md5(path):
    return hashlib.md5(open(path, "rb").read()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("spec")
    p.add_argument("--root", default="results")
    p.add_argument("--guide", default="CLAUDE.md")
    a = p.parse_args()

    spec = json.load(open(a.spec))
    cohort, task = spec["cohort"], spec["task"]
    design = spec.get("design", "free")
    reps = int(spec.get("replicates", 1))
    pins = spec.get("pins") or {}
    if not re.fullmatch(r"[A-Za-z0-9._-]+", cohort):
        sys.exit("cohort must be a plain folder name")
    if design == "free" and pins:
        sys.exit("a free design has no pins; use design factorial")
    if design == "factorial" and not pins:
        sys.exit("a factorial design needs pins")
    root = os.path.join(a.root, cohort)
    if os.path.exists(root):
        sys.exit(f"{root} exists; choose a new cohort name")

    guide_text = open(a.guide, encoding="utf-8").read()
    m = re.search(r"[Gg]uide version (\d+)", guide_text)
    guide_version = m.group(1) if m else ""

    axes = list(pins)
    combos = list(itertools.product(*[list(pins[x]) for x in axes])) if axes else [()]
    rows, prompts = [], []
    n = 0
    for combo in combos:
        for rep in range(1, reps + 1):
            n += 1
            uid = f"u{n:02d}"
            run_dir = os.path.join(root, uid)
            os.makedirs(run_dir)
            shutil.copyfile(a.guide, os.path.join(run_dir, "GUIDE.md"))
            choices = dict(zip(axes, combo))
            lines = [task.strip(), ""]
            if spec.get("scope"):
                lines += ["Scope (cover exactly this, no more and no less):"] + [f"- {x}" for x in spec["scope"]] + [""]
            lines += [
                     f"Read {os.path.join(run_dir, 'GUIDE.md')} for the resources available and work only from it and the resources it lists.",
                     f"Keep all your outputs in {run_dir}/ and record your work with tools/runlog.py {run_dir} (actions, decisions with reasons, outputs, results).",
                      "Log every choice that could change the results as its own decision entry: runlog.py <run_dir> decision \"what you chose\" --why \"alternatives and reason\" --axis \"the dimension decided\" --choice \"the value\"."]
            if choices:
                lines += ["", "Requirements for this analysis:"] + [f"- {pins[x][v]}" for x, v in choices.items()]
            prompt = "\n".join(lines) + "\n"
            open(os.path.join(run_dir, "PROMPT.txt"), "w").write(prompt)
            json.dump({"cohort": cohort, "universe": uid, "replicate": rep, "choices": choices,
                       "guide_version": guide_version, "guide_md5": md5(a.guide)},
                      open(os.path.join(run_dir, "spec.json"), "w"), indent=2)
            rows.append([uid, str(rep)] + [choices.get(x, "") for x in axes])
            prompts.append((uid, run_dir, prompt))

    json.dump({**spec, "guide_version": guide_version, "guide_md5": md5(a.guide), "n_universes": n},
              open(os.path.join(root, "cohort.json"), "w"), indent=2)
    with open(os.path.join(root, "cohort.tsv"), "w") as f:
        f.write("\t".join(["universe", "replicate"] + axes) + "\n")
        f.writelines("\t".join(r) + "\n" for r in rows)
    print(f"{cohort}: {n} universes, design {design}, guide version {guide_version or '?'} -> {root}/")
    for uid, run_dir, _ in prompts:
        print(f"  {uid}  prompt: {run_dir}/PROMPT.txt")


if __name__ == "__main__":
    main()
