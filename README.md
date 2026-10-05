# AgentMultiverse

A workbench for Agentic Multiverse Analysis (AMA): many blinded AI analysts run the same analysis their own way, every run is recorded, and the results are compared to see what the analytic choices change and which choices work best.

The framework is domain-neutral.
The first installed domain is GWAS fine-mapping from summary statistics (SuSiE, FINEMAP, CARMA, PolyFun priors, UK Biobank, Pan-UKBB and 1000 Genomes LD).
See `framework/AMA.md` for the definition, glossary, modes and what is still missing.

## Layout

| Path | Role | Domain |
|---|---|---|
| `CLAUDE.md` | The only document agents read: a resource list with the facts needed to use each resource | per domain |
| `FINEMAPPING.md`, `FINEMAPPING.v*.md` | Copy of the live guide and frozen earlier versions | fine-mapping |
| `tools/runlog.py` | Run record: actions, decisions with reasons (optionally `--axis`/`--choice`), outputs with md5 | any |
| `tools/memgate.py`, `tools/fetch.py` | Shared RAM budget; deduplicated resumable downloads | any |
| `tools/fmcache.R`, `tools/fmcache.py` | Content-keyed cache for slow method calls | any |
| `tools/archive_study.sh` | Zip a finished cohort into the blinded archive | any |
| `orchestrate/universe.py` | Expand a cohort spec into universe folders (frozen guide, prompt, spec) | any |
| `orchestrate/compare.py` | Spread, agreement, attribution to pinned axes, score against truth | any (defaults read `pip.tsv`) |
| `tools/prep/`, `prep/<study>/` | Shared decision-free inputs per study | fine-mapping |
| `tools/load_region.*`, `tools/check_outputs.py`, `tools/polyfun/` | Loaders, output checker, PolyFun | fine-mapping |
| `bin/`, `env/`, `env_hail/`, `ref/` | Software and reference data | fine-mapping |
| `results/<cohort>/<universe>/` | Work in progress | any |
| `.claude/` | Blinding rules (deny rules and a Bash hook) and resource caps | any |

## Run a cohort

```bash
export PATH="$PWD/bin:$PWD/env/bin:$PATH"
python orchestrate/universe.py my_cohort.json        # writes results/<cohort>/u01.../{GUIDE.md,PROMPT.txt,spec.json}
# start one blinded agent per universe from its PROMPT.txt, save its transcript in the universe folder
python orchestrate/compare.py results/<cohort> [--truth truth.tsv]
tools/archive_study.sh <cohort>
```

Spec example (factorial over the LD panel, 3 replicates each):

```json
{"cohort": "eadb2026_v4_panel", "design": "factorial", "replicates": 3,
 "task": "Fine-map the genome-wide significant loci of study eadb2026.",
 "pins": {"panel": {"ukb": "Use UK Biobank LD.", "1kg": "Use 1000 Genomes EUR LD."}}}
```

## Blinding

Agents work in this folder, so anything readable here can influence them.
They are meant to see only `CLAUDE.md` and the resources it lists.
Deny rules and a Bash hook in `.claude/` block the archive, session transcripts, the result cache, this README and `framework/`.
The block applies to the orchestrator too, by design.
To edit these documents, lift the rules for the edit or work from a session opened in another folder.
`orchestrate/` is not blocked so the orchestrator can run it, so it should stay generic and hold no findings.

## Machine limits

Built for one laptop: 10 cores, 24 GB RAM, about 3-4 MB/s network.
BLAS and OpenMP are capped at 2 threads per process in `.claude/settings.json`.
Route heavy steps through `tools/memgate.py`.

## Rebuild

The `Rebuild` section of `CLAUDE.md` recreates binaries, environments and reference data (about 20 GB, 20-30 minutes).
