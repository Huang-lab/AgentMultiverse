# Agentic Multiverse Analysis (AMA)

## Definition

Agentic Multiverse Analysis runs many independent, blinded AI analysts on the same data with the same resources.
Each run is fully recorded.
The spread of results across runs is the measurement.
That spread is attributed to the analytic choices that produced it.
Given ground truth, it can also be used to find which choices work best.

It extends two existing ideas to AI analysts.
Multiverse analysis (Steegen et al. 2016) enumerates the defensible analytic choices and reports results across all of them.
Many-analysts studies (Silberzahn et al. 2018) give one dataset to independent teams and compare their conclusions.
Both are limited by human effort.
Agents make a cohort of dozens of recorded analysts cheap enough to be a routine experiment.

## Glossary

| Term | Meaning |
|---|---|
| Universe | One recorded agent run: prompt, frozen guide, log, outputs, transcript |
| Cohort | A set of universes started from the same task and the same guide version |
| Guide | The only document the agent reads; it lists resources and does not prescribe methods |
| Workbench | The software, reference data and shared prep that every universe can use |
| Domain pack | The part of the workbench and guide specific to one field (today: GWAS fine-mapping) |
| Choice axis | A dimension on which analysts can differ, such as reference panel, variant filter, method, prior |
| Pin | A requirement in the prompt that fixes one axis; unpinned axes stay free |
| Decision ledger | The decisions in `log.jsonl`, each with its reason and optionally `axis` and `choice` |
| Spread | Variation of an output (for example a PIP) across the universes of a cohort |
| Attribution | Share of the spread explained by each pinned axis (eta squared) |
| Benchmark | A task with known truth, used to score universes |

## Modes

1. **Free.** No pins. Blinded agents choose everything. This measures the natural variation between analysts.
2. **Directed (factorial).** Some axes are pinned in a full factorial design, replicated. The rest stay free. This attributes spread to the pinned axes.
3. **Searched.** Run directed cohorts on a benchmark, keep the winning levels, pin them, and open new axes. This finds the best way to do the analysis.

Only mode 3 needs ground truth.
Agreement between universes is not accuracy: every universe can share the same bias.

## Principles

1. **Blind.** An agent sees only the guide and the resources. Earlier results, prompts, logs and the experimental design stay hidden.
2. **Record.** Every action, decision (with reason) and output is logged in a fixed format, so universes compare mechanically.
3. **Share the mechanical work.** Precompute everything decision-free once (`prep/`), keep all of it, and flag instead of filtering. Parallelize only the decisions.
4. **Choices are arguments.** Shared tools take every analytic choice as a required argument with no default.
5. **Freeze the guide.** Each universe keeps a byte-identical copy of the guide version it used.
6. **Same output tables.** Universes write the same result tables, so they can be joined.

## Lifecycle of a cohort

1. Write a cohort spec (`orchestrate/universe.py` documents the fields).
2. `orchestrate/universe.py SPEC.json` creates one folder per universe with `GUIDE.md`, `PROMPT.txt` and `spec.json`.
3. Start one blinded agent per universe from its `PROMPT.txt`, respecting the machine limits (RAM, network, concurrency).
4. Each agent logs with `tools/runlog.py` and writes the standard tables.
5. Save each agent's transcript into its universe folder.
6. `orchestrate/compare.py COHORT_DIR` writes `spread.tsv`, `agreement.tsv`, `attribution.tsv` and, with `--truth`, `score.tsv`.
7. `tools/archive_study.sh COHORT` zips the cohort into the blinded archive.

## Status

Built: blinding, run record, shared prep, result cache, resource gate, cohort launcher, comparator, archive.
Not yet built:
- A scored benchmark (known causal variants or simulated loci). Without it, mode 3 is not possible.
- A controlled vocabulary of choice axes. Today axes are optional free text in the decision ledger.
- Automatic agent launching. The orchestrator starts agents from the prompts by hand.
- Tests for how much the pinned wording of an axis biases the agent beyond the pin itself.

## Adding a domain

1. Add the domain's software and reference data under the workbench layout and list them in `CLAUDE.md`.
2. Put decision-free precomputation in `prep/<study>/` and make shared loaders take choices as required arguments.
3. Define the standard result tables and a checker, as `tools/check_outputs.py` does for fine-mapping.
4. Provide a benchmark with truth if you want searched mode.
5. Reuse `tools/runlog.py`, `tools/memgate.py`, `tools/fetch.py`, `tools/fmcache.*`, `orchestrate/` and the blinding hooks unchanged.
