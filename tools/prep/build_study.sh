#!/usr/bin/env bash
# Build the shared, decision-free inputs for one study once, so parallel analyses reuse them.
# Usage (from the project root): tools/prep/build_study.sh STUDY [FLANK_BP]
#   reads data/sumstats/STUDY/*.tsv.gz (GWAS-SSF, GRCh38) and writes prep/STUDY/.
# Steps: harmonize (all variants kept, flagged) -> regions and leads -> PolyFun priors -> LD per region for each
# panel (1000G hg38 EUR, UK Biobank PolyFun, Pan-UKBB EUR; the last needs env_hail/ and downloads for hours).
# Each step skips work that is already done, so an interrupted build can be rerun.
set -euo pipefail
study=${1:?usage: tools/prep/build_study.sh STUDY [FLANK_BP]}
flank=${2:-500000}
export PATH="$PWD/bin:$PWD/env/bin:$PATH"
out=prep/$study
mkdir -p "$out/sumstats" "$out/.work"

for f in data/sumstats/$study/*.tsv.gz; do
  s=$(basename "$f" .tsv.gz)
  [ -e "$out/sumstats/$s.parquet" ] || python tools/prep/harmonize.py "$out/sumstats" "$f"
done
[ -e "$out/regions.tsv" ] || python tools/prep/regions.py "$out" "$flank"
[ -e ref/gencode/genes.v50.tsv.gz ] || python tools/prep/genes.py ref/gencode/genes.v50.tsv.gz
[ -e "$out/priors/polyfun_snpvar_bin.parquet" ] || python tools/prep/priors.py "$out"
python tools/prep/ld.py 1kg "$out" > "$out/.work/ld_1kg.log" 2>&1 &
python tools/prep/ld.py panukb "$out" > "$out/.work/ld_panukb.log" 2>&1 &
python tools/prep/ld.py ukb "$out" > "$out/.work/ld_ukb.log" 2>&1
wait
(cd "$out" && find sumstats ld priors regions.tsv leads.tsv -type f ! -name '.*' ! -name '*.lock' | sort \
  | xargs md5 -r > MD5SUMS)
echo "built $out: $(du -sh "$out" | cut -f1)"
