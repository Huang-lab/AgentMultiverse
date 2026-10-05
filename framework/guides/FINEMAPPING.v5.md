# Analysis workbench: fine-mapping

Local software and reference data for analysing data in this folder; the installed domain is fine-mapping GWAS loci from summary statistics.
The run record, memory gate, download and result-cache tools below do not depend on the domain.
Everything lives in this folder, and all paths below are relative to the project root.
Built and verified on macOS arm64 on 2026-09-26; the Rebuild section recreates it on any machine.
Guide version 5 (2026-10-02); an analysis that relies on this guide keeps a copy of the version it used with its results.

```bash
export PATH="$PWD/bin:$PWD/env/bin:$PATH"   # python, Rscript, plink, plink2, finemap, liftOver, bcftools, tabix
```

## What is here

| Path | Contents |
|---|---|
| `bin/plink2`, `bin/plink` | PLINK 2.0 alpha 7.8 (arm64) and PLINK 1.9 (universal) |
| `bin/finemap` + `bin/lib/` | FINEMAP 1.4.2, macOS x86_64 build run under Rosetta; its GCC, OpenMP and zstd libraries are bundled in `bin/lib`, so keep the two together |
| `bin/liftOver` | UCSC liftOver (arm64) |
| `env/` | conda env: Python 3.12 (numpy, scipy, pandas 2, pyarrow, scikit-learn, pandas-plink, matplotlib), R 4.5 (susieR 0.14.2 with Rfast, CARMA 1.0, data.table with R.utils so `fread` reads `.gz`, dplyr, jsonlite), bcftools and htslib 1.24 |
| `env_hail/` | Separate env for Hail 0.2.139 (Python 3.11, Java 11), only needed to read Pan-UKBB LD: `env_hail/bin/python` |
| `tools/polyfun/` | PolyFun: `snpvar_meta.chr*.parquet` (precomputed functional priors, hg19, ~19M UKB SNPs), `ukb_regions.tsv.gz` (index of UK Biobank LD windows), `finemapper.py` (`load_ld_npz`, SuSiE/FINEMAP wrapper) |
| `tools/runlog.py` | Run record: `runlog.py <run_dir> {init,action,decision,issue,result,output} "..."` appends to `LOG.md`, `log.jsonl` and `manifest.tsv` (outputs with md5); a decision needs `--why`, and may name the dimension decided and the value chosen with `--axis` and `--choice` |
| `tools/memgate.py` | Shared RAM budget for analyses running side by side: `memgate.py GB -- COMMAND` waits until GB is free, then runs COMMAND; `memgate.py status` shows what is held |
| `tools/fetch.py` | `fetch.py URL DEST`: downloads once, resumes, and makes concurrent requests for the same file wait for one download |
| `tools/load_region.py`, `tools/load_region.R` | `load_region(prep, sumstats, region, panel, ambiguous, min_n_frac, min_maf)`: z aligned to the panel, the LD submatrix and a variant table from `prep/`; every choice is a required argument, and `counts` says how many variants each step removed |
| `tools/fmcache.R`, `tools/fmcache.py` | Shared result cache (see Working alongside other analyses) |
| `tools/check_outputs.py` | `check_outputs.py <run_dir> --prep prep/<study>` validates the standard result tables |
| `tools/prep/` | Builds the shared inputs under `prep/` (see Shared prep); `build_study.sh STUDY` runs it all |
| `prep/<study>/` | Shared, precomputed inputs per study: harmonized sumstats, regions, leads, and per-region LD for each panel |
| `ref/1kg/hg19/` | 1000 Genomes phase 3, GRCh37, 2504 samples: `all_phase3.{pgen,pvar.zst,psam}`, IDs are rsIDs |
| `ref/1kg/hg38/` | 1000 Genomes NYGC 30x, GRCh38, 3202 samples including relatives: `all_hg38.{pgen,pvar.zst,psam}`, IDs are dbSNP 156 rsIDs |
| `ref/1kg/*/{EUR,AFR,EAS,SAS,AMR}.keep` | Unrelated samples per superpopulation (hg19: 503/652/504/484/347, hg38: 525/680/506/511/351) |
| `ref/ukbb_ld/` | Cache for UK Biobank LD windows (`<name>.npz` + `<name>.gz`); fetch missing ones with `tools/fetch.py` |
| `ref/panukb_ld/` | Pan-UKBB EUR LD (about 420K European-ancestry UK Biobank samples, GRCh37, released 2021): `UKBB.EUR.variants.parquet` maps all 23,960,350 matrix indices (`idx`) to GRCh37 `chr_hg19`, `bp_hg19`, `ref`, `alt`, `alt_freq`, `rsid`; the block matrix itself stays on S3 (`https://pan-ukb-us-east-1.s3.amazonaws.com/ld_release/UKBB.EUR.ldadj.bm`, 4096 x 4096 float64 blocks, upper triangle, signs relative to `alt`) and `tools/prep/ld.py panukb` downloads only the blocks a region needs |
| `ref/ld_blocks/` | Approximately independent LD blocks: hg19 Berisa-Pickrell (EUR/AFR/ASN), hg38 pyrho (EUR/AFR/EAS/SAS) |
| `ref/liftover/` | UCSC chains `hg19ToHg38`, `hg38ToHg19` |
| `ref/gencode/` | GENCODE v50 basic GTF for GRCh38 and GRCh37 (`lift37`), and `genes.v50.tsv.gz`: one row per gene with GRCh38 and GRCh37 coordinates (1-based, inclusive) |
| `data/sumstats/`, `results/` | Inputs and outputs |

UK Biobank LD (337K British-ancestry samples, GRCh37, 3 Mb windows every 1 Mb) is the best reference for European hg19 summary stats.
A window is typically about 1 GB (up to about 3 GB), so fetch it per locus rather than in bulk.
Pick the row of `tools/polyfun/ukb_regions.tsv.gz` whose window contains the locus near its center, download `<URL_PREFIX>.npz` and `<URL_PREFIX>.gz` into `ref/ukbb_ld/`, and load them with `load_ld_npz("ref/ukbb_ld/<name>")` from `tools/polyfun/finemapper.py`.
That returns the full matrix plus a SNP table, with signs relative to `allele1`.

## Shared prep

For a study with a `prep/<study>/` folder, the expensive and decision-free work is already done once for all analyses: start from these files instead of redoing harmonization, liftOver, LD extraction or downloads.
Nothing in them is filtered beyond what is stated, so every QC choice is still yours.
Treat `prep/` as read-only and write your own outputs to `results/`.

| File | Contents |
|---|---|
| `sumstats/<GCST>.parquet` | One per input file, every row kept, sorted by `chr`, `bp` (GRCh38). Columns: `vkey`, `chr`, `bp`, `ea`, `oa`, `rsid`, `beta`, `se`, `z` (= beta/se), `p` (as published), `log10p` (from z, so it does not underflow), `eaf`, `n` (cases + controls), `n_cases`, `n_controls`, `neff`, `neff_cases`, `neff_controls`, `het_isq`, `het_p`, flags `ambiguous` (A/T or C/G SNP) and `dup_vkey` (vkey repeated in this file), and GRCh37 liftOver `chr_hg19`, `bp_hg19`, `hg19_minus` (mapped to the minus strand; empty when unmapped) |
| `regions.tsv` | Candidate regions: every variant with P < 5e-8 in any of the study's files ±500 kb, overlaps merged, named `chr<c>_<start>_<end>` (GRCh38), with GRCh37 span, `mhc` flag (overlaps chr6:25-34 Mb), strongest `log10p` and where it is |
| `leads.tsv` | One row per file and region: that file's top variant there, `n_sig` (variants with P < 5e-8) and `n_variants` |
| `ld/<region>/<panel>.bin` | Signed r, float32, m x m, row-major, full square (diagonal 1); load as in step 3 below, or `np.memmap(f, np.float32, "r", shape=(m, m))` to read only the rows you need |
| `ld/<region>/<panel>.vars.tsv` | One row per matrix row: `vkey`, `rsid`, `chr`, `bp`, `bp_hg19`, `a1`, `a2`, `a1_freq` (panel frequency of a1; empty for UK Biobank) |
| `priors/polyfun_snpvar_bin.parquet` | `vkey`, `polyfun_snp`, `snpvar_bin`: PolyFun's precomputed prior for every study variant that matches a PolyFun SNP by GRCh37 position and alleles (about 15M), as published, with no flooring or normalization |
| `ld/index.tsv` | Per region and panel: `m`, source window, how many sumstats variants the region has, and a note where the panel does not cover the whole region |
| `MD5SUMS` | Checksums of all of the above |

Panels: `1kg_hg38_EUR` is 1000G NYGC 30x, the 525 unrelated EUR samples, biallelic, MAF ≥ 0.001 in those samples; `ukb_hg19` is the PolyFun UK Biobank window whose span overlaps the region most (all of its variants); `panukb_hg19_EUR` is Pan-UKBB EUR over the region's GRCh37 span (all of its variants, `a1` = GRCh37 `alt`).
The two UK Biobank panels come from largely the same people and agree closely (mean |Δr| about 0.001 on a test region), but differ in samples, variant sets and covariate adjustment.
Each matrix holds only panel variants whose chromosome, position and allele pair match a sumstats variant in any of the study's files, so `m` is at most the region's variant count.
`vkey` is `chr:bp:A:B` in GRCh38 with the two alleles sorted, identical in the sumstats and in every `vars.tsv`, so join on it.
Signs are relative to `a1`, and `a1`/`a2` are given on the GRCh38 forward strand (UK Biobank alleles at minus-strand liftOver sites are complemented back; minus-strand indels are not matched), so a file's z is aligned by multiplying by +1 where `ea == a1` and by -1 where `ea == a2`.
Strand-ambiguous and duplicated variants are in the matrices; drop or keep them as you decide.
PolyFun UK Biobank windows are 3 Mb, so a region wider than about 2 Mb may be only partly covered (see `note`), and there are no windows over most of the MHC; Pan-UKBB covers every region in full.
To use a different region, flank, panel or MAF threshold, subset these matrices or build your own from `ref/` as described below.

## Working alongside other analyses

Several analyses usually run on this machine at once (10 cores, 24 GB RAM, about 3-4 MB/s of network), so share it:

- BLAS and OpenMP are capped at 2 threads per process through the project settings; give plink2 `--threads 2` as well.
  Route steps that need more than about 2 GB of RAM through `tools/memgate.py`.
- Read LD from `prep/` in place (memory-mapped or via `load_region`); do not copy matrices into `results/`, and record the path and md5 of what you used instead.
- Slow method calls go through the shared result cache, so an identical configuration is computed once across all analyses:
  in R, `source("tools/fmcache.R")` and call `cached_call("susieR::susie_rss", list(z = z, R = R, n = n, L = 10))` or `cached_call("CARMA::CARMA", list(...))` exactly as you would `do.call`;
  for FINEMAP, `tools/fmcache.py finemap --z X.z --ld X.ld --n-samples N --out PREFIX -- <finemap args>`.
  The key covers the function, its version and every argument, so a hit is always the result for exactly your inputs; pass `replicate = 2, 3, ...` (`--replicate`) for independent runs of a stochastic method such as FINEMAP.
  Hashing costs about 1 s per 150 MB of matrix, so small SuSiE fits gain little. The cache is reached only through these wrappers.
- Write results in the standard tables so runs can be compared by joining them: `configs.tsv` (`config`, `description`), `pip.tsv` (`config`, `study`, `region`, `method`, `vkey`, `rsid`, `pip`, `cs`) and `cs.tsv` (`config`, `study`, `region`, `method`, `cs`, `size`, `coverage`, `min_abs_r`, `top_vkey`, `top_pip`, `vkeys`), where `config` labels one combination of your choices; `tools/check_outputs.py` documents and checks them.

## Choices are yours

Nothing here prescribes how to analyse the data, and no setting below is recommended.
Each of these is your judgment, so use whatever you think is best for the data in front of you, and record each one with `tools/runlog.py ... decision` (`--axis`, `--choice`):

- Which files and loci to analyse, how a locus is bounded (the prep region, a narrower or wider flank, LD blocks), and how the MHC, inversions and very strong signals are handled.
- Genome build, and which variants enter: alleles, indels, strand-ambiguous SNPs, duplicates, MAF, per-SNP N, imputation quality, and which lead variants QC may remove.
- Which z and which N: beta/se or p, N or effective N, one file at a time or combined.
- LD: the panel (UK Biobank, Pan-UKBB, 1000 Genomes, other ancestry), the window, regularization or shrinkage of R, and in-sample LD where available.
- Method and settings: SuSiE-RSS, FINEMAP, CARMA, or several; number of effects or causal variants, priors (flat or PolyFun), residual variance, credible-set coverage and purity, estimating or fixing parameters, and how many runs to make.
- Diagnostics and what to do with them: `estimate_s_rss`, `kriging_rss`, convergence, spurious sets next to strong signals, and loci you decide not to trust.
- How results are reported: which method or panel is primary, how disagreements between runs are handled, and what is written to the standard tables.

## Approach

These are defaults and checks, not a fixed pipeline: adapt each step to the data, and record every choice that could change the results with `tools/runlog.py`.
For a study under `prep/`, the mechanical parts of steps 1-3 are already done (see Shared prep); the choices in them are not.

1. Harmonize the sumstats to CHR, BP, effect and other allele, BETA and SE (or Z), P, per-SNP N, and EAF.
   Take z from BETA/SE, because published P values underflow at very strong signals.
   Detect the build by looking up a few rsIDs in the `.pvar` files (rs12740374 is chr1:109817590 in hg19 and chr1:109274968 in hg38), and liftOver if needed.
2. Define loci, for example lead variants with P < 5e-8 plus ±0.5-1 Mb (merging overlaps), or the LD blocks.
   Decide early how to treat regions where fine-mapping with reference LD breaks down: the MHC (chr6:25-34 Mb), large inversions, and signals too strong for any reference panel to match.
3. Choose LD: in-sample if available, otherwise the largest ancestry-matched panel (UK Biobank for European hg19 data), otherwise the matching 1000G superpopulation.
   Panel size limits how well strong signals fit, and the panel's MAF threshold decides which rare variants can be fine-mapped at all.
   ```bash
   plink2 --pfile ref/1kg/hg19/all_phase3 vzs --keep ref/1kg/hg19/EUR.keep --chr 1 --from-bp 109317590 --to-bp 110317590 \
          --maf 0.01 --max-alleles 2 --rm-dup exclude-all --make-pgen --out results/locus      # ~3 s
   plink2 --pfile results/locus --r-unphased square bin4 ref-based --out results/locus         # float32 m x m, rows = results/locus.pvar
   ```
   Load it with `np.fromfile(f, np.float32).reshape(m, m)` in Python or `matrix(readBin(f, "double", m*m, size = 4), m, m)` in R.
4. Match variants to the LD panel by position and allele pair (not rsID), align every z to the panel's counted allele (1000G `ALT`, UK Biobank `allele1`), and drop allele mismatches and duplicates.
   A/T and C/G SNPs are only ambiguous when their frequency cannot resolve the strand: dropping all of them can remove the causal variant, so prefer keeping those whose frequency matches the panel and is not close to 0.5.
   When per-SNP N varies, restrict to variants with similar N (see Pitfalls), and look at strong variants that fall just below the cut.
   Keep track of which lead variants the QC removed and why.
5. Check the fit before trusting any result: `susieR::estimate_s_rss` (near 0 means z and LD agree), `kriging_rss` outliers, and convergence.
   Comparing `estimate_s_rss` for the aligned z with the same z under random sign flips is a quick test of both the alignment and the LD matrix.
6. Fine-map, for example with `susie_rss(z = z, R = R, n = n, L = 10)` and 95% credible sets from `susie_get_cs(fit, Xcorr = R, min_abs_corr = 0.5)`.
   Cross-check with FINEMAP (`finemap --sss`, `.z` columns `rsid chromosome position allele1 allele2 maf beta se` with beta = z and se = 1) and with CARMA (`CARMA(list(z), list(R), lambda.list = list(1/sqrt(p)), outlier.switch = TRUE)`), which models reference-LD mismatch.
   For functional priors (hg19, EUR), pass PolyFun `snpvar_bin` as `prior_weights`.
7. Interpret credible sets rather than single PIPs: the top variant moves with the method, prior and panel much more than the set does.
   Low-|z| sets next to a very strong signal usually absorb LD mismatch, and a confident set at a locus whose lead failed QC may only tag the removed variant.
   Report per locus: PIPs, credible sets (size, purity, top variant), `estimate_s_rss`, convergence, agreement between methods, removed leads, and nearby genes from GENCODE.

Pitfalls found while testing this setup:

- `plink2 --r-unphased` signs are relative to the major allele unless you pass `ref-based`; without it the signs will not match the z orientation.
- Meta-analyses have variable per-SNP N, and mixing SNPs with very different N breaks z/LD consistency.
  On the SORT1 test, keeping N ≥ 0.9·max(N) moved `estimate_s_rss` from 0.51 to 0.06 and made SuSiE converge.
  Without such a filter, very strong loci split into many spurious single-variant sets, and the lead can lose its PIP to a neighbour.
- 1000G has only 347-680 unrelated samples per superpopulation, so with very strong signals (|z| > ~30) expect spurious extra signals (FINEMAP piling up at the maximum number of causal SNPs, many singleton SuSiE sets).
  On the SORT1 test, 1000G EUR gave 6 SuSiE credible sets while UK Biobank LD gave 1, so prefer UK Biobank LD for EUR, or in-sample LD.
  At 15 strong loci of a meta-analysis of about a million people, UK Biobank LD also lowered `estimate_s_rss` about 20-fold and removed most low-|z| sets.
- UK Biobank `.npz` files store one triangle of the matrix; `load_ld_npz` mirrors it, but mirror it yourself if you slice the file directly.
  There are no UK Biobank windows over the MHC (roughly chr6:28-33 Mb).
- `susie_rss` holds about 10 dense copies of R even with Rfast installed (1.5 GB at 4,300 variants; 2.2 GB without Rfast).
  Read large LD files only for the rows you need (for example with a memory map), and plan for dense regions: 2 Mb of the MHC holds about 34,000 EUR variants with MAF ≥ 0.01 in 1000G hg38.
- plink2 fails with `--memory 4000` on the 75M-variant hg38 file even for one chromosome (9000 works, as does the default of half the RAM); for many parallel jobs, make per-chromosome subsets first.
  Its rsIDs are shared by different alleles at the same site, so `--rm-dup exclude-all` drops about 7% of records; when matching by position and alleles, reset the IDs with `--set-all-var-ids '@:#:$r:$a'` first.
- FINEMAP 1.4.2 cannot be seeded, and repeated runs can disagree at hard loci, so run it twice and report the stability.
  CARMA 1.0 defaults took from about 1 to over 16 minutes and up to about 3 GB per locus at 1,000-3,000 variants, so budget the time or restrict the variants.
- PolyFun `snpvar_bin` is zero for about half of the variants; floor it at max/100 before normalizing, as PolyFun's `extract_snpvar.py` does, or those variants can never be selected.
- `polyfun/finemapper.py --method susie` needs susieR 0.11.92 because the API changed; with this env use `--method finemap --finemap-exe bin/finemap`, or call susieR directly.
- SuSiE-inf and FINEMAP-inf (github.com/FinucaneLab/fine-mapping-inf) assume in-sample LD and are not installed.

Sanity check: GLGC 2013 LDL (hg19, `https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST002001-GCST003000/GCST002222/jointGwasMc_LDL%20(1).txt.gz`), SORT1 ±500 kb, 1000G EUR LD, N ≥ 0.9·max(N).
Expected: SuSiE puts rs12740374, the experimentally validated causal SNP, alone in CS1 with PIP 0.98, and FINEMAP and CARMA give 0.99 and 1.00.
With UK Biobank LD (`chr1_108000001_111000001`, z aligned to `allele1`), `estimate_s_rss` is 0.04 and SuSiE returns a single credible set, rs12740374 (0.93) plus rs660240.

## Rebuild

Run from the project root; it needs about 20 GB of disk and 20-30 minutes on a fast connection.
On Linux, take the Linux builds from the same pages, and FINEMAP needs no patching.

```bash
mkdir -p bin tools ref/{1kg/hg19,1kg/hg38,ukbb_ld,ld_blocks,liftover,gencode} data/sumstats results
dl() { curl -fL --retry 5 -C - -o "$1" "$2"; }

# Binaries: newest builds from https://www.cog-genomics.org/plink/2.0/ and https://www.cog-genomics.org/plink/,
# FINEMAP from http://www.christianbenner.com/finemap_v1.4.2_{MacOSX,x86_64}.tgz,
# liftOver from https://hgdownload.soe.ucsc.edu/admin/exe/{macOSX.arm64,linux.x86_64}/liftOver

# Python/R env (compilers are only needed to build CARMA)
conda create -y -p env --override-channels -c conda-forge -c bioconda python=3.12 numpy scipy "pandas<3" pyarrow \
  scikit-learn tqdm bitarray networkx pandas-plink matplotlib r-base=4.5 r-susier r-data.table r-r.utils r-matrix r-remotes \
  r-rcpp r-rcpparmadillo r-rcppgsl r-dplyr r-glmnet r-mass r-rfast r-jsonlite gsl bcftools htslib compilers
conda create -y -p env_hail --override-channels -c conda-forge python=3.11 openjdk=11 pip && env_hail/bin/pip install hail pyarrow
PATH="$PWD/env/bin:$PATH" env/bin/Rscript -e 'remotes::install_github("ZikunY/CARMA", upgrade = "never", dependencies = FALSE)'
git clone --depth 1 https://github.com/omerwe/polyfun tools/polyfun

# 1000 Genomes in PLINK 2 format (https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg)
dl ref/1kg/hg19/all_phase3.pgen.zst 'https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst?dl=1'
dl ref/1kg/hg19/all_phase3.pvar.zst 'https://www.dropbox.com/s/c95n8quqwqww4s0/all_phase3_noannot.pvar.zst?dl=1'
dl ref/1kg/hg19/all_phase3.psam 'https://www.dropbox.com/scl/fi/haqvrumpuzfutklstazwk/phase3_corrected.psam?rlkey=0yyifzj2fb863ddbmsv4jkeq6&dl=1'
dl ref/1kg/hg19/deg2_phase3.king.cutoff.out.id 'https://www.dropbox.com/s/zj8d14vv9mp6x3c/deg2_phase3.king.cutoff.out.id?dl=1'
dl ref/1kg/hg38/all_hg38.pgen.zst 'https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1'
dl ref/1kg/hg38/all_hg38.pvar.zst 'https://www.dropbox.com/scl/fi/id642dpdd858uy41og8qi/all_hg38_rs_noannot.pvar.zst?rlkey=sskyiyam1bsqweujjmxqv1h55&dl=1'
dl ref/1kg/hg38/all_hg38.psam 'https://www.dropbox.com/scl/fi/u5udzzaibgyvxzfnjcvjc/hg38_corrected.psam?rlkey=oecjnk4vmbhc8b1p202l0ih4x&dl=1'
dl ref/1kg/hg38/deg2_hg38.king.cutoff.out.id 'https://www.dropbox.com/s/4zhmxpk5oclfplp/deg2_hg38.king.cutoff.out.id?dl=1'
for f in ref/1kg/*/*.pgen.zst; do env/bin/zstd -d --rm "$f"; done
for b in hg19 hg38; do for pop in EUR AFR EAS SAS AMR; do
  awk -v p=$pop 'NR==FNR{r[$1]; next} FNR>1 && $5==p && !($1 in r){print $1}' \
    ref/1kg/$b/deg2_*.king.cutoff.out.id ref/1kg/$b/*.psam > ref/1kg/$b/$pop.keep
done; done

# LD blocks, liftOver chains, genes
for p in EUR AFR ASN; do dl ref/ld_blocks/hg19_${p}_berisa_pickrell.bed https://bitbucket.org/nygcresearch/ldetect-data/raw/master/$p/fourier_ls-all.bed; done
for p in EUR AFR EAS SAS; do dl ref/ld_blocks/hg38_${p}_pyrho.bed https://raw.githubusercontent.com/jmacdon/LDblocks_GRCh38/master/data/pyrho_${p}_LD_blocks.bed; done
dl ref/liftover/hg19ToHg38.over.chain.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
dl ref/liftover/hg38ToHg19.over.chain.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
G=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50
dl ref/gencode/gencode.v50.basic.annotation.gtf.gz $G/gencode.v50.basic.annotation.gtf.gz
dl ref/gencode/gencode.v50lift37.basic.annotation.gtf.gz $G/GRCh37_mapping/gencode.v50lift37.basic.annotation.gtf.gz

# Pan-UKBB EUR variant index (about 1.5 GB of Hail tables), then a parquet of it
env/bin/python tools/prep/s3_mirror.py https://pan-ukb-us-east-1.s3.amazonaws.com ld_release/UKBB.EUR.ldadj.variant.ht/ ref/panukb_ld --jobs 8
printf 'ld_release/UKBB.EUR.ldadj.bm/metadata.json\n' > /tmp/k && env/bin/python tools/prep/s3_mirror.py https://pan-ukb-us-east-1.s3.amazonaws.com x ref/panukb_ld --keys /tmp/k
env_hail/bin/python tools/prep/panukb_index.py

# Shared prep per study (GWAS-SSF files in data/sumstats/<study>/): for eadb2026 about 1 hour and 40 GB,
# plus about 10 hours of Pan-UKBB block downloads (roughly 1 GB per region, deleted after use) and 20 GB
tools/prep/build_study.sh eadb2026
```

FINEMAP on macOS: the official Mac build is x86_64 and links Intel Homebrew's GCC 10 runtime, zstd and libomp from `/usr/local`.
To run it without Intel Homebrew, put x86_64 copies in `bin/lib`: `libstdc++.6`, `libgfortran.5`, `libquadmath.0` and `libgcc_s.1.1` from the MacPorts archive `https://packages.macports.org/libgcc15/libgcc15-15.2.0_0+stdlib_flag.darwin_25.x86_64.tbz2`, plus `libomp`, `libzstd.1` and `libiconv.2` from `CONDA_SUBDIR=osx-64 conda create -p /tmp/x86 -c conda-forge llvm-openmp zstd libiconv`.
Then repoint every `/usr/local/...` dependency listed by `otool -L bin/finemap` to `@executable_path/lib/<name>`, repoint `libstdc++`'s libiconv and `libgfortran`'s libquadmath to `@loader_path/`, and ad-hoc sign everything with `codesign -f -s -`.

## References

- SuSiE: Wang et al., JRSS-B, doi:10.1111/rssb.12388
- SuSiE-RSS: Zou et al., PLoS Genet, doi:10.1371/journal.pgen.1010299
- FINEMAP: Benner et al., Bioinformatics, doi:10.1093/bioinformatics/btw018
- PolyFun and UK Biobank LD: Weissbrod et al., Nat Genet, doi:10.1038/s41588-020-00735-5
- CARMA: Yang et al., Nat Genet, doi:10.1038/s41588-023-01392-0
- Reference-LD pitfalls: Benner et al., AJHG, doi:10.1016/j.ajhg.2017.08.012; Kanai et al., Cell Genomics, doi:10.1016/j.xgen.2022.100210
- Infinitesimal effects (SuSiE-inf/FINEMAP-inf): Cui et al., Nat Genet, doi:10.1038/s41588-023-01597-3
- Summary statistics: GWAS Catalog, https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/
