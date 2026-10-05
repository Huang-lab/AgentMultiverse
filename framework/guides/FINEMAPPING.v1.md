# Fine-mapping toolkit

Local software and reference data for fine-mapping GWAS loci from summary statistics.
Everything lives in this folder, and all paths below are relative to the project root.
Built and verified on macOS arm64 on 2026-09-26; the Rebuild section recreates it on any machine.

```bash
export PATH="$PWD/bin:$PWD/env/bin:$PATH"   # python, Rscript, plink, plink2, finemap, liftOver, bcftools, tabix
```

## What is here

| Path | Contents |
|---|---|
| `bin/plink2`, `bin/plink` | PLINK 2.0 alpha 7.8 (arm64) and PLINK 1.9 (universal) |
| `bin/finemap` + `bin/lib/` | FINEMAP 1.4.2, macOS x86_64 build run under Rosetta; its GCC, OpenMP and zstd libraries are bundled in `bin/lib`, so keep the two together |
| `bin/liftOver` | UCSC liftOver (arm64) |
| `env/` | conda env: Python 3.12 (numpy, scipy, pandas 2, pyarrow, scikit-learn, pandas-plink, matplotlib), R 4.5 (susieR 0.14.2, CARMA 1.0, data.table with R.utils so `fread` reads `.gz`, dplyr), bcftools and htslib 1.24 |
| `tools/polyfun/` | PolyFun: `snpvar_meta.chr*.parquet` (precomputed functional priors, hg19, ~19M UKB SNPs), `ukb_regions.tsv.gz` (index of UK Biobank LD windows), `finemapper.py` (`load_ld_npz`, SuSiE/FINEMAP wrapper) |
| `tools/runlog.py` | Run record: `runlog.py <run_dir> {init,action,decision,issue,result,output} "..."` appends to `LOG.md`, `log.jsonl` and `manifest.tsv` (outputs with md5); a decision needs `--why` |
| `ref/1kg/hg19/` | 1000 Genomes phase 3, GRCh37, 2504 samples: `all_phase3.{pgen,pvar.zst,psam}`, IDs are rsIDs |
| `ref/1kg/hg38/` | 1000 Genomes NYGC 30x, GRCh38, 3202 samples including relatives: `all_hg38.{pgen,pvar.zst,psam}`, IDs are dbSNP 156 rsIDs |
| `ref/1kg/*/{EUR,AFR,EAS,SAS,AMR}.keep` | Unrelated samples per superpopulation (hg19: 503/652/504/484/347, hg38: 525/680/506/511/351) |
| `ref/ukbb_ld/` | Cache for UK Biobank LD windows, empty until needed |
| `ref/ld_blocks/` | Approximately independent LD blocks: hg19 Berisa-Pickrell (EUR/AFR/ASN), hg38 pyrho (EUR/AFR/EAS/SAS) |
| `ref/liftover/` | UCSC chains `hg19ToHg38`, `hg38ToHg19` |
| `ref/gencode/` | GENCODE v50 basic GTF for GRCh38 and GRCh37 (`lift37`) |
| `data/sumstats/`, `results/` | Inputs and outputs |

UK Biobank LD (337K British-ancestry samples, GRCh37, 3 Mb windows every 1 Mb) is the best reference for European hg19 summary stats.
A window is typically about 1 GB (up to about 3 GB), so fetch it per locus rather than in bulk.
Pick the row of `tools/polyfun/ukb_regions.tsv.gz` whose window contains the locus near its center, download `<URL_PREFIX>.npz` and `<URL_PREFIX>.gz` into `ref/ukbb_ld/`, and load them with `load_ld_npz("ref/ukbb_ld/<name>")` from `tools/polyfun/finemapper.py`.
That returns the full matrix plus a SNP table, with signs relative to `allele1`.

## Recipe

1. Harmonize the sumstats to CHR, BP, effect and other allele, BETA and SE (or Z), P, per-SNP N, and EAF.
   Detect the build by looking up a few rsIDs in the `.pvar` files (rs12740374 is chr1:109817590 in hg19 and chr1:109274968 in hg38), and liftOver if needed.
2. Define loci as lead variants with P < 5e-8 plus ±0.5-1 Mb (merge overlaps), or use the LD blocks; flag the MHC (chr6:25-34 Mb).
3. Get LD: use in-sample LD if available, else UK Biobank LD for EUR hg19, else the 1000G superpopulation that matches the GWAS.
   ```bash
   plink2 --pfile ref/1kg/hg19/all_phase3 vzs --keep ref/1kg/hg19/EUR.keep --chr 1 --from-bp 109317590 --to-bp 110317590 \
          --maf 0.01 --max-alleles 2 --rm-dup exclude-all --make-pgen --out results/locus      # ~3 s
   plink2 --pfile results/locus --r-unphased square bin4 ref-based --out results/locus         # float32 m x m, rows = results/locus.pvar
   ```
   Load it with `np.fromfile(f, np.float32).reshape(m, m)` in Python or `matrix(readBin(f, "double", m*m, size = 4), m, m)` in R.
4. Align every z to the reference's counted allele (1000G `ALT`, UK Biobank `allele1`), flip the sign when alleles are swapped, and drop mismatches, duplicates and strand-ambiguous A/T and C/G SNPs.
5. Run QC before trusting any result: `susieR::estimate_s_rss` (near 0 means z and LD agree), `kriging_rss` outliers, and convergence.
6. Fine-map with `susie_rss(z = z, R = R, n = n, L = 10)` and take 95% credible sets from `susie_get_cs(fit, Xcorr = R, min_abs_corr = 0.5)`.
   Cross-check with FINEMAP (`finemap --sss`, `.z` columns `rsid chromosome position allele1 allele2 maf beta se` with beta = z and se = 1) and with CARMA (`CARMA(list(z), list(R), lambda.list = list(1/sqrt(p)), outlier.switch = TRUE)`), which models reference-LD mismatch.
   For functional priors (hg19, EUR), pass PolyFun `snpvar_bin` as `prior_weights`.
7. Report per locus: PIPs, credible sets (size, purity, top variant), `estimate_s_rss`, convergence, agreement between methods, and nearby genes from GENCODE.
   Record actions, decisions and outputs with `tools/runlog.py` as you go.

Pitfalls found while testing this setup:

- `plink2 --r-unphased` signs are relative to the major allele unless you pass `ref-based`; without it the signs will not match the z orientation.
- Meta-analyses have variable per-SNP N, and mixing SNPs with very different N breaks z/LD consistency.
  On the SORT1 test, keeping N ≥ 0.9·max(N) moved `estimate_s_rss` from 0.51 to 0.06 and made SuSiE converge.
- 1000G has only 347-680 unrelated samples per superpopulation, so with very strong signals (|z| > ~30) expect spurious extra signals (FINEMAP piling up at the maximum number of causal SNPs, many singleton SuSiE sets).
  On the SORT1 test, 1000G EUR gave 6 SuSiE credible sets while UK Biobank LD gave 1, so prefer UK Biobank LD for EUR, or in-sample LD.
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
  r-rcpp r-rcpparmadillo r-rcppgsl r-dplyr r-glmnet r-mass gsl bcftools htslib compilers
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
