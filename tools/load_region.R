# Load one region of a shared-prep study in R: z aligned to the panel's a1, the LD submatrix, and a table.
# Wraps tools/load_region.py (see its header for the exact filters); every choice is a required argument.
#
#   source("tools/load_region.R")
#   d <- load_region("prep/eadb2026", "GCST90704647", "chr19_...", "ukb_hg19",
#                    ambiguous = "drop", min_n_frac = 0.9, min_maf = 0)
#   d$z, d$R (m x m), d$table (data.table; `row` indexes the panel matrix), d$counts

load_region <- function(prep, sumstats, region, panel, ambiguous, min_n_frac, min_maf) {
  root <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "tools/x"), ".."), mustWork = FALSE)
  out <- tempfile("region_")
  args <- c(file.path(root, "tools/load_region.py"), prep, sumstats, region, panel, "--ambiguous", ambiguous,
            "--min-n-frac", format(min_n_frac, digits = 17), "--min-maf", format(min_maf, digits = 17), "--out", out)
  status <- system2(file.path(root, "env/bin/python"), shQuote(args))
  if (status != 0) stop("load_region.py failed")
  info <- jsonlite::read_json(paste0(out, ".counts.json"))
  table <- data.table::fread(paste0(out, ".tsv"))
  m <- info$m
  R <- matrix(readBin(paste0(out, ".ld.bin"), "double", m * m, size = 4), m, m)
  unlink(paste0(out, c(".tsv", ".ld.bin", ".counts.json")))
  list(z = table$z_aligned, R = R, table = table, counts = info$counts, choices = info$choices)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
