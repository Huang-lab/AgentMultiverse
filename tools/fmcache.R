# Shared result cache for slow fine-mapping calls, so identical configurations are computed once
# across all analyses on this machine.
#
#   source("tools/fmcache.R")
#   fit <- cached_call("susieR::susie_rss", list(z = z, R = R, n = n, L = 10))
#   res <- cached_call("CARMA::CARMA", list(list(z), list(R), lambda.list = list(1 / sqrt(p)), outlier.switch = TRUE))
#
# The call is do.call(fun, args), and the cache key covers everything that can change the result:
# the function name, its package version, every argument (numeric vectors and matrices by the md5 of
# their float64 bytes, column-major, plus dimensions; other values verbatim) and `replicate`.
# Because the key is built from the same `args` that are run, a hit can never belong to other inputs.
# Pass replicate = 2, 3, ... to get fresh runs of a stochastic method.
# attr(result, "fmcache") tells whether it was a hit, the key, and the original compute time;
# log it with tools/runlog.py. If another analysis is computing the same key, this waits for it.

.fmcache_root <- local({
  here <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) "tools")
  normalizePath(file.path(here, "..", ".fmcache"), mustWork = FALSE)
})

.fmcache_md5 <- function(x) {
  tmp <- tempfile()
  con <- pipe(paste("openssl dgst -md5 -r >", shQuote(tmp)), "wb")
  n <- length(x)
  step <- 2^26
  for (s in seq(1, max(n, 1), by = step)) {
    if (n == 0) break
    writeBin(as.double(x[s:min(n, s + step - 1)]), con, size = 8, endian = "little")
  }
  close(con)
  h <- strsplit(readLines(tmp, warn = FALSE), " ")[[1]][1]
  unlink(tmp)
  h
}

.fmcache_spec <- function(x) {
  if (is.list(x)) {
    nm <- names(x)
    out <- lapply(x, .fmcache_spec)
    if (!is.null(nm) && all(nzchar(nm))) out <- out[order(nm)]
    return(out)
  }
  if (is.function(x)) stop("fmcache: functions cannot be part of a cache key; pass their result instead")
  if (is.numeric(x) || is.logical(x)) {
    if (length(x) > 1 || !is.null(dim(x))) {
      return(list(md5 = .fmcache_md5(x), length = length(x), dim = if (is.null(dim(x))) NULL else as.list(dim(x))))
    }
    return(x)
  }
  if (is.null(x)) return("NULL")
  x
}

.fmcache_alive <- function(pid) isTRUE(tools::pskill(pid, 0L))

cached_call <- function(fun, args, replicate = 1L) {
  stopifnot(is.character(fun), length(fun) == 1, is.list(args))
  pkg <- if (grepl("::", fun)) sub("::.*", "", fun) else NA
  f <- if (is.na(pkg)) get(fun, mode = "function") else getExportedValue(pkg, sub(".*::", "", fun))
  spec <- list(fun = fun, version = if (is.na(pkg)) "" else as.character(utils::packageVersion(pkg)),
               args = .fmcache_spec(args), replicate = as.integer(replicate))
  json <- as.character(jsonlite::toJSON(spec, auto_unbox = TRUE, digits = NA, null = "null"))
  jf <- tempfile()
  writeLines(json, jf, sep = "", useBytes = TRUE)
  key <- unname(tools::md5sum(jf))
  unlink(jf)
  dir <- file.path(.fmcache_root, gsub("[^A-Za-z0-9_.]", "_", fun), substr(key, 1, 2), key)
  result_file <- file.path(dir, "result.rds")
  lock <- paste0(dir, ".lock")
  dir.create(dirname(dir), recursive = TRUE, showWarnings = FALSE)

  repeat {
    if (file.exists(result_file)) {
      res <- readRDS(result_file)
      meta <- jsonlite::read_json(file.path(dir, "meta.json"))
      attr(res, "fmcache") <- list(hit = TRUE, key = key, seconds = meta$seconds)
      message(sprintf("fmcache: hit %s %s (computed once in %.0f s)", fun, key, meta$seconds))
      return(res)
    }
    if (dir.create(lock, showWarnings = FALSE)) break
    pid <- suppressWarnings(as.integer(readLines(file.path(lock, "pid"), warn = FALSE)))
    if (length(pid) == 1 && !is.na(pid) && !.fmcache_alive(pid)) {
      unlink(lock, recursive = TRUE)
    } else {
      Sys.sleep(2)
    }
  }
  writeLines(as.character(Sys.getpid()), file.path(lock, "pid"))
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)

  t0 <- Sys.time()
  res <- do.call(f, args)
  seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  dir.create(dir, showWarnings = FALSE)
  tmp <- file.path(dir, ".result.rds.part")
  saveRDS(res, tmp)
  jsonlite::write_json(list(key = key, spec = spec, seconds = seconds, created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
                            r_version = R.version.string), file.path(dir, "meta.json"), auto_unbox = TRUE, digits = NA,
                       null = "null")
  file.rename(tmp, result_file)
  attr(res, "fmcache") <- list(hit = FALSE, key = key, seconds = seconds)
  message(sprintf("fmcache: computed %s %s in %.0f s", fun, key, seconds))
  res
}
