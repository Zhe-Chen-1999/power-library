## =====================================================================
## load.R -- single entry point for the library
##
##   source("R/load.R")           # from the library root
##   source("../../R/load.R")     # from analyses/<name>/
##
## Sources the engine and every design module, and makes the library root
## available as POWERLIB_ROOT so analyses can find jags/ and registry/
## regardless of their own working directory.
## =====================================================================

## Candidate starting points, most reliable first:
##   1. the directory of this file, if R told us (source(), knitr)
##   2. the script directory, if we were launched by Rscript
##   3. the working directory
## From each we walk up looking for R/engine.R, so the library loads the
## same way from the root, from analyses/<name>/, or from tests/.
.powerlib_find_root <- function() {
  starts <- character()

  for (i in seq_len(sys.nframe())) {
    f <- sys.frame(i)$ofile
    if (!is.null(f) && is.character(f)) {
      starts <- c(starts, dirname(normalizePath(f, mustWork = FALSE)))
    }
  }

  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)
  if (length(fa)) {
    starts <- c(starts, dirname(normalizePath(sub("^--file=", "", fa[1]), mustWork = FALSE)))
  }
  starts <- c(starts, getwd())

  for (s in starts) {
    for (up in c(".", "..", "../..", "../../..")) {
      cand <- normalizePath(file.path(s, up), mustWork = FALSE)
      if (file.exists(file.path(cand, "R", "engine.R"))) return(cand)
    }
  }
  stop("could not locate the power-library root (no R/engine.R found at or above ",
       paste(unique(starts), collapse = ", "), ")")
}

POWERLIB_ROOT <- .powerlib_find_root()

for (f in c("engine.R", "dgp_binary.R", "dgp_weibull_bayes.R", "registry.R")) {
  source(file.path(POWERLIB_ROOT, "R", f))
}

#' Path to a file inside the library, from anywhere
powerlib_path <- function(...) file.path(POWERLIB_ROOT, ...)

## Designs available in this library, keyed by the name used in the registry.
POWERLIB_DESIGNS <- list(
  binary_parallel = design_binary_parallel,
  binary_cluster  = design_binary_cluster,
  weibull_bayes   = design_weibull_bayes,
  weibull_cox     = design_weibull_cox
)

#' List the designs the library currently provides
list_designs <- function() {
  do.call(rbind, lapply(names(POWERLIB_DESIGNS), function(nm) {
    m <- POWERLIB_DESIGNS[[nm]]()$meta
    data.frame(design_fn = nm, design = m$design %||% "",
               outcome = m$outcome %||% "", features = m$features %||% "")
  }))
}

`%||%` <- function(a, b) if (is.null(a)) b else a
