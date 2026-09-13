## =====================================================================
## engine.R -- the simulation core shared by every design in the library
##
## The contract is deliberately small. A *design* is three things:
##
##   dgp(params)          -> a data.frame: one simulated trial
##   analyze(data, params)-> a named numeric vector: one trial's results,
##                           which MUST contain an element called `reject`
##                           (1 = "declare success", 0 = not)
##   defaults             -> a list of parameters the design assumes unless
##                           a scenario row overrides them
##
## Whatever `reject` means is the design's business: p < alpha for a
## frequentist test, Pr(HR > 1) > 0.90 for a Bayesian one. The engine only
## averages it. Everything else analyze() returns is averaged too and
## carried along, which is how we get posterior means, event counts,
## convergence rates etc. for free.
##
## Adding a new design = write those three things. Nothing here changes.
## =====================================================================

suppressPackageStartupMessages({
  library(future.apply)
  library(tibble)
  library(dplyr)
})

## ---------------------------------------------------------------------
## Design constructor
## ---------------------------------------------------------------------

#' Define a power-simulation design
#'
#' @param name      short identifier, used in output and the registry
#' @param dgp       function(params) -> data.frame
#' @param analyze   function(data, params) -> named numeric vector with `reject`
#' @param defaults  named list of default parameters
#' @param analytic  optional function(params) -> closed-form power, used by
#'                  check_against_analytic() to validate the simulation
#' @param meta      named list of registry metadata (design, outcome, etc.)
new_design <- function(name, dgp, analyze, defaults = list(),
                       analytic = NULL, meta = list()) {
  stopifnot(is.function(dgp), is.function(analyze))
  structure(
    list(name = name, dgp = dgp, analyze = analyze, defaults = defaults,
         analytic = analytic, meta = meta),
    class = "power_design"
  )
}

print.power_design <- function(x, ...) {
  cat("<power_design>", x$name, "\n")
  cat("  defaults:", paste(names(x$defaults), collapse = ", "), "\n")
  cat("  analytic check:", if (is.null(x$analytic)) "no" else "yes", "\n")
  invisible(x)
}

## ---------------------------------------------------------------------
## Parameter resolution
##
## A scenario grid is a data.frame with one row per scenario. Its columns
## override design defaults. List-columns are supported so that a scenario
## can carry a vector-valued parameter (e.g. p = c(.40, .55, .55) for a
## three-arm trial) in a single cell.
## ---------------------------------------------------------------------

resolve_params <- function(design, grid_row) {
  params <- design$defaults
  for (nm in names(grid_row)) {
    v <- grid_row[[nm]]
    if (is.list(v)) v <- v[[1]]          # unwrap a list-column cell
    params[[nm]] <- v
  }
  params
}

#' A design's defaults with some parameters overridden
#'
#' Use this rather than `c(design$defaults, list(...))`: `c()` on two lists
#' with a shared name keeps *both*, and `params$x` then silently returns the
#' default instead of your override.
design_params <- function(design, ...) {
  modifyList(design$defaults, list(...))
}

#' Build a scenario grid, keeping vector-valued parameters intact
#'
#' Plain values are crossed as in expand.grid(). To vary a vector-valued
#' parameter, pass a list of vectors -- each element becomes one scenario
#' level and is stored in a list-column.
scenarios <- function(...) {
  args <- list(...)
  idx <- do.call(expand.grid, c(lapply(args, seq_along),
                                list(KEEP.OUT.ATTRS = FALSE)))
  out <- lapply(names(args), function(nm) {
    vals <- args[[nm]]
    if (is.list(vals)) I(vals[idx[[nm]]]) else vals[idx[[nm]]]
  })
  names(out) <- names(args)
  tibble::as_tibble(out)
}

## ---------------------------------------------------------------------
## Single replicate (also the debugging entry point)
## ---------------------------------------------------------------------

#' Run one replicate. Useful on its own when a scenario misbehaves and you
#' want to see the simulated data rather than an averaged number.
replicate_one <- function(design, params, seed = NULL, return_data = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  dat <- design$dgp(params)
  res <- design$analyze(dat, params)
  if (return_data) list(data = dat, result = res) else res
}

## ---------------------------------------------------------------------
## The engine
## ---------------------------------------------------------------------

#' Run a power simulation over a grid of scenarios
#'
#' @param design   a power_design
#' @param grid     data.frame of scenarios; NULL runs the defaults alone
#' @param n_sim    replicates per scenario
#' @param seed     integer; fixes the parallel RNG streams, so results are
#'                 reproducible regardless of how many workers are used
#' @param workers  number of parallel workers; 1 forces sequential
#' @param verbose  print per-scenario progress
#'
#' @return tibble: one row per scenario, with `power`, its Monte Carlo SE,
#'   a 95% MC interval, the count of failed replicates, and the mean of
#'   every other quantity analyze() returned.
run_power <- function(design, grid = NULL, n_sim = 1000,
                      seed = 20260913, workers = max(1, parallel::detectCores() - 1),
                      verbose = TRUE) {

  stopifnot(inherits(design, "power_design"))
  if (is.null(grid)) grid <- tibble::tibble(.scenario = 1L)
  grid <- tibble::as_tibble(grid)
  n_scen <- nrow(grid)

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  if (workers > 1) future::plan(future::multisession, workers = workers)
  else future::plan(future::sequential)

  if (verbose) {
    message(sprintf("[%s] %d scenario%s x %d replicates on %d worker%s",
                    design$name, n_scen, if (n_scen == 1) "" else "s",
                    n_sim, workers, if (workers == 1) "" else "s"))
  }

  t0 <- Sys.time()
  rows <- vector("list", n_scen)

  for (i in seq_len(n_scen)) {
    params <- resolve_params(design, grid[i, , drop = FALSE])
    dgp <- design$dgp; analyze <- design$analyze

    reps <- future.apply::future_lapply(
      seq_len(n_sim),
      function(r) {
        tryCatch({
          d <- dgp(params)
          out <- analyze(d, params)
          if (!"reject" %in% names(out)) {
            stop("analyze() must return an element named 'reject'")
          }
          out
        }, error = function(e) {
          structure(NA_real_, names = "reject", failed = conditionMessage(e))
        })
      },
      # One stream per scenario, derived from `seed`, so scenario i always
      # sees the same draws no matter what else is in the grid.
      future.seed = seed + 1000L * i
    )

    rows[[i]] <- summarise_reps(reps, n_sim)
    if (verbose) {
      message(sprintf("  scenario %d/%d: power = %.3f (SE %.3f)%s",
                      i, n_scen, rows[[i]]$power, rows[[i]]$mc_se,
                      if (rows[[i]]$n_failed > 0)
                        sprintf("  [%d failed]", rows[[i]]$n_failed) else ""))
    }
  }

  res <- dplyr::bind_rows(rows)
  out <- dplyr::bind_cols(grid, res)
  attr(out, "design")  <- design$name
  attr(out, "n_sim")   <- n_sim
  attr(out, "seed")    <- seed
  attr(out, "runtime") <- difftime(Sys.time(), t0, units = "secs")
  class(out) <- c("power_result", class(out))
  out
}

#' Wilson interval for the simulated power
#'
#' Used instead of the Wald interval because power estimates near 0 or 1 are
#' exactly where check_against_analytic() gets used, and there the Wald
#' interval collapses to a point and would reject a correct design.
wilson_ci <- function(p, n, z = 1.96) {
  if (!is.finite(p) || n <= 0) return(c(NA_real_, NA_real_))
  d      <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / d
  half   <- z / d * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  c(max(0, centre - half), min(1, centre + half))
}

## Collapse the replicate list into one summary row.
summarise_reps <- function(reps, n_sim) {
  # A replicate counts as failed if it errored *or* if it came back with a
  # non-finite `reject`. Averaging with na.rm would otherwise quietly drop
  # non-converged fits, which flatters exactly the scenarios that struggle.
  ok <- vapply(reps, function(x) {
    r <- suppressWarnings(as.numeric(x["reject"]))
    length(r) == 1L && !is.na(r) && is.finite(r)
  }, logical(1))
  n_failed <- sum(!ok)
  if (!any(ok)) {
    return(tibble::tibble(power = NA_real_, mc_se = NA_real_,
                          power_lo = NA_real_, power_hi = NA_real_,
                          n_sim = n_sim, n_failed = n_failed))
  }

  keep <- reps[ok]
  nms <- unique(unlist(lapply(keep, names)))
  mat <- vapply(keep, function(x) {
    v <- rep(NA_real_, length(nms)); names(v) <- nms
    v[names(x)] <- as.numeric(x)
    v
  }, numeric(length(nms)))
  if (length(nms) == 1L) mat <- matrix(mat, nrow = 1, dimnames = list(nms, NULL))

  means <- rowMeans(mat, na.rm = TRUE)
  power <- unname(means["reject"])
  n_ok  <- sum(ok)
  # MC SE of a proportion; the engine's own precision, not the trial's.
  mc_se <- sqrt(power * (1 - power) / n_ok)
  ci    <- wilson_ci(power, n_ok)

  extra <- means[setdiff(nms, "reject")]
  out <- tibble::tibble(
    power    = power,
    mc_se    = mc_se,
    power_lo = ci[1],
    power_hi = ci[2],
    n_sim    = n_sim,
    n_failed = n_failed
  )
  if (length(extra)) out <- dplyr::bind_cols(out, tibble::as_tibble(as.list(extra)))
  out
}

print.power_result <- function(x, ...) {
  cat(sprintf("<power_result> %s | %d sim/scenario | seed %s | %.1f s\n",
              attr(x, "design"), attr(x, "n_sim"), attr(x, "seed"),
              as.numeric(attr(x, "runtime"))))
  NextMethod()
}

## ---------------------------------------------------------------------
## Validation helpers
## ---------------------------------------------------------------------

#' Compare simulated power against a design's closed-form power
#'
#' Every simulation in this library that *has* an analytic counterpart
#' should be checked against it once. It is the cheapest available guard
#' against a data-generating bug, and it is how we know the engine itself
#' is wired correctly.
#'
#' @return the result tibble plus `power_analytic`, `abs_diff`, and
#'   `within_mc` (is the analytic value inside the MC interval?)
check_against_analytic <- function(design, grid = NULL, n_sim = 2000, ...) {
  if (is.null(design$analytic)) stop("design '", design$name, "' has no analytic form")
  res <- run_power(design, grid, n_sim = n_sim, ...)
  if (is.null(grid)) grid <- tibble::tibble(.scenario = 1L)
  ana <- vapply(seq_len(nrow(grid)), function(i) {
    design$analytic(resolve_params(design, grid[i, , drop = FALSE]))
  }, numeric(1))
  res$power_analytic <- ana
  res$abs_diff  <- abs(res$power - ana)
  res$within_mc <- ana >= res$power_lo & ana <= res$power_hi
  res
}

#' Smallest n per arm reaching `target` power, by bisection on the simulation
#'
#' Simulation-based sample size search. Each candidate n costs a full
#' `n_sim` run, so keep n_sim modest while searching and confirm the answer
#' with one high-precision run at the end.
solve_n <- function(design, params = list(), target = 0.80,
                    n_range = c(20, 5000), n_sim = 500, tol = 5,
                    n_arg = "n_per_arm", verbose = TRUE, ...) {
  lo <- n_range[1]; hi <- n_range[2]
  power_at <- function(n) {
    cells <- list(n)
    names(cells) <- n_arg
    for (nm in names(params)) {
      cells[[nm]] <- if (length(params[[nm]]) > 1) I(list(params[[nm]])) else params[[nm]]
    }
    p <- run_power(design, tibble::as_tibble(cells), n_sim = n_sim,
                   verbose = FALSE, ...)$power
    if (verbose) message(sprintf("  n = %d -> power %.3f", n, p))
    p
  }
  if (power_at(hi) < target) {
    warning("target power not reached at n = ", hi)
    return(NA_integer_)
  }
  while (hi - lo > tol) {
    mid <- floor((lo + hi) / 2)
    if (power_at(mid) >= target) hi <- mid else lo <- mid
  }
  hi
}
