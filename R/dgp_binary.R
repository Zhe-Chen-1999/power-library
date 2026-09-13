## =====================================================================
## dgp_binary.R -- binary-outcome designs
##
##   design_binary_parallel()  individually randomised, 2+ arms
##   design_binary_cluster()   cluster randomised, 2+ arms, ICC
##
## Both take a vector `p` of arm-specific event probabilities, so the arms
## are never assumed symmetric -- a three-arm trial with two *different*
## active arms is just p = c(0.40, 0.55, 0.50).
##
## Multiplicity is a parameter, not a hard-coded constant: pass
## `n_comparisons` and `adjust = "bonferroni"` and the engine reports power
## at the adjusted level.
## =====================================================================

suppressPackageStartupMessages({
  library(lme4)
})

## ---------------------------------------------------------------------
## shared helpers
## ---------------------------------------------------------------------

#' Testing level after multiplicity adjustment
alpha_adjusted <- function(alpha, adjust = c("none", "bonferroni"),
                           n_comparisons = 1) {
  adjust <- match.arg(adjust)
  switch(adjust, none = alpha, bonferroni = alpha / n_comparisons)
}

#' Beta parameters matching a target mean and ICC
#'
#' For clustered binary data the ICC on the proportion scale is
#' 1 / (1 + a + b), so a + b = (1 - rho) / rho pins the ICC exactly and the
#' mean pins the split. This gives the intended ICC by construction rather
#' than by an approximation, which matters when rho is large.
beta_from_icc <- function(p, icc) {
  if (icc <= 0) return(NULL)
  s <- (1 - icc) / icc
  list(a = p * s, b = (1 - p) * s)
}

#' Two-proportion test on a simulated arm pair
#'
#' `test = "chisq"` is the uncorrected two-sample z/chi-square test, which
#' is what power.prop.test() assumes -- keeping them aligned is what lets
#' the analytic check be meaningful.
test_two_props <- function(y1, y0, test = c("chisq", "fisher", "logistic")) {
  test <- match.arg(test)
  x <- c(sum(y1), sum(y0)); n <- c(length(y1), length(y0))
  pv <- switch(test,
    chisq    = suppressWarnings(prop.test(x, n, correct = FALSE)$p.value),
    fisher   = fisher.test(matrix(c(x[1], n[1] - x[1], x[2], n[2] - x[2]), 2))$p.value,
    logistic = {
      d <- data.frame(y = c(y1, y0), g = rep(c(1, 0), n))
      suppressWarnings(coef(summary(glm(y ~ g, binomial, d)))["g", "Pr(>|z|)"])
    })
  # Degenerate tables (no events, or no non-events, in either arm) give a
  # 0/0 statistic. There is no evidence against the null there, so treat it
  # as a non-rejection rather than letting NaN propagate into the average.
  if (!is.finite(pv)) 1 else pv
}

#' Exact power of the uncorrected two-proportion test, by enumeration
#'
#' Sums the joint binomial probability over every 2x2 table the trial could
#' produce and adds up the ones that reject. No normal approximation, no
#' Monte Carlo error -- this is the true power of the test that
#' `test_two_props(test = "chisq")` actually performs, which is what makes
#' it the right yardstick for validating the simulation.
#'
#' Falls back to the normal approximation above `max_n`, where the
#' enumeration grid gets expensive and the two agree to ~1e-3 anyway.
power_two_props_exact <- function(n1, n0, p1, p0, alpha, max_n = 2500) {
  if (max(n1, n0) > max_n) return(power_two_props_normal(n1, n0, p1, p0, alpha))
  x1 <- 0:n1; x0 <- 0:n0
  ph1 <- x1 / n1; ph0 <- x0 / n0
  pbar <- outer(x1, x0, "+") / (n1 + n0)
  se   <- sqrt(pbar * (1 - pbar) * (1 / n1 + 1 / n0))
  z    <- (outer(ph1, ph0, "-")) / se
  rej  <- is.finite(z) & abs(z) > qnorm(1 - alpha / 2)
  sum(outer(dbinom(x1, n1, p1), dbinom(x0, n0, p0))[rej])
}

#' Normal-approximation power, i.e. what power.prop.test() reports
#'
#' Kept because it is the convention in protocols and grant applications.
#' It is mildly conservative relative to the exact power of the test.
power_two_props_normal <- function(n1, n0, p1, p0, alpha) {
  n <- harmonic_n(c(n1, n0))
  power.prop.test(n = n, p1 = p0, p2 = p1, sig.level = alpha,
                  alternative = "two.sided")$power
}

## ---------------------------------------------------------------------
## Design 1: individually randomised parallel-arm trial, binary outcome
## ---------------------------------------------------------------------

#' @param p             vector of event probabilities, one per arm
#' @param n_per_arm     scalar (equal allocation) or vector of arm sizes
#' @param contrast      c(treat_index, control_index) into `p`
#' @param alpha         family-wise level before adjustment
#' @param adjust        "none" or "bonferroni"
#' @param n_comparisons number of comparisons the adjustment protects
#' @param test          "chisq", "fisher", or "logistic"
design_binary_parallel <- function() {
  new_design(
    name = "binary_parallel",

    defaults = list(
      p = c(0.40, 0.55), n_per_arm = 200, contrast = c(2, 1),
      alpha = 0.05, adjust = "none", n_comparisons = 1, test = "chisq",
      analytic_method = "exact"
    ),

    dgp = function(params) {
      k <- length(params$p)
      n <- params$n_per_arm
      if (length(n) == 1) n <- rep(n, k)
      data.frame(
        arm = rep(seq_len(k), times = n),
        y   = unlist(lapply(seq_len(k), function(j) rbinom(n[j], 1, params$p[j])))
      )
    },

    analyze = function(data, params) {
      a  <- alpha_adjusted(params$alpha, params$adjust, params$n_comparisons)
      i1 <- params$contrast[1]; i0 <- params$contrast[2]
      y1 <- data$y[data$arm == i1]; y0 <- data$y[data$arm == i0]
      pv <- test_two_props(y1, y0, params$test)
      c(reject     = as.numeric(pv < a),
        risk_diff  = mean(y1) - mean(y0),
        p_treat    = mean(y1),
        p_control  = mean(y0),
        alpha_used = a)
    },

    # Default is the exact enumerated power of the test that analyze()
    # actually runs, so check_against_analytic() is a genuine correctness
    # check. Set analytic_method = "normal" to reproduce power.prop.test()
    # instead, which is the number protocols usually quote.
    analytic = function(params) {
      a <- alpha_adjusted(params$alpha, params$adjust, params$n_comparisons)
      n <- params$n_per_arm
      if (length(n) == 1) n <- rep(n, length(params$p))
      n1 <- n[params$contrast[1]]; n0 <- n[params$contrast[2]]
      p1 <- params$p[params$contrast[1]]; p0 <- params$p[params$contrast[2]]
      if (identical(params$analytic_method, "normal")) {
        power_two_props_normal(n1, n0, p1, p0, a)
      } else {
        power_two_props_exact(n1, n0, p1, p0, a)
      }
    },

    meta = list(
      design = "Parallel-arm RCT (individually randomised)",
      outcome = "Binary",
      features = "2+ arms; unequal arm effects; Bonferroni; subgroup power via reduced n"
    )
  )
}

## power.prop.test takes a single n; for unequal arms use the harmonic mean,
## which is the usual equal-information equivalent.
harmonic_n <- function(n) 2 / (1 / n[1] + 1 / n[2])

## ---------------------------------------------------------------------
## Design 2: cluster-randomised trial, binary outcome
## ---------------------------------------------------------------------

#' @param p            vector of event probabilities, one per arm
#' @param n_clusters   clusters per arm (scalar or vector)
#' @param m            cluster size; scalar, or c(mean, cv) when
#'                     `m_varies = TRUE` for unequal cluster sizes
#' @param icc          intracluster correlation on the proportion scale
#' @param test         "cluster_ttest" (t-test on cluster proportions, the
#'                     standard small-k CRT analysis) or "glmm"
design_binary_cluster <- function() {
  new_design(
    name = "binary_cluster",

    defaults = list(
      p = c(0.40, 0.55), n_clusters = 20, m = 30, icc = 0.05,
      m_varies = FALSE, m_cv = 0.5, contrast = c(2, 1),
      alpha = 0.05, adjust = "none", n_comparisons = 1, test = "cluster_ttest"
    ),

    dgp = function(params) {
      k  <- length(params$p)
      nc <- params$n_clusters; if (length(nc) == 1) nc <- rep(nc, k)
      out <- vector("list", k)
      cid <- 0L
      for (j in seq_len(k)) {
        ab <- beta_from_icc(params$p[j], params$icc)
        # Cluster-level risks; icc = 0 collapses to the constant p.
        pj <- if (is.null(ab)) rep(params$p[j], nc[j]) else rbeta(nc[j], ab$a, ab$b)
        sizes <- if (isTRUE(params$m_varies)) {
          # Gamma sizes with the requested CV, floored at 2 so a cluster
          # always contributes a proportion.
          pmax(2, round(rgamma(nc[j], shape = 1 / params$m_cv^2,
                               scale = params$m * params$m_cv^2)))
        } else rep(params$m, nc[j])
        out[[j]] <- data.frame(
          arm     = j,
          cluster = rep(cid + seq_len(nc[j]), times = sizes),
          y       = unlist(mapply(function(pp, ss) rbinom(ss, 1, pp),
                                  pj, sizes, SIMPLIFY = FALSE))
        )
        cid <- cid + nc[j]
      }
      do.call(rbind, out)
    },

    analyze = function(data, params) {
      a  <- alpha_adjusted(params$alpha, params$adjust, params$n_comparisons)
      i1 <- params$contrast[1]; i0 <- params$contrast[2]
      d  <- data[data$arm %in% c(i1, i0), ]

      if (params$test == "cluster_ttest") {
        agg <- aggregate(y ~ cluster + arm, d, mean)
        pv  <- t.test(agg$y[agg$arm == i1], agg$y[agg$arm == i0])$p.value
      } else {
        d$g <- as.numeric(d$arm == i1)
        fit <- suppressWarnings(suppressMessages(
          lme4::glmer(y ~ g + (1 | cluster), data = d, family = binomial,
                      control = lme4::glmerControl(calc.derivs = FALSE))))
        pv <- coef(summary(fit))["g", "Pr(>|z|)"]
      }

      y1 <- d$y[d$arm == i1]; y0 <- d$y[d$arm == i0]
      c(reject     = as.numeric(pv < a),
        risk_diff  = mean(y1) - mean(y0),
        n_total    = nrow(d),
        alpha_used = a)
    },

    # Standard design-effect inflation: n_eff = n / (1 + (m - 1) * icc).
    # Exact only for equal cluster sizes, so treat it as a sanity check
    # rather than ground truth when m_varies = TRUE.
    analytic = function(params) {
      a  <- alpha_adjusted(params$alpha, params$adjust, params$n_comparisons)
      nc <- params$n_clusters; if (length(nc) > 1) nc <- harmonic_n(nc[params$contrast])
      deff <- 1 + (params$m - 1) * params$icc
      power.prop.test(n = nc * params$m / deff,
                      p1 = params$p[params$contrast[2]],
                      p2 = params$p[params$contrast[1]],
                      sig.level = a, alternative = "two.sided")$power
    },

    meta = list(
      design = "Cluster-randomised trial",
      outcome = "Binary",
      features = "ICC (exact, beta-binomial); unequal cluster sizes; cluster-level t-test or GLMM"
    )
  )
}
