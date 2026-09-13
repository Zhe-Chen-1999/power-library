## =====================================================================
## dgp_weibull_bayes.R -- Bayesian Weibull time-to-event design
##
## Ported from the PARMA sample size work (Weibull survival model in JAGS,
## time to resolution of hypoxemia, with death as a competing event).
##
## The only thing that makes this "Bayesian" from the engine's point of
## view is what `reject` means: instead of p < alpha it is
##
##     Pr(beta.trt > 0 | data) >= threshold
##
## i.e. the posterior probability of *any* benefit clears a pre-specified
## bar. The engine averages that over replicates exactly as it averages a
## frequentist rejection, and the average is the Bayesian power.
##
## design_weibull_bayes()  JAGS Weibull PH model
## design_weibull_cox()    same DGP, Cox score test -- the frequentist
##                         comparator, on identical data
## =====================================================================

suppressPackageStartupMessages(library(survival))

## rjags needs JAGS itself installed (brew install jags), which not everyone
## has. Load it lazily so the rest of the library -- and the binary designs
## in particular -- still work on a machine without it. The error then
## arrives when you actually ask for a Bayesian design, and says what to do.
require_rjags <- function() {
  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop("design_weibull_bayes() needs the 'rjags' package and a JAGS install.\n",
         "  macOS:  brew install jags && Rscript -e 'install.packages(\"rjags\")'",
         call. = FALSE)
  }
}

## Absolute path to a bundled JAGS model.
##
## Resolved in the parent process at design construction, not inside the
## worker: parallel workers do not reliably inherit the working directory,
## and a relative path that happens to work when knitting from
## analyses/<name>/ would fail from anywhere else. The absolute string then
## travels to each worker as an ordinary parameter.
jags_model_path <- function(file = "weibull_ph.jags") {
  if (exists("POWERLIB_ROOT", inherits = TRUE)) {
    p <- file.path(get("POWERLIB_ROOT", inherits = TRUE), "jags", file)
    if (file.exists(p)) return(normalizePath(p))
  }
  cands <- file.path(c("jags", "../jags", "../../jags"), file)
  hit <- cands[file.exists(cands)]
  if (!length(hit)) stop("cannot locate ", file, "; pass jags_file = <path>")
  normalizePath(hit[1])
}

## ---------------------------------------------------------------------
## Parameterisation helpers
##
## R's rweibull(shape = a, scale = b) has S(t) = exp(-(t/b)^a), so it maps
## to the JAGS form with alpha = a and lambda = b^(-a). Under proportional
## hazards with log-HR beta, the treated arm therefore needs
##
##     b_trt = b_ctrl * exp(-beta / a)
##
## Getting this exponent wrong is the classic way these simulations end up
## quietly targeting the wrong effect size, so it lives in one function.
## ---------------------------------------------------------------------

weibull_scale_for_arm <- function(scale_ctrl, shape, beta, trt) {
  scale_ctrl * exp(-beta * trt / shape)
}

#' beta0 implied by the data-generating control arm, on the JAGS scale
#'
#' Useful for centring the baseline prior honestly: if you want a prior
#' "centred at a 5-day mean in the control arm", this is the number.
weibull_beta0_true <- function(scale_ctrl, shape) -shape * log(scale_ctrl)

#' Control-arm summaries implied by (shape, scale), for reporting
weibull_control_summary <- function(scale_ctrl, shape) {
  c(mean   = scale_ctrl * gamma(1 + 1 / shape),
    median = scale_ctrl * log(2)^(1 / shape))
}

## ---------------------------------------------------------------------
## Shared DGP: Weibull primary event + competing event + admin censoring
## ---------------------------------------------------------------------

weibull_dgp <- function(params) {
  n   <- params$n
  trt <- rbinom(n, 1, params$p_treat)

  b   <- weibull_scale_for_arm(params$scale_event, params$shape_event,
                               params$beta, trt)
  t_event <- rweibull(n, shape = params$shape_event, scale = b)

  t_comp <- if (identical(params$competing_dist, "none")) {
    rep(Inf, n)
  } else if (identical(params$competing_dist, "weibull")) {
    rweibull(n, shape = params$shape_comp, scale = params$scale_comp)
  } else {
    rexp(n, rate = 1 / params$scale_comp)
  }

  # The competing event and the administrative horizon both act as right
  # censoring for the cause-specific hazard of the primary event.
  t_cen  <- pmin(t_comp, params$horizon)
  status <- as.numeric(t_event <= t_cen)
  t_obs  <- pmin(t_event, t_cen)

  data.frame(id = seq_len(n), trt = trt, time = t_obs,
             status = status, t_cen = t_cen)
}

## ---------------------------------------------------------------------
## Design 1: Bayesian Weibull PH model in JAGS
## ---------------------------------------------------------------------

#' @param n              total randomised
#' @param beta           true log hazard ratio (positive = treated reach the
#'                       event sooner; for a *good* event that is benefit)
#' @param shape_event    Weibull shape for the primary event
#' @param scale_event    Weibull scale, control arm
#' @param competing_dist "exp", "weibull", or "none"
#' @param scale_comp     mean (exp) or scale (weibull) of the competing event
#' @param horizon        administrative censoring time
#' @param threshold      posterior probability bar for declaring success
#' @param prior_*        JAGS hyperparameters; dnorm takes a *precision*
#' @param center_beta0_at_truth  if TRUE, centre the baseline prior on the
#'                       value implied by the DGP instead of a fixed number
## DGP parameters shared by the Bayesian and frequentist versions, so the
## two designs cannot drift apart in what they simulate.
weibull_dgp_defaults <- function() {
  list(
    n = 160, p_treat = 0.5, beta = log(1.6),
    shape_event = 0.9, scale_event = 4,
    competing_dist = "exp", scale_comp = 21, shape_comp = 1,
    horizon = 28
  )
}

design_weibull_bayes <- function() {
  require_rjags()
  new_design(
    name = "weibull_bayes",

    defaults = c(weibull_dgp_defaults(), list(
      threshold = 0.90,
      prior_trt_mean = 0, prior_trt_prec = 0.2,
      prior_beta0_mean = NA, prior_beta0_prec = 0.5,
      center_beta0_at_truth = TRUE,
      prior_alpha_shape = 1.1, prior_alpha_rate = 1.1,
      n_adapt = 500, n_burn = 1000, n_iter = 5000,
      jags_file = jags_model_path("weibull_ph.jags")
    )),

    dgp = weibull_dgp,

    analyze = function(data, params) {
      n <- nrow(data)

      # JAGS wants the event time as NA wherever the observation is
      # censored; dinterval() then samples it above t.cen.
      t <- ifelse(data$status == 1, data$time, NA_real_)

      # horizon = Inf with no competing event gives an infinite censoring
      # time. Nothing is censored in that case, but JAGS will not accept Inf
      # in its data, so substitute a bound that is finite and never binding.
      t_cen <- data$t_cen
      if (any(!is.finite(t_cen))) t_cen[!is.finite(t_cen)] <- 10 * max(data$time)

      t_init <- ifelse(data$status == 0, t_cen + 1, NA_real_)

      # Centre the baseline prior on the value the DGP implies, unless an
      # explicit mean is supplied. Overriding therefore needs *both*
      # center_beta0_at_truth = FALSE and a value for prior_beta0_mean --
      # clearing the flag alone leaves nothing to centre on.
      b0_mean <- if (isTRUE(params$center_beta0_at_truth) || is.na(params$prior_beta0_mean)) {
        weibull_beta0_true(params$scale_event, params$shape_event)
      } else params$prior_beta0_mean

      jdata <- list(
        N = n, t = t, t.cen = t_cen, trt = data$trt,
        is.censored = 1 - data$status,
        prior_beta0_mean  = b0_mean,
        prior_beta0_prec  = params$prior_beta0_prec,
        prior_trt_mean    = params$prior_trt_mean,
        prior_trt_prec    = params$prior_trt_prec,
        prior_alpha_shape = params$prior_alpha_shape,
        prior_alpha_rate  = params$prior_alpha_rate
      )
      jinits <- list(list(t = t_init, beta0 = b0_mean, beta.trt = 0))

      jm <- rjags::jags.model(params$jags_file, data = jdata,
                              inits = jinits, n.chains = 1,
                              n.adapt = params$n_adapt, quiet = TRUE)
      update(jm, params$n_burn, progress.bar = "none")
      s <- rjags::coda.samples(jm, c("beta0", "beta.trt", "alpha"),
                               n.iter = params$n_iter, progress.bar = "none")

      post_beta  <- as.numeric(s[[1]][, "beta.trt"])
      post_alpha <- as.numeric(s[[1]][, "alpha"])

      p_benefit <- mean(post_beta > 0)
      # unname(): quantile() labels its result "2.5%"/"97.5%", and c() would
      # paste those onto the element names, giving columns called
      # `post_lo_beta.2.5%` downstream.
      q <- unname(quantile(post_beta, c(0.025, 0.5, 0.975)))

      c(reject = as.numeric(p_benefit >= params$threshold),
        # Same rule applied against the *assumed true* effect: how often we
        # would claim more benefit than the design actually built in. The
        # PARMA memo reports this as its "type I error"; it is a
        # calibration diagnostic, not a frequentist type I error rate.
        reject_vs_truth = as.numeric(mean(post_beta > params$beta) >= params$threshold),
        # Frequentist-style error control on the same posterior: 95%
        # credible interval excluding the null.
        reject_ci95     = as.numeric(q[1] > 0),
        prob_benefit    = p_benefit,
        post_mean_beta  = mean(post_beta),
        post_lo_beta    = q[1],
        post_hi_beta    = q[3],
        post_mean_HR    = mean(exp(post_beta)),
        post_mean_alpha = mean(post_alpha),
        n_events        = sum(data$status),
        event_frac      = mean(data$status))
    },

    meta = list(
      design = "Parallel-arm RCT, Bayesian analysis",
      outcome = "Time-to-event (Weibull PH)",
      features = "Competing risk as cause-specific censoring; admin censoring; posterior-probability decision rule; prior sensitivity"
    )
  )
}

## ---------------------------------------------------------------------
## Design 2: frequentist comparator on the identical DGP
##
## Sharing weibull_dgp() means any difference in power is attributable to
## the analysis, not to a differently simulated trial.
## ---------------------------------------------------------------------

design_weibull_cox <- function() {
  new_design(
    name = "weibull_cox",
    defaults = c(weibull_dgp_defaults(), list(alpha = 0.05)),
    dgp = weibull_dgp,
    analyze = function(data, params) {
      fit <- survival::coxph(survival::Surv(time, status) ~ trt, data = data)
      sm  <- summary(fit)
      c(reject   = as.numeric(sm$coefficients["trt", "Pr(>|z|)"] < params$alpha),
        loghr    = unname(sm$coefficients["trt", "coef"]),
        se_loghr = unname(sm$coefficients["trt", "se(coef)"]),
        n_events = sum(data$status))
    },
    meta = list(
      design = "Parallel-arm RCT, frequentist comparator",
      outcome = "Time-to-event (Cox PH)",
      features = "Same DGP as weibull_bayes; isolates the effect of the analysis model"
    )
  )
}
