## =====================================================================
## test_engine.R -- run with:  Rscript tests/test_engine.R
##
## These are the checks that make the library trustworthy to reuse. The
## important ones compare simulated power against an exactly computed
## value: if the engine or a DGP is wrong, that is where it shows up.
## =====================================================================

## Works from the library root (Rscript tests/test_engine.R) or from tests/.
source(if (file.exists("R/load.R")) "R/load.R" else "../R/load.R")

pass <- 0L; fail <- 0L
ok <- function(label, cond, detail = "") {
  if (isTRUE(cond)) { pass <<- pass + 1L; cat(sprintf("  PASS  %s\n", label)) }
  else { fail <<- fail + 1L; cat(sprintf("  FAIL  %s  %s\n", label, detail)) }
}

cat("\n== 1. binary_parallel: simulated power vs EXACT enumerated power ==\n")
cat("   (exact = the true power of the test analyze() runs, so the only\n")
cat("    expected discrepancy is Monte Carlo error)\n")
d <- design_binary_parallel()
g <- scenarios(p = list(c(0.40, 0.55), c(0.30, 0.40), c(0.50, 0.60)),
               n_per_arm = c(150, 300))
r <- check_against_analytic(d, g, n_sim = 4000, seed = 101, verbose = FALSE)
print(as.data.frame(r[, c("n_per_arm", "power", "mc_se", "power_analytic",
                          "abs_diff", "within_mc")]), row.names = FALSE)
ok("all scenarios within Monte Carlo error of the exact power", all(r$within_mc),
   sprintf("max |diff| = %.4f", max(r$abs_diff)))

cat("\n== 1b. exact vs power.prop.test: quantifying the approximation gap ==\n")
cat("   (power.prop.test is mildly conservative; this is why the library\n")
cat("    validates against the exact value, not against the approximation)\n")
cmp <- do.call(rbind, lapply(seq_len(nrow(g)), function(i) {
  pp <- resolve_params(d, g[i, , drop = FALSE])
  n <- pp$n_per_arm
  data.frame(n_per_arm = n,
             p_ctrl = pp$p[1], p_trt = pp$p[2],
             exact  = power_two_props_exact(n, n, pp$p[2], pp$p[1], 0.05),
             normal = power_two_props_normal(n, n, pp$p[2], pp$p[1], 0.05))
}))
cmp$gap <- cmp$exact - cmp$normal
print(cmp, row.names = FALSE, digits = 4)

cat("\n== 2. Bonferroni actually lowers the testing level ==\n")
r2 <- run_power(d, data.frame(adjust = c("none", "bonferroni"), n_comparisons = 3),
                n_sim = 3000, seed = 102, verbose = FALSE)
print(as.data.frame(r2[, c("adjust", "alpha_used", "power")]), row.names = FALSE)
ok("alpha divided by n_comparisons", isTRUE(all.equal(r2$alpha_used[2], 0.05 / 3)))
ok("adjusted power is lower", r2$power[2] < r2$power[1])

cat("\n== 3. type I error matches the test's own exact size ==\n")
cat("   (the chi-square test on discrete data is not exactly 0.05; the\n")
cat("    right target is its exact size, which we can enumerate)\n")
r3 <- run_power(d, data.frame(p = I(list(c(0.40, 0.40))), n_per_arm = 300),
                n_sim = 8000, seed = 103, verbose = FALSE)
size_exact <- power_two_props_exact(300, 300, 0.40, 0.40, 0.05)
cat(sprintf("  simulated = %.4f (MC 95%% CI %.4f-%.4f);  exact size = %.4f\n",
            r3$power, r3$power_lo, r3$power_hi, size_exact))
ok("exact size inside the MC interval",
   size_exact >= r3$power_lo && size_exact <= r3$power_hi)
ok("exact size is close to nominal 0.05", abs(size_exact - 0.05) < 0.01,
   sprintf("exact size = %.4f", size_exact))

cat("\n== 4. clustered binary: ICC recovered from the simulated data ==\n")
set.seed(7)
dc <- design_binary_cluster()
dat <- dc$dgp(design_params(dc, icc = 0.10, n_clusters = 400, m = 30))
a1 <- dat[dat$arm == 1, ]                      # ANOVA estimator within one arm
ms <- summary(aov(y ~ factor(cluster), a1))[[1]][, "Mean Sq"]
icc_hat <- (ms[1] - ms[2]) / (ms[1] + (30 - 1) * ms[2])
cat(sprintf("  target ICC = 0.100, recovered = %.3f\n", icc_hat))
ok("recovered ICC within 0.02 of target", abs(icc_hat - 0.10) < 0.02)

cat("\n== 4b. ICC = 0 collapses to independent Bernoulli ==\n")
set.seed(8)
dat0 <- dc$dgp(design_params(dc, icc = 0, n_clusters = 400, m = 30))
a0 <- dat0[dat0$arm == 1, ]
ms0 <- summary(aov(y ~ factor(cluster), a0))[[1]][, "Mean Sq"]
icc0 <- (ms0[1] - ms0[2]) / (ms0[1] + 29 * ms0[2])
cat(sprintf("  recovered ICC = %.4f\n", icc0))
ok("no clustering induced when icc = 0", abs(icc0) < 0.01)

cat("\n== 5. clustered binary power vs the design-effect approximation ==\n")
cat("   (the design effect is itself an approximation, so this is a sanity\n")
cat("    check on the order of magnitude, not an exactness claim)\n")
rc <- check_against_analytic(dc, data.frame(n_clusters = c(20, 40), m = 30, icc = 0.05),
                             n_sim = 1500, seed = 104, verbose = FALSE)
print(as.data.frame(rc[, c("n_clusters", "power", "mc_se", "power_analytic", "abs_diff")]),
      row.names = FALSE)
ok("within 0.06 of the design-effect approximation", all(rc$abs_diff < 0.06),
   sprintf("max |diff| = %.3f", max(rc$abs_diff)))

cat("\n== 5b. unequal cluster sizes: mean and CV are as requested ==\n")
set.seed(9)
dv <- dc$dgp(design_params(dc, n_clusters = 3000, m = 30, m_varies = TRUE, m_cv = 0.5))
sz <- as.numeric(table(dv$cluster))
cat(sprintf("  target mean 30 / CV 0.50  ->  simulated mean %.1f / CV %.2f\n",
            mean(sz), sd(sz) / mean(sz)))
ok("mean cluster size within 5% of target", abs(mean(sz) - 30) / 30 < 0.05)
ok("CV within 0.05 of target", abs(sd(sz) / mean(sz) - 0.5) < 0.05)
ok("unequal sizes reduce power vs equal sizes at the same mean", {
  re <- run_power(dc, data.frame(m_varies = c(FALSE, TRUE)), n_sim = 1500,
                  seed = 106, verbose = FALSE)
  cat(sprintf("  power: equal sizes %.3f, unequal sizes %.3f\n", re$power[1], re$power[2]))
  re$power[2] <= re$power[1]
})

cat("\n== 6. Weibull DGP targets the intended hazard ratio ==\n")
cat("   (the one place a parameterisation slip would silently change the\n")
cat("    effect size being powered for)\n")
set.seed(11)
dw <- design_weibull_bayes()
big <- dw$dgp(design_params(dw, n = 200000, horizon = Inf, competing_dist = "none"))
cf <- unname(coef(survival::coxph(survival::Surv(time, status) ~ trt, data = big)))
cat(sprintf("  true log HR = %.4f, Cox estimate on n = %d: %.4f\n",
            log(1.6), nrow(big), cf))
ok("Cox recovers the true log HR", abs(cf - log(1.6)) < 0.02,
   sprintf("|diff| = %.4f", abs(cf - log(1.6))))

cat("\n== 6b. Weibull shape and control-arm median are as specified ==\n")
ctrl <- big$time[big$trt == 0]
med_target <- unname(weibull_control_summary(4, 0.9)["median"])
cat(sprintf("  target control median = %.3f, simulated = %.3f\n",
            med_target, median(ctrl)))
ok("control-arm median matches", abs(median(ctrl) - med_target) / med_target < 0.03)

cat("\n== 6c. censoring reduces events without changing the estimand ==\n")
set.seed(12)
cen <- dw$dgp(design_params(dw, n = 200000))
cf2 <- unname(coef(survival::coxph(survival::Surv(time, status) ~ trt, data = cen)))
cat(sprintf("  event fraction = %.3f, Cox log HR = %.4f\n", mean(cen$status), cf2))
ok("events lost to competing risk / horizon", mean(cen$status) < 0.95)
ok("log HR still recovered under censoring", abs(cf2 - log(1.6)) < 0.03,
   sprintf("|diff| = %.4f", abs(cf2 - log(1.6))))

cat("\n== 7. reproducibility: same seed, same answer, any worker count ==\n")
ra <- run_power(d, NULL, n_sim = 500, seed = 999, workers = 1, verbose = FALSE)
rb <- run_power(d, NULL, n_sim = 500, seed = 999, workers = 4, verbose = FALSE)
cat(sprintf("  1 worker: %.4f   4 workers: %.4f\n", ra$power, rb$power))
ok("identical across worker counts", isTRUE(all.equal(ra$power, rb$power)))

cat("\n== 8. failed replicates are counted, not silently dropped ==\n")
dbad <- new_design("always_errors", dgp = function(p) data.frame(x = 1),
                   analyze = function(d, p) stop("boom"))
rbad <- run_power(dbad, NULL, n_sim = 20, workers = 1, verbose = FALSE)
ok("errored replicates reported as failed", rbad$n_failed == 20)

# The subtler case: analyze() returns without erroring but `reject` is NA,
# as a non-converged model would. Averaging with na.rm would hide these.
dna <- new_design("half_na", dgp = function(p) data.frame(i = 1),
                  analyze = function(d, p) {
                    if (runif(1) < 0.5) c(reject = NA_real_, x = 1)
                    else c(reject = 1, x = 1)
                  })
rna <- run_power(dna, NULL, n_sim = 400, seed = 55, workers = 1, verbose = FALSE)
cat(sprintf("  n_failed = %d of 400, power = %.3f\n", rna$n_failed, rna$power))
ok("NA rejections counted as failures, not dropped", rna$n_failed > 150)
ok("power computed only over usable replicates", isTRUE(all.equal(rna$power, 1)))

cat("\n== 8b. solve_n() finds the n where exact power crosses the target ==\n")
n_sim_solve <- 3000
n_hat <- solve_n(d, params = list(p = c(0.40, 0.55)), target = 0.80,
                 n_range = c(50, 600), n_sim = n_sim_solve, tol = 4,
                 seed = 33, verbose = FALSE)
# Exact crossing point, by stepping the enumerated power
n_true <- Find(function(n) power_two_props_exact(n, n, 0.55, 0.40, 0.05) >= 0.80,
               50:600)
cat(sprintf("  solve_n = %d (%d sims/step), exact crossing = %d\n",
            n_hat, n_sim_solve, n_true))
ok("solve_n within 15 of the exact crossing", abs(n_hat - n_true) <= 15,
   sprintf("difference = %d", abs(n_hat - n_true)))

cat("\n== 9. Bayesian Weibull: JAGS design runs and behaves sensibly ==\n")
cat("   (small n_sim -- this is a smoke test, the real run is in analyses/)\n")
rj <- run_power(design_weibull_bayes(),
                data.frame(beta = c(0, log(1.6))),
                n_sim = 60, seed = 105, workers = 4, verbose = FALSE)
print(as.data.frame(rj[, c("beta", "power", "mc_se", "prob_benefit",
                           "post_mean_beta", "event_frac", "n_failed")]),
      row.names = FALSE, digits = 3)
ok("no failed JAGS fits", all(rj$n_failed == 0))
ok("posterior mean tracks the truth under the null", abs(rj$post_mean_beta[1]) < 0.15)
ok("power is higher under the alternative than the null",
   rj$power[2] > rj$power[1])

cat(sprintf("\n%d passed, %d failed\n", pass, fail))
if (fail > 0) quit(status = 1)
