# Power library

Shared, reusable power and sample size analyses for the group.

Two problems this is trying to solve. First, we keep re-deriving the same
calculations from scratch — the PCORI three-arm binary tables and the PARMA
Bayesian survival simulation had almost nothing in common as code, but about
80% in common as *structure*. Second, when a new project starts, nobody can
easily find the closest thing we've already done. The index below is the answer
to the second problem; the engine in `R/` is the answer to the first.

## The index

Scan this to find the nearest precedent for a new project. Generated from
[`registry/analyses.yml`](registry/analyses.yml) — edit that file, then run
`source("R/load.R"); build_index()`.

<!-- INDEX:START -->

| Study | Outcome | Framework | Arms | Design | Key design features | Code |
|---|---|---|---|---|---|---|
| **PCORI LOI** (2026)<br><sub>Sample size and power for a three-arm pragmatic trial, plus rural/urban subgroup power</sub> | Binary | frequentist | 3 | Parallel-arm RCT, individually randomised | Two non-symmetric active arms; Bonferroni across 2 primary and 3 pairwise comparisons; prespecified subgroup power at a 1/3 vs 2/3 split; fixed N = 3000 budget | [analysis.Rmd](analyses/2026-09_pcori-3arm-binary/analysis.Rmd) |
| **PARMA (pediatric ARDS, high vs low driving pressure)** (2025)<br><sub>Bayesian power to detect a >90% posterior probability of benefit in time to hypoxemia resolution</sub> | Time-to-event | Bayesian | 2 | Parallel-arm RCT, Bayesian analysis | Weibull PH model in JAGS; death as competing event handled as cause-specific censoring; 28-day administrative censoring; posterior-probability decision rule; prior sensitivity on the treatment-effect precision | [analysis.Rmd](analyses/2025-09_parma-weibull-bayes/analysis.Rmd) |

<!-- INDEX:END -->

## Quick start

```r
source("R/load.R")

d <- design_binary_parallel()

# Power across a grid of scenarios
run_power(d,
          scenarios(p = list(c(0.40, 0.50), c(0.40, 0.55)),
                    n_per_arm = c(200, 400)),
          n_sim = 2000)

# Smallest n per arm reaching 80% power
solve_n(d, params = list(p = c(0.40, 0.55)), target = 0.80)

# Validate the simulation against the exact power of the test it runs
check_against_analytic(d, data.frame(n_per_arm = c(150, 300)))
```

## How it works

The engine knows nothing about any particular design. A **design** is three
things:

```r
new_design(
  name     = "my_design",
  dgp      = function(params) { ... },        # -> a data.frame: one trial
  analyze  = function(data, params) { ... },  # -> named vector incl. `reject`
  defaults = list(...),                       # parameters, overridable per scenario
  analytic = function(params) { ... },        # optional closed form, for validation
  meta     = list(design = , outcome = , features = )
)
```

`run_power()` replicates that pair over a grid of scenarios and averages
`reject`. **What `reject` means is the design's business.** For a frequentist
test it is `p < alpha`. For the Bayesian survival design it is
`Pr(HR > 1 | data) >= 0.90`. The engine doesn't care, which is why one piece of
machinery covers both, and why adding a stepped-wedge or ordinal design means
writing two functions rather than another standalone script.

Everything else `analyze()` returns is averaged and carried along, so posterior
means, event counts, realised effect sizes and convergence diagnostics come out
of the same call as the power.

### What you get back

`run_power()` returns the scenario grid plus:

| Column | Meaning |
|---|---|
| `power` | proportion of replicates with `reject == 1` |
| `mc_se`, `power_lo`, `power_hi` | Monte Carlo error of that estimate — the *simulation's* precision, not the trial's |
| `n_sim`, `n_failed` | replicates attempted, and how many errored |
| *(others)* | mean of every other quantity `analyze()` returned |

`n_failed` matters. Failed replicates are counted and reported, never silently
dropped — a design that only converges on the easy scenarios would otherwise
look better than it is.

## Validation

Every design with a closed-form counterpart is checked against it:

```r
check_against_analytic(design_binary_parallel(),
                       data.frame(n_per_arm = c(150, 300)))
```

For the binary designs the comparison is against the **exact** power, computed
by enumerating every 2×2 table the trial could produce
(`power_two_props_exact()`) — not against `power.prop.test()`, which is a normal
approximation and runs about 1–2 points conservative at moderate n. The exact
value is the true power of the test the simulation actually performs, so
agreement is a real correctness check rather than two approximations agreeing
with each other.

`power_two_props_normal()` reproduces `power.prop.test()` when you want the
conventional number for a protocol.

Run the full suite:

```
Rscript tests/test_engine.R
```

## Layout

```
R/
  engine.R              run_power(), scenarios(), solve_n(), check_against_analytic()
  dgp_binary.R          binary outcomes: parallel-arm, cluster-randomised
  dgp_weibull_bayes.R   Bayesian Weibull survival (JAGS) + Cox comparator
  registry.R            the index table
  load.R                source this
jags/
  weibull_ph.jags       Weibull PH model; N and priors passed as data
analyses/
  <year-month>_<slug>/  one directory per analysis, with analysis.Rmd
registry/
  analyses.yml          the index, one entry per analysis
tests/
  test_engine.R
```

## Adding an analysis

1. `mkdir analyses/2026-11_my-trial` and write `analysis.Rmd` there. Start from
   whichever existing analysis is closest — that's what the index is for.
2. Source the library with `source("../../R/load.R")`.
3. If the design is new, add a `design_*()` constructor to `R/` rather than
   defining it inline in the Rmd. That is the difference between a library and
   a folder of scripts.
4. Register it and refresh the index:

```r
source("R/load.R")
add_entry(
  id        = "2026-11_my-trial",
  title     = "Power for the primary endpoint",
  study     = "MY-TRIAL",
  analyst   = "...",
  design    = "Stepped-wedge cluster-randomised",
  outcome   = "Binary",
  framework = "frequentist",
  n_arms    = "2",
  features  = "12 clusters, 6 steps, ICC 0.03, CAC 0.8"
)
build_index()
```

Keep the `features` field specific and keep filling in every column. Once there
are twenty entries, the value of the index is being able to answer "show me
everything we've done with an ICC" or "show me every Bayesian time-to-event
analysis" — and that only works if the fields are populated consistently.

## Designs currently available

| Constructor | Design | Outcome | Notes |
|---|---|---|---|
| `design_binary_parallel()` | Parallel-arm RCT | Binary | 2+ arms, arm-specific event rates, Bonferroni, exact power available |
| `design_binary_cluster()` | Cluster-randomised | Binary | ICC exact via beta-binomial, unequal cluster sizes, cluster-level t-test or GLMM |
| `design_weibull_bayes()` | Parallel-arm, Bayesian | Time-to-event | Weibull PH in JAGS, competing risks, posterior-probability decision rule |
| `design_weibull_cox()` | Parallel-arm, frequentist | Time-to-event | Same DGP as above; isolates the effect of the analysis model |

Obvious gaps, roughly in order of how often we hit them: **stepped-wedge**
(Yingying has run these), **continuous outcomes with repeated measures**,
**ordinal outcomes** (proportional odds), and **group-sequential / interim
analyses**. Each is a `design_*()` constructor away.

## Requirements

R ≥ 4.4, and:

```r
install.packages(c("future.apply", "tibble", "dplyr", "yaml",
                   "lme4", "survival", "ggplot2", "knitr", "rmarkdown"))
install.packages("rjags")   # requires JAGS: brew install jags
```

`rjags` is only needed for the Bayesian survival designs; everything else loads
without it.
