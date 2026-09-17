# Power library

A record of the power and sample size analyses the group has run, so that when a
new project starts you can find the closest thing we've already done and work
from it rather than from scratch.

## Index

| Study | Outcome | Framework | Arms | Design | Key design features | Code |
|---|---|---|---|---|---|---|
| **PCORI LOI** (2026)<br><sub>Three-arm pragmatic trial </sub><br><sub>Zhe Chen | Binary | frequentist | 3 | Parallel-arm, individually randomised | Two active arms vs usual care; pairwise comparisons tested at Bonferroni-adjusted significance level; rural/urban subgroups | [analysis.Rmd](analyses/2026-09_pcori-3arm-binary/analysis.Rmd) |
| **PARMA** (2025)<br><sub>Pediatric ARDS, high vs low driving pressure; Bayesian time-to-event model for time to hypoxemia resolution</sub><br><sub>Zhe Chen | Time-to-event | Bayesian | 2 | Parallel-arm, Weibull PH fitted in JAGS | Bayesian survival model with Weibull distribution; Death as a competing event; 28-day administrative censoring; posterior probability decision rule | [analysis.Rmd](analyses/2025-09_parma-weibull-bayes/analysis.Rmd) |

### Column conventions

| Column | What goes in it |
|---|---|
| **Study** | Trial or grant name, year, a one-line description, and analyst on a second line. |
| **Outcome** | Binary / Continuous / Time-to-event / Ordinal / Count |
| **Framework** | frequentist / Bayesian |
| **Arms** | Number of arms |
| **Design** | Parallel-arm / Cluster-randomised / Stepped-wedge / Crossover, and how randomisation works |
| **Key design features** | The things that determine whether this analysis is worth copying, e.g., clustering and ICC, repeated measures, multiplicity, competing risks, sample size constraints |
| **Code** | Link to the analysis document |

## How the documents are built

**Each analysis is a single self-contained R Markdown file.** Nothing is
`source()`d — the simulation loop, the data-generating model, the analysis, and
(for the Bayesian one) the JAGS model all live in the file you're reading. 
Every function is commented in detail.

They share a common structure:

| Section | What it does |
|---|---|
| The question | Trial overview and design features |
| Assumptions | Every parameter the simulation depends on |
| Machinery | The commented helper functions and model specification |
| Results | Power tables and plots |

## Adding an analysis

1. Find the nearest row in the index and copy that folder. Use `YYYY-MM_slug` so
   it sorts chronologically.

   ```
   cp -r analyses/2026-09_pcori-3arm-binary analyses/2026-11_my-trial
   ```

2. Work through it top to bottom, replacing the numbers. The Assumptions chunk
   is where most of the editing happens; `simulate_one()` (or
   `simulate_trial()`) is the only function you are likely to rewrite.

3. Add a row to the index table above.


