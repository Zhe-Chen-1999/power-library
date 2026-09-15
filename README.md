# Power library

A record of the power and sample size analyses the group has run, so that when a
new project starts you can find the closest thing we've already done and work
from it rather than from scratch.

## Index

| Study | Outcome | Framework | Arms | Design | Key design features | Code |
|---|---|---|---|---|---|---|
| **PCORI LOI** (2026)<br><sub>Three-arm pragmatic trial; sample size, power at budget, rural/urban subgroups</sub><br><sub>Zhe Chen | Binary | frequentist | 3 | Parallel-arm, individually randomised | Two non-symmetric active arms; Bonferroni across 2 primary and 3 pairwise comparisons; prespecified subgroup power at a 1/3 vs 2/3 split; fixed N = 3000 | [analysis.Rmd](analyses/2026-09_pcori-3arm-binary/analysis.Rmd) |
| **PARMA** (2025)<br><sub>Pediatric ARDS, high vs low driving pressure; Bayesian power for time to hypoxemia resolution</sub><br><sub>Zhe Chen | Time-to-event | Bayesian | 2 | Parallel-arm, Weibull PH fitted in JAGS | Death as a competing event, handled as cause-specific censoring; 28-day administrative censoring; posterior-probability decision rule; prior sensitivity | [analysis.Rmd](analyses/2025-09_parma-weibull-bayes/analysis.Rmd) |

### Column conventions

Keep these consistent — the value of the index is being able to scan it for
"anything with an ICC" or "any Bayesian time-to-event", and that only works if
entries are described the same way.

| Column | What goes in it |
|---|---|
| **Study** | Trial or grant name, year, a one-line description, then analyst and status on a second line. Status is `draft`, `final`, or `used in submission`. |
| **Outcome** | Binary / Continuous / Time-to-event / Ordinal / Count |
| **Framework** | frequentist / Bayesian |
| **Arms** | Number of arms |
| **Design** | Parallel-arm / Cluster-randomised / Stepped-wedge / Crossover, and how randomisation works |
| **Key design features** | The things that determine whether this analysis is worth copying: clustering and ICC, repeated measures, multiplicity, competing risks, interim analyses, sample size constraints |
| **Code** | Link to the analysis document |

## How the documents are built

**Each analysis is a single self-contained R Markdown file.** Nothing is
`source()`d — the simulation loop, the data-generating model, the analysis, and
(for the Bayesian one) the JAGS model all live in the file you're reading. There
is no package to install and no framework to learn. You can read one top to
bottom, and you can copy one somewhere else and it still runs.

Every function is commented in detail, including why it is written the way it
is, so the parts that are easy to get wrong — the parallel random number seeds,
the Weibull parameterisation, the handling of censored observations — say so
where you meet them.

They share a common shape:

| Section | What it does |
|---|---|
| The question | The design in plain sentences, including where each assumed number came from |
| Assumptions | Every parameter the simulation depends on, in one chunk |
| Machinery | The commented helper functions |
| Validation | Checks the simulation against something computed a different way |
| Results | Power tables and plots |
| Sensitivity | Varies whatever was least certain in "The question" |
| Assumptions and limitations | What the simulation does *not* represent |

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

If nothing in the index is close — a stepped-wedge design, say — start from
whichever is nearest structurally and replace the data-generating function. The
surrounding machinery does not care what the outcome type is.

Two conventions worth keeping, both of which the existing analyses follow:

- **Validate against something.** If the design has a closed form, check the
  simulation against it and show the comparison. If it doesn't, verify the
  data-generating model recovers the effect you asked for, and report the type I
  error under the null. A simulation that produces a number is not the same as
  one that produces the right number.
- **Keep `kable()` out of cached chunks.** Put the computation in a chunk with
  `cache = TRUE` and the table in a separate chunk. Otherwise fixing a typo
  re-runs the simulation, which for a JAGS design means waiting fifteen minutes.

## Runtime

The binary designs are effectively free — thousands of replicates a second, and
the PCORI document knits in about 15 seconds.

JAGS is not. One fit is roughly 60 ms, so a 1000-replicate scenario takes about a
minute on 8 cores, and the PARMA document is around 13,000 fits — about 15
minutes from a cold cache. While exploring, drop `n_sim` and `n_iter`; raise them
only for the numbers that go into the application. The reported `mc_se` tells you
when they're high enough.

## Requirements

R ≥ 4.4, and:

```r
install.packages(c("future.apply", "knitr", "ggplot2", "rmarkdown"))

# time-to-event analyses only
install.packages(c("survival", "rjags"))   # rjags needs JAGS: brew install jags
```

Rendering to HTML also needs `pandoc` (`brew install pandoc`). Without it you can
still run a document end to end with `knitr::knit("analysis.Rmd")`.

## Gaps

Designs we've hit before and haven't written up yet, roughly in order of how
often they come up: **stepped-wedge**, **continuous outcomes with repeated
measures**, **cluster-randomised with an ICC**, **ordinal outcomes**, and
**group-sequential designs with interim analyses**. Each is a new analysis
folder, started from whichever existing one is closest.
