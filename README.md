# Power library

A record of the power and sample size analyses the group has run, so that when a
new project starts you can find the closest thing we've already done and work
from it rather than from scratch.

Two parts: an **index** of completed analyses, and a **template** for each design
type to start a new one from.

## Index

| Study | Outcome | Framework | Arms | Design | Key design features | Code |
|---|---|---|---|---|---|---|
| **PCORI LOI** (2026)<br><sub>Three-arm pragmatic trial; sample size, power at budget, rural/urban subgroups</sub><br><sub>Zhe Chen · final</sub> | Binary | frequentist | 3 | Parallel-arm, individually randomised | Two non-symmetric active arms; Bonferroni across 2 primary and 3 pairwise comparisons; prespecified subgroup power at a 1/3 vs 2/3 split; fixed N = 3000 | [analysis.Rmd](analyses/2026-09_pcori-3arm-binary/analysis.Rmd) |
| **PARMA** (2025)<br><sub>Pediatric ARDS, high vs low driving pressure; Bayesian power for time to hypoxemia resolution</sub><br><sub>Zhe Chen, Nadir Yehya · used in submission</sub> | Time-to-event | Bayesian | 2 | Parallel-arm, Weibull PH fitted in JAGS | Death as a competing event, handled as cause-specific censoring; 28-day administrative censoring; posterior-probability decision rule; prior sensitivity | [analysis.Rmd](analyses/2025-09_parma-weibull-bayes/analysis.Rmd) |

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

## Templates

| Template | Use it for |
|---|---|
| [`templates/binary-parallel.Rmd`](templates/binary-parallel.Rmd) | Parallel-arm trial, binary outcome. Multiplicity, subgroup power, exact power validation. |
| [`templates/bayesian-survival.Rmd`](templates/bayesian-survival.Rmd) | Parallel-arm trial, time-to-event outcome, Bayesian analysis. Weibull PH in JAGS, competing risks, posterior-probability decision rule. |

Both run as-is with placeholder numbers, so knit one before changing anything and
confirm you get output.

**Each document is self-contained.** Nothing is `source()`d — the simulation
loop, the analysis, and (for the Bayesian one) the JAGS model all live in the
file you're reading. You can copy a single `.Rmd` somewhere else and it still
works, and you can read one top to bottom without following anything into a
library.

## Adding an analysis

1. Copy the closest template — or the closest completed analysis, if one of them
   is nearer to your design.

   ```
   mkdir analyses/2026-11_my-trial
   cp templates/binary-parallel.Rmd analyses/2026-11_my-trial/analysis.Rmd
   ```

   Use `YYYY-MM_slug` so the folder sorts chronologically.

2. Work through it top to bottom, replacing the placeholder numbers.

3. Add a row to the index table above.

Two conventions worth keeping, both of which the templates already follow:

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

# time-to-event templates only
install.packages(c("survival", "rjags"))   # rjags needs JAGS: brew install jags
```

Rendering to HTML also needs `pandoc` (`brew install pandoc`). Without it you can
still run a document end to end with `knitr::knit("analysis.Rmd")`.

## Gaps

Designs we've hit before and haven't written up yet, roughly in order of how
often they come up: **stepped-wedge**, **continuous outcomes with repeated
measures**, **cluster-randomised with an ICC**, **ordinal outcomes**, and
**group-sequential designs with interim analyses**. Each is a new template.
