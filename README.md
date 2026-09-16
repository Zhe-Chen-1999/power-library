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
(for the Bayesian one) the JAGS model all live in the file you're reading. 
Every function is commented in detail.

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

Conventions worth keeping:

- **Validate against something.** If the design has a closed form, check the
  simulation against it and show the comparison. If it doesn't, verify the
  data-generating model recovers the effect you asked for, and report the type I
  error under the null. A simulation that produces a number is not the same as
  one that produces the right number.

