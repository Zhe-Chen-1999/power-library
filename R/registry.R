## =====================================================================
## registry.R -- the index Yingying asked for
##
## registry/analyses.yml is the single source of truth: one entry per
## analysis, describing the study, the outcome, the design features that
## determine which past analysis is worth copying, and where the code is.
##
##   build_index()   regenerates the table in README.md in place
##   add_entry()     appends a stub entry for a new analysis
##
## Keep the columns stable. Once there are twenty entries, the value of
## this file is being able to sort it -- "show me everything with an ICC",
## "show me every Bayesian time-to-event" -- and that only works if the
## same fields are filled in every time.
## =====================================================================

suppressPackageStartupMessages(library(yaml))

REGISTRY_FIELDS <- c(
  id          = "short slug, matches the analyses/ subdirectory",
  title       = "what the analysis was for",
  study       = "trial or grant name",
  year        = "when it was run",
  analyst     = "who ran it",
  design      = "e.g. parallel-arm RCT, cluster-randomised, stepped wedge",
  outcome     = "binary / continuous / time-to-event / ordinal / count",
  framework   = "frequentist / Bayesian",
  n_arms      = "number of arms",
  features    = "key design features: ICC, repeated measures, multiplicity, competing risks, interim analyses",
  design_fn   = "design constructor in R/ that implements it, if any",
  path        = "link to the Markdown or code",
  status      = "draft / final / used-in-submission"
)

registry_file <- function() file.path(POWERLIB_ROOT, "registry", "analyses.yml")

read_registry <- function(path = registry_file()) {
  if (!file.exists(path)) return(list())
  yaml::read_yaml(path)$analyses %||% list()
}

#' Registry as a data.frame, one row per analysis
registry_table <- function(path = registry_file()) {
  entries <- read_registry(path)
  if (!length(entries)) return(data.frame())
  cols <- names(REGISTRY_FIELDS)
  do.call(rbind, lapply(entries, function(e) {
    row <- lapply(cols, function(k) as.character(e[[k]] %||% ""))
    names(row) <- cols
    as.data.frame(row, stringsAsFactors = FALSE)
  }))
}

#' Append a stub entry, so adding an analysis is a one-liner rather than
#' a hand-edit of YAML that someone will forget to do.
add_entry <- function(id, title, study, year = format(Sys.Date(), "%Y"),
                      analyst = "", design = "", outcome = "",
                      framework = "frequentist", n_arms = "", features = "",
                      design_fn = "", path = "", status = "draft",
                      file = registry_file()) {
  doc <- if (file.exists(file)) yaml::read_yaml(file) else list(analyses = list())
  new <- list(id = id, title = title, study = study, year = year,
              analyst = analyst, design = design, outcome = outcome,
              framework = framework, n_arms = n_arms, features = features,
              design_fn = design_fn,
              path = if (nzchar(path)) path else paste0("analyses/", id, "/analysis.Rmd"),
              status = status)
  doc$analyses <- c(doc$analyses, list(new))
  yaml::write_yaml(doc, file)
  message("added '", id, "' to the registry; run build_index() to refresh README.md")
  invisible(new)
}

#' Render the registry as a Markdown table
#'
#' Deliberately narrow: study / outcome / framework / arms / key features /
#' link. The full YAML has everything else -- this is the thing you scan.
registry_markdown <- function(path = registry_file()) {
  tab <- registry_table(path)
  if (!nrow(tab)) return("_No analyses registered yet._")

  hdr <- c("Study", "Outcome", "Framework", "Arms", "Design", "Key design features", "Code")
  rows <- apply(tab, 1, function(r) {
    link <- if (nzchar(r[["path"]])) sprintf("[%s](%s)", basename(r[["path"]]), r[["path"]]) else ""
    # Analyst and status ride along under the study name rather than as their
    # own columns: you need both before reusing something (who to ask, and
    # whether it was ever finished), but neither is what you scan the table
    # for, and seven columns is already as wide as this reads well.
    byline <- paste(c(r[["analyst"]], r[["status"]])[nzchar(c(r[["analyst"]], r[["status"]]))],
                    collapse = " &middot; ")
    study <- paste0("**", r[["study"]], "** (", r[["year"]], ")",
                    "<br><sub>", r[["title"]], "</sub>",
                    if (nzchar(byline)) paste0("<br><sub>", byline, "</sub>") else "")
    sprintf("| %s | %s | %s | %s | %s | %s | %s |",
            study, r[["outcome"]], r[["framework"]], r[["n_arms"]], r[["design"]],
            r[["features"]], link)
  })

  paste(c(
    paste0("| ", paste(hdr, collapse = " | "), " |"),
    paste0("|", paste(rep("---", length(hdr)), collapse = "|"), "|"),
    rows
  ), collapse = "\n")
}

#' Regenerate the index table inside README.md
#'
#' Rewrites only what is between the marker comments, so prose around it
#' survives.
build_index <- function(readme = file.path(POWERLIB_ROOT, "README.md"),
                        path = registry_file()) {
  start <- "<!-- INDEX:START -->"
  end   <- "<!-- INDEX:END -->"
  lines <- readLines(readme, warn = FALSE)
  i <- which(lines == start); j <- which(lines == end)
  if (length(i) != 1 || length(j) != 1 || j <= i) {
    stop("README.md must contain exactly one ", start, " ... ", end, " block")
  }
  new <- c(lines[seq_len(i)], "", registry_markdown(path), "", lines[j:length(lines)])
  writeLines(new, readme)
  message("README.md index rebuilt (", nrow(registry_table(path)), " analyses)")
  invisible(TRUE)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
