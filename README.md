# Group Sequential Design Explorer

A Shiny web application for designing group sequential clinical trials,
built on the open-source R package **gsDesign** (GPL >= 3). It replicates
the feature set of the gsDesign Explorer web app: create a group sequential
design for time-to-event, binomial, or normal endpoints (or extend a
fixed-design sample size), inspect boundary plots, a tabular summary, a
textual summary, and a re-runnable R Markdown report, update bounds at
the time of analysis, and save / restore designs.

The app is an original implementation: it contains no code copied from any
other web application. All statistical calculations call the gsDesign
package through its public API.

## Credits

gsDesign was created by **Keaven Anderson** and is maintained by
**Nan Xiao**. See https://cran.r-project.org/package=gsDesign and
https://keaven.github.io/gsDesign/.

## Features

- **Create Design**
  - Design parameters in tabs: Sample Size, Timing, Boundaries,
    Hazard Ratios, Failure Rates, Enrollment.
  - Endpoints: time-to-event (`gsSurv` / `gsSurvCalendar`), binomial
    (`nBinomial` + `gsDesign`), normal (`nNormal` + `gsDesign`), and
    fixed-design N.
  - Time-to-event: variable to adjust to obtain desired power — sample
    size (gsSurv scales the enrollment rate), number of events
    (`nSurvival` + `gsDesign`), or study duration (gsSurv solves for T);
    calendar-based timing (`gsSurvCalendar`); enrollment ramp-up
    (piecewise rates); control failure rate specified by median, failure
    rate, or survival rate; censoring (dropout) rate.
  - Timing: information-based (equal or custom information fractions) or
    calendar-based; any number of interim analyses.
  - Boundaries: one-sided, two-sided symmetric, or two-sided asymmetric
    (binding / non-binding futility); upper and lower spending functions —
    O'Brien-Fleming and Pocock (Lan-DeMets), Hwang-Shih-DeCani, power
    family, logistic, and custom spending points; total H0 spending
    summary.
  - Design output: boundary plots (Z-scale, nominal p-value, cumulative
    power, enrollment/events), tabular summary (`gsBoundSummary`,
    downloadable as CSV), textual summary, and an R Markdown report
    generator producing re-runnable gsDesign code (downloadable as
    .Rmd or rendered .html).
- **Update Bounds**: update bounds at the time of analysis by re-running
  `gsDesign()` with the observed information (`n.I`) and the planned
  maximum information (`maxn.IPlan`); flexible total number of analyses;
  final-analysis alpha-spending strategy (spending function vs. spend all
  remaining alpha via `usTime`); tabular and R Markdown report output.
- **Save / Restore**: download a design as .rds (parameters + computed
  design object) or .json (parameters only); restore by uploading an .rds
  file.
- **About**: feature summary, contributors, license, references.

## Run locally

```r
# Install dependencies (once)
install.packages(c("shiny", "gsDesign", "ggplot2", "DT", "jsonlite",
                   "rmarkdown"))

# From R, in this directory:
shiny::runApp()
```

Requires R >= 4.0 and gsDesign >= 3.0 (developed and validated against
R 4.3.3 and gsDesign 3.11.0).

## Validation

`test_validation.R` checks the gsDesign API usage and the app's own
helpers, including that every design type attains its target power and
that the generated R code reproduces the design exactly:

```sh
Rscript test_validation.R
```

Expected result: 68 checks, all PASS (see `validation_output.txt`).
Key numbers checked:

- 2-look O'Brien-Fleming (Lan-DeMets), alpha 0.025 one-sided, power 0.9:
  critical values **2.9626 / 1.9686**.
- gsDesign documented gsSurv example: N = 460, events = 97.
- Fixed-design N: `nNormal(delta1 = 1, sd = 1)` = 42.03 total;
  `nBinomial(p1 = 0.3, p2 = 0.45)` = 433.64 total.
- All endpoint / spending-function / adjust-mode combinations built
  through `build_design()` attain power 0.9 (within 0.01).
- Bound update: re-running `gsDesign(n.I, maxn.IPlan)` keeps the first
  bound and shifts later ones; the spend-all final-analysis strategy
  (final bound 1.977) differs from the spending-function strategy
  (2.123).

## Deploy to shinyapps.io

1. Create a free account at https://www.shinyapps.io/ and install the
   `rsconnect` package:

   ```r
   install.packages("rsconnect")
   ```

2. Get your account token and secret from the shinyapps.io dashboard
   (Account → Tokens). **Do not share them in chat or commit them to
   the repository.** Set them from R (once per machine):

   ```r
   rsconnect::setAccountInfo(name = "<account>",
                             token = "<token>",
                             secret = "<secret>")
   ```

3. Deploy this directory:

   ```r
   rsconnect::deployApp(appDir = ".", appName = "gsdesign-explorer")
   ```

   The deployed app will be at
   `https://<account>.shinyapps.io/gsdesign-explorer/`.

   Alternatively, use GitHub Actions to deploy on push (see
   `.github/workflows/shinyapps-deploy.yml` if present in the repo).
   The workflow needs the repository secrets `SHINYAPPS_TOKEN` and
   `SHINYAPPS_SECRET`, set by the account owner in the GitHub
   repository settings — never in code.

## Repository

Source lives at https://github.com/hutianle/githubtest (the app is
`app.R` in this directory). Push and deployment are done by the
repository owner; this directory is the working copy only.
