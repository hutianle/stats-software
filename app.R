# ============================================================================
# Group Sequential Design Explorer
# -----------------------------------------------------------------------------
# An original Shiny web application for designing group sequential clinical
# trials. All statistical calculations are performed by the open-source R
# package gsDesign (GPL >= 3) through its public API; this app contains no
# code copied from any other web application.
#
# Endpoints supported: time-to-event (gsSurv / gsSurvCalendar), binomial
# (nBinomial + gsDesign), normal (nNormal + gsDesign), and fixed-design N.
#
# Credits: gsDesign was created by Keaven Anderson and is maintained by
# Nan Xiao. See the About tab for details.
# ============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(gsDesign)
  library(ggplot2)
  library(DT)
  library(jsonlite)
})

# ---------------------------------------------------------------------------
# Helper: spending functions
# ---------------------------------------------------------------------------
# Each spending function is identified by a short key. `spending_fn()`
# returns a list with the function itself (`fn`) and its parameter (`par`)
# as expected by gsDesign()/gsSurv() through the sfu/sfupar and sfl/sflpar
# arguments.
#
# Notes on gsDesign 3.x:
# - There is no sfOF or sfPocock in gsDesign 3.11.0. The classic
#   O'Brien-Fleming and Pocock boundaries are the Lan-DeMets spending
#   functions sfLDOF and sfLDPocock.
# - sfPoints takes a vector of cumulative spending proportions, one per
#   interim analysis (length k - 1); the final analysis spends the rest.

SPENDING_CHOICES <- c(
  "O'Brien-Fleming (Lan-DeMets)" = "LDOF",
  "Pocock (Lan-DeMets)"          = "LDPocock",
  "Hwang-Shih-DeCani"            = "HSD",
  "Power family"                 = "Power",
  "Logistic"                     = "Logistic",
  "Custom points"                = "Points"
)

spending_fn <- function(name, param = list()) {
  switch(name,
    "LDOF"     = list(fn = sfLDOF,     par = NULL),
    "LDPocock" = list(fn = sfLDPocock, par = NULL),
    "HSD"      = list(fn = sfHSD,      par = param$gamma),
    "Power"    = list(fn = sfPower,    par = param$rho),
    "Logistic" = list(fn = sfLogistic, par = param$nu),
    "Points"   = list(fn = sfPoints,   par = param$points),
    list(fn = NULL, par = NULL)
  )
}

spending_label <- function(name) {
  nm <- names(SPENDING_CHOICES)[SPENDING_CHOICES == name]
  if (length(nm) == 0) name else nm
}

# Map a spending-function object back to its name in the gsDesign package,
# for generating re-runnable R code. (deparse() on a function object would
# expand the whole function body, not its name.)
sf_fn_name <- function(fn) {
  if (is.null(fn)) return("NULL")
  ns <- asNamespace("gsDesign")
  for (nm in c("sfLDOF", "sfLDPocock", "sfHSD", "sfPower", "sfLogistic",
               "sfPoints", "sfLinear", "sfExponential", "sfTDist",
               "sfTruncated", "sfXG1", "sfBetaDist")) {
    if (exists(nm, envir = ns, inherits = FALSE) &&
        identical(fn, get(nm, envir = ns))) {
      return(nm)
    }
  }
  "NULL"
}

# ---------------------------------------------------------------------------
# Helper: parsing
# ---------------------------------------------------------------------------
# Parse a comma-separated string of numbers, e.g. "0.5, 0.75" -> c(0.5, 0.75).
parse_num_vector <- function(x) {
  if (is.null(x) || length(x) == 0) return(numeric(0))
  parts <- strsplit(as.character(x), ",")[[1]]
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  suppressWarnings(as.numeric(parts))
}

# Lower bounds for gsProbability(): designs with no lower bound (e.g.
# test.type 1) have d$lower$bound == NULL; gsProbability needs a finite
# vector, so use a very low value that is effectively never crossed.
lower_bounds <- function(d) {
  lb <- d$lower$bound
  if (is.null(lb)) rep(-10, d$k) else lb
}

# ---------------------------------------------------------------------------
# Helper: power of a design at its design alternative
# ---------------------------------------------------------------------------
# Cumulative probability of crossing the upper bound at or before the final
# analysis, under theta = design$delta. gsProbability() returns a list whose
# $upper$prob is a k x 1 matrix of per-analysis crossing probabilities, so
# total power is the sum of its entries.
design_power <- function(d) {
  p <- tryCatch(
    gsProbability(k = d$k, theta = d$delta, n.I = d$n.I,
                  a = lower_bounds(d), b = d$upper$bound),
    error = function(e) NULL
  )
  if (is.null(p)) return(NA_real_)
  sum(p$upper$prob)
}

# Per-analysis (not cumulative) probability of crossing the upper bound.
crossing_probs <- function(d, theta) {
  p <- tryCatch(
    gsProbability(k = d$k, theta = theta, n.I = d$n.I,
                  a = lower_bounds(d), b = d$upper$bound),
    error = function(e) NULL
  )
  if (is.null(p)) return(rep(NA_real_, d$k))
  as.numeric(p$upper$prob)
}

# ---------------------------------------------------------------------------
# Fixed-design sample sizes
# ---------------------------------------------------------------------------
# nBinomial() and nNormal() return a numeric scalar (total sample size),
# not a list.
#
# Binomial fixed-design N (total) and the standardised effect (delta) used
# by gsDesign, for the three scales supported by nBinomial() ("Difference",
# "RR", "OR"). delta0 is on the natural scale of the comparison: a rate
# difference for "Difference", log relative risk for "RR", and log odds
# ratio for "OR" (0 = no difference in all cases).
binomial_fixed <- function(p1, p2, alpha, beta, ratio, sided, scale, delta0) {
  nb <- nBinomial(p1 = p1, p2 = p2, alpha = alpha, beta = beta,
                  ratio = ratio, sided = sided, scale = scale, delta0 = delta0)
  n_total <- as.numeric(nb)
  pbar <- (p1 + ratio * p2) / (1 + ratio)
  delta <- switch(scale,
    "Difference" = (p2 - p1 - delta0) * sqrt(ratio) /
                     ((1 + ratio) * sqrt(pbar * (1 - pbar))),
    "RR" = (log(p2 / p1) - delta0) /
             sqrt((1 + ratio) * ((1 - p1) / p1 + (1 - p2) / (ratio * p2))),
    "OR" = (log((p2 / (1 - p2)) / (p1 / (1 - p1))) - delta0) /
             sqrt((1 + ratio) * (1 / (p1 * (1 - p1)) +
                                   1 / (ratio * p2 * (1 - p2))))
  )
  list(n = n_total, delta = delta)
}

# Normal fixed-design N (total) and standardised effect (delta).
normal_fixed <- function(effect, sd_c, sd_e, alpha, beta, ratio, sided, delta0) {
  nn <- nNormal(delta1 = effect, sd = sd_c, sd2 = sd_e, alpha = alpha,
                beta = beta, ratio = ratio, sided = sided)
  n_total <- as.numeric(nn)
  delta <- (effect - delta0) / sqrt((1 + ratio) * (sd_c^2 + sd_e^2 / ratio))
  list(n = n_total, delta = delta)
}

# Time-to-event: control hazard from the chosen specification (exponential).
control_hazard <- function(spec, median, fail_rate, fail_time,
                           surv_rate, surv_time) {
  switch(spec,
    "median" = log(2) / median,
    "failrate" = -log(1 - fail_rate) / fail_time,
    "survrate" = -log(surv_rate) / surv_time
  )
}

# Piecewise enrollment (ramp-up): half the steady-state rate during the
# ramp-up period, then the steady-state rate. gsSurv scales both rates by
# the same factor to obtain the desired power.
enrollment_gamma_R <- function(rate, rampup, duration) {
  if (!is.null(rampup) && rampup > 0) {
    list(gamma = c(rate / 2, rate), R = c(rampup, duration))
  } else {
    list(gamma = rate, R = duration)
  }
}

# ---------------------------------------------------------------------------
# Design builder
# ---------------------------------------------------------------------------
# `p` is a named list of design parameters (mirroring the UI inputs).
# Returns a list with:
#   $design  - the gsDesign object (class gsDesign, or gsSurv for TTE)
#   $info    - extra information for display (fixed-design N, events, ...)
#   $error   - error message if the design could not be built, else NULL
#
# Note: gsDesign() has no `ratio` argument (the randomization ratio only
# enters through n.fix / delta and, for TTE, through gsSurv()/nSurvival()),
# so ratio is never passed to gsDesign().
build_design <- function(p) {
  k <- p$n_interim + 1L
  if (p$equal_spacing) {
    timing <- (1:k) / k
  } else {
    fr <- parse_num_vector(p$info_frac)
    fr <- fr[!is.na(fr)]
    if (length(fr) != k - 1L || any(fr <= 0) || any(fr >= 1) || is.unsorted(fr)) {
      return(list(design = NULL, info = NULL,
                  error = paste0("Information fractions must be ", k - 1L,
                                 " increasing values between 0 and 1.")))
    }
    timing <- c(fr, 1)
  }

  sfu <- spending_fn(p$sfu_name,
                     list(gamma = p$sfu_gamma, rho = p$sfu_rho, nu = p$sfu_nu,
                          points = parse_num_vector(p$sfu_points_a)))
  if (p$sfl_name == "None") {
    sfl <- list(fn = NULL, par = NULL)
  } else {
    sfl <- spending_fn(p$sfl_name,
                       list(gamma = p$sfl_gamma, rho = p$sfl_rho,
                            nu = p$sfl_nu,
                            points = parse_num_vector(p$sfl_points_a)))
  }

  test.type <- switch(p$testing,
    "one-sided" = 1L,
    "two-sided symmetric" = 2L,
    "two-sided asymmetric, non-binding" = 6L,
    "two-sided asymmetric, binding" = 5L,
    1L
  )
  # For a two-sided symmetric design the lower bound mirrors the upper bound.
  if (test.type == 2L) sfl <- sfu

  alpha <- p$alpha
  beta <- 1 - p$power
  sided <- if (test.type == 1L) 1L else 2L

  info <- list(endpoint = p$endpoint, k = k, timing = timing,
               test.type = test.type, alpha = alpha, beta = beta,
               ratio = p$ratio, time_units = p$time_units,
               sfu_name = spending_label(p$sfu_name),
               sfl_name = if (p$sfl_name == "None") "None"
                          else spending_label(p$sfl_name),
               timing_type = p$timing_type, adjust_var = p$adjust_var)

  d <- tryCatch({
    if (p$endpoint == "Time-to-event") {
      lambdaC <- control_hazard(p$control_spec, p$control_median,
                                p$control_fail_rate, p$control_fail_time,
                                p$control_surv_rate, p$control_surv_time)
      info$lambdaC <- lambdaC
      enr <- enrollment_gamma_R(p$enroll_rate, p$rampup, p$enroll_duration)
      if (identical(p$timing_type, "Calendar-based")) {
        # Analyses at fixed calendar times; spending by information
        # fraction at those times (gsSurvCalendar default).
        cal_times <- parse_num_vector(p$calendar_times)
        cal_times <- cal_times[!is.na(cal_times)]
        if (length(cal_times) < 1L || any(cal_times <= 0) ||
            is.unsorted(cal_times, strictly = TRUE)) {
          stop("Calendar times must be increasing positive values.")
        }
        info$k <- length(cal_times)
        info$timing <- NULL
        info$calendar_mode <- TRUE
        d <- gsSurvCalendar(test.type = test.type, alpha = alpha, beta = beta,
                            sfu = sfu$fn, sfupar = sfu$par,
                            sfl = sfl$fn, sflpar = sfl$par,
                            lambdaC = lambdaC, hr = p$hr, hr0 = p$hr0,
                            eta = p$eta, gamma = enr$gamma, R = enr$R,
                            minfup = p$minfup, ratio = p$ratio,
                            calendarTime = cal_times)
        info$n_total <- max(d$N)
        info$n_events <- d$n.I[d$k]
        info$gamma_derived <- d$gamma
        info$calendar_times <- d$T
        d
      } else if (identical(p$adjust_var, "Number of events")) {
        # Fixed-design events and N from nSurvival(), then a group
        # sequential design on the event (information) scale, following the
        # documented gsDesign pattern: delta is left at 0 so gsDesign
        # computes the standardised effect giving the desired power.
        ss <- nSurvival(lambda1 = lambdaC, lambda2 = lambdaC * p$hr,
                        Ts = p$study_duration, Tr = p$enroll_duration,
                        eta = p$eta, ratio = p$ratio, alpha = alpha,
                        beta = beta, sided = sided)
        info$events_mode <- TRUE
        info$n_total <- ss$n
        info$n_events <- ss$nEvents
        d <- gsDesign(k = k, test.type = test.type, alpha = alpha,
                      beta = beta, timing = timing, sfu = sfu$fn,
                      sfupar = sfu$par, sfl = sfl$fn, sflpar = sfl$par,
                      n.fix = ss$nEvents, nFixSurv = ss$n,
                      delta1 = log(p$hr))
        # gsBoundSummary reads x$hr0 for HR-at-bound calculations when
        # logdelta = TRUE; a plain gsDesign object has no hr0 field.
        d$hr0 <- p$hr0
        info$delta <- d$delta
        d
      } else if (identical(p$adjust_var, "Study duration")) {
        # T = NULL: gsSurv solves for the study duration with the
        # enrollment rate held fixed (and extends R as needed).
        d <- gsSurv(k = k, test.type = test.type, alpha = alpha, beta = beta,
                    timing = timing, sfu = sfu$fn, sfupar = sfu$par,
                    sfl = sfl$fn, sflpar = sfl$par,
                    lambdaC = lambdaC, hr = p$hr, hr0 = p$hr0, eta = p$eta,
                    gamma = enr$gamma, R = enr$R,
                    T = NULL, minfup = p$minfup, ratio = p$ratio)
        info$n_total <- max(d$N)
        info$n_events <- d$n.I[d$k]
        info$gamma_derived <- d$gamma
        info$study_duration_solved <- d$T[d$k]
        info$enroll_duration_solved <- d$R
        d
      } else {
        # "Sample size": gsSurv scales the enrollment rate (gamma) to
        # obtain the desired power with the given enrollment duration.
        d <- gsSurv(k = k, test.type = test.type, alpha = alpha, beta = beta,
                    timing = timing, sfu = sfu$fn, sfupar = sfu$par,
                    sfl = sfl$fn, sflpar = sfl$par,
                    lambdaC = lambdaC, hr = p$hr, hr0 = p$hr0, eta = p$eta,
                    gamma = enr$gamma, R = enr$R,
                    T = p$study_duration, minfup = p$minfup, ratio = p$ratio)
        info$n_total <- max(d$N)
        info$n_events <- d$n.I[d$k]
        info$gamma_derived <- d$gamma
        info$study_duration <- d$T[d$k]
        d
      }
    } else if (p$endpoint == "Binomial") {
      bf <- binomial_fixed(p$p1, p$p2, alpha, beta, p$ratio, sided,
                           p$bin_scale, p$delta0)
      info$n_total <- bf$n
      info$delta <- bf$delta
      gsDesign(k = k, test.type = test.type, alpha = alpha, beta = beta,
               timing = timing, sfu = sfu$fn, sfupar = sfu$par,
               sfl = sfl$fn, sflpar = sfl$par,
               n.fix = bf$n, delta = bf$delta, delta0 = 0, delta1 = bf$delta)
    } else if (p$endpoint == "Normal") {
      nf <- normal_fixed(p$effect, p$sd_c, p$sd_e, alpha, beta,
                         p$ratio, sided, p$delta0)
      info$n_total <- nf$n
      info$delta <- nf$delta
      gsDesign(k = k, test.type = test.type, alpha = alpha, beta = beta,
               timing = timing, sfu = sfu$fn, sfupar = sfu$par,
               sfl = sfl$fn, sflpar = sfl$par,
               n.fix = nf$n, delta = nf$delta, delta0 = 0, delta1 = nf$delta)
    } else { # "Fixed design N"
      info$n_total <- p$n_fixed
      info$delta <- p$delta_std
      gsDesign(k = k, test.type = test.type, alpha = alpha, beta = beta,
               timing = timing, sfu = sfu$fn, sfupar = sfu$par,
               sfl = sfl$fn, sflpar = sfl$par,
               n.fix = p$n_fixed, delta = p$delta_std,
               delta0 = 0, delta1 = p$delta_std)
    }
  }, error = function(e) e)

  if (inherits(d, "error")) {
    return(list(design = NULL, info = info,
                error = paste("Could not build design:", conditionMessage(d))))
  }
  list(design = d, info = info, error = NULL)
}

# ---------------------------------------------------------------------------
# Bound update
# ---------------------------------------------------------------------------
# Update bounds at the time of analysis by re-running gsDesign() with the
# observed information (n.I) and the planned maximum information
# (maxn.IPlan). This is the documented gsDesign approach to bound updating.
# `fracs` are the information fractions at each analysis of the updated
# design (observed fractions first, then planned future fractions, final = 1).
# `spend_all_final`: if TRUE, the final analysis spends all remaining alpha
# (usTime = 1 at the final analysis); otherwise the spending function is
# evaluated at the actual final information fraction.
update_design_bounds <- function(d, fracs, spend_all_final = FALSE) {
  k_new <- length(fracs)
  max_info <- d$n.I[d$k]
  I_abs <- fracs * max_info
  usTime <- if (isTRUE(spend_all_final)) c(fracs[-k_new], 1) else fracs
  gsDesign(k = k_new, test.type = d$test.type, alpha = d$alpha,
           beta = d$beta, n.fix = d$n.fix, n.I = I_abs,
           maxn.IPlan = max_info, nFixSurv = d$nFixSurv,
           sfu = d$upper$sf, sfupar = d$upper$param,
           sfl = d$lower$sf, sflpar = d$lower$param,
           delta = d$delta, delta0 = d$delta0, delta1 = d$delta1,
           usTime = usTime)
}

# ---------------------------------------------------------------------------
# Design outputs
# ---------------------------------------------------------------------------
# Tabular summary (gsBoundSummary) as a data frame. gsBoundSummary already
# returns a data frame in gsDesign 3.x.
design_table <- function(d, info) {
  logdelta <- inherits(d, "gsSurv") || isTRUE(info$events_mode)
  if (logdelta) {
    gsBoundSummary(d, logdelta = TRUE, ratio = info$ratio,
                   timename = info$time_units, tdigits = 1)
  } else {
    gsBoundSummary(d, ratio = info$ratio, timename = info$time_units,
                   tdigits = 1)
  }
}

# Short textual summary of a design.
design_text <- function(d, info) {
  # d$en is a length-2 vector: expected sample size under H0, then under H1.
  # d$n.I[d$k] is the maximum sample size (or total events for TTE).
  gs_n <- if (inherits(d, "gsSurv")) {
    max(d$N)
  } else if (isTRUE(info$events_mode)) {
    d$nFixSurv
  } else {
    d$n.I[d$k]
  }
  lines <- c(
    paste0("Endpoint type: ", info$endpoint),
    paste0("Number of analyses: ", d$k, " (",
           d$k - 1, " interim ", if (d$k - 1 == 1) "analysis" else "analyses",
           " + final)"),
    paste0("One-sided alpha: ", format(d$alpha, digits = 4),
           "; power: ", format(1 - d$beta, digits = 4)),
    paste0("Test type: ", d$test.type,
           switch(as.character(d$test.type),
             "1" = " (one-sided)",
             "2" = " (two-sided symmetric)",
             "3" = " (two-sided asymmetric, non-binding futility)",
             "4" = " (two-sided asymmetric, binding futility)",
             "5" = " (two-sided asymmetric, beta-spending, binding futility)",
             "6" = " (two-sided asymmetric, beta-spending, non-binding futility)",
             "")),
    paste0("Upper spending function: ", info$sfu_name,
           if (!is.null(d$upper$param)) paste0(" (parameter: ",
             paste(d$upper$param, collapse = ", "), ")") else ""),
    if (!is.null(d$lower$bound) && any(is.finite(d$lower$bound)))
      paste0("Lower spending function: ", info$sfl_name,
             if (!is.null(d$lower$param)) paste0(" (parameter: ",
               paste(d$lower$param, collapse = ", "), ")") else "")
    else "Lower bound: none (efficacy only)",
    paste0("Fixed-design sample size: ", round(info$n_total, 1)),
    paste0("Group-sequential sample size: ", round(gs_n, 1)),
    if (inherits(d, "gsSurv") || isTRUE(info$events_mode))
      paste0("Total events: ", round(d$n.I[d$k], 1))
    else NULL,
    if (inherits(d, "gsSurv") && !is.null(info$gamma_derived))
      paste0("Enrollment rate (adjusted to obtain power): ",
             paste(round(info$gamma_derived, 2), collapse = ", "),
             " per ", info$time_units)
    else NULL,
    if (isTRUE(info$calendar_mode))
      paste0("Calendar times of analyses: ",
             paste(round(info$calendar_times, 1), collapse = ", "),
             " ", info$time_units)
    else NULL,
    if (!is.null(info$study_duration_solved))
      paste0("Study duration (solved): ",
             round(info$study_duration_solved, 1), " ", info$time_units,
             "; enrollment duration extended to ",
             round(info$enroll_duration_solved, 1))
    else NULL,
    paste0("Upper (efficacy) bounds (Z): ",
           paste(round(d$upper$bound, 3), collapse = ", ")),
    if (!is.null(d$lower$bound) && any(is.finite(d$lower$bound)))
      paste0("Lower (futility) bounds (Z): ",
             paste(round(d$lower$bound, 3), collapse = ", "))
    else NULL,
    paste0("Cumulative alpha spent: ",
           paste(round(cumsum(d$upper$spend), 4), collapse = ", ")),
    if (!is.null(d$lower$spend) && any(d$lower$spend > 0))
      paste0("Cumulative beta spent: ",
             paste(round(cumsum(d$lower$spend), 4), collapse = ", "))
    else NULL,
    paste0("Expected sample size under H0: ", round(d$en[1], 1)),
    paste0("Expected sample size under H1: ", round(d$en[2], 1))
  )
  paste(lines[!vapply(lines, is.null, logical(1))], collapse = "\n")
}

# R code that reproduces the design (used by the Report tab).
design_code <- function(p) {
  sf_call <- function(name) {
    if (name == "None") return("NULL")
    switch(name,
      "LDOF" = "sfLDOF", "LDPocock" = "sfLDPocock", "HSD" = "sfHSD",
      "Power" = "sfPower", "Logistic" = "sfLogistic", "Points" = "sfPoints")
  }
  sf_par <- function(name, gamma, rho, nu, points_a) {
    if (name %in% c("HSD", "Power", "Logistic")) {
      switch(name, "HSD" = gamma, "Power" = rho, "Logistic" = nu)
    } else if (name == "Points") {
      sprintf("c(%s)", points_a)
    } else "NULL"
  }
  k <- p$n_interim + 1L
  timing <- if (p$equal_spacing) {
    sprintf("timing = (1:%d) / %d", k, k)
  } else {
    sprintf("timing = c(%s, 1)", p$info_frac)
  }
  test.type <- switch(p$testing,
    "one-sided" = 1L, "two-sided symmetric" = 2L,
    "two-sided asymmetric, non-binding" = 6L,
    "two-sided asymmetric, binding" = 5L, 1L)
  common <- sprintf(paste0(
    "  k = %d, test.type = %d, alpha = %s, beta = %s,\n",
    "  %s,\n",
    "  sfu = %s, sfupar = %s,\n",
    "  sfl = %s, sflpar = %s"),
    k, test.type, p$alpha, 1 - p$power, timing,
    sf_call(p$sfu_name),
    sf_par(p$sfu_name, p$sfu_gamma, p$sfu_rho, p$sfu_nu, p$sfu_points_a),
    if (p$sfl_name == "None") "NULL" else sf_call(p$sfl_name),
    if (p$sfl_name == "None") "NULL" else
      sf_par(p$sfl_name, p$sfl_gamma, p$sfl_rho, p$sfl_nu, p$sfl_points_a))
  # gsSurvCalendar takes no k or timing argument (k comes from the length
  # of calendarTime).
  common_cal <- sprintf(paste0(
    "  test.type = %d, alpha = %s, beta = %s,\n",
    "  sfu = %s, sfupar = %s,\n",
    "  sfl = %s, sflpar = %s"),
    test.type, p$alpha, 1 - p$power,
    sf_call(p$sfu_name),
    sf_par(p$sfu_name, p$sfu_gamma, p$sfu_rho, p$sfu_nu, p$sfu_points_a),
    if (p$sfl_name == "None") "NULL" else sf_call(p$sfl_name),
    if (p$sfl_name == "None") "NULL" else
      sf_par(p$sfl_name, p$sfl_gamma, p$sfl_rho, p$sfl_nu, p$sfl_points_a))
  lambdaC <- control_hazard(p$control_spec, p$control_median,
                            p$control_fail_rate, p$control_fail_time,
                            p$control_surv_rate, p$control_surv_time)
  gamma_R <- function() {
    if (!is.null(p$rampup) && p$rampup > 0) {
      sprintf("gamma = c(%s, %s), R = c(%s, %s)",
              p$enroll_rate / 2, p$enroll_rate, p$rampup, p$enroll_duration)
    } else {
      sprintf("gamma = %s, R = %s", p$enroll_rate, p$enroll_duration)
    }
  }
  body <- switch(p$endpoint,
    "Time-to-event" = {
      if (identical(p$timing_type, "Calendar-based")) {
        sprintf(paste0("x <- gsSurvCalendar(\n%s,\n",
                       "  lambdaC = %s, hr = %s, hr0 = %s, eta = %s,\n",
                       "  %s, minfup = %s, ratio = %s,\n",
                       "  calendarTime = c(%s))"),
                common_cal, format(lambdaC, digits = 6), p$hr, p$hr0, p$eta,
                gamma_R(), p$minfup, p$ratio, p$calendar_times)
      } else if (identical(p$adjust_var, "Number of events")) {
        sprintf(paste0(
          "ss <- nSurvival(lambda1 = %s, lambda2 = %s * %s,\n",
          "  Ts = %s, Tr = %s, eta = %s, ratio = %s,\n",
          "  alpha = %s, beta = %s, sided = %d)\n\n",
          "x <- gsDesign(\n%s,\n",
          "  n.fix = ss$nEvents, nFixSurv = ss$n, delta1 = log(%s))\n",
          "x$hr0 <- %s  # for HR-at-bound calculations in gsBoundSummary()"),
          format(lambdaC, digits = 6), format(lambdaC, digits = 6), p$hr,
          p$study_duration, p$enroll_duration, p$eta, p$ratio,
          p$alpha, 1 - p$power, if (test.type == 1L) 1L else 2L,
          common, p$hr, p$hr0)
      } else if (identical(p$adjust_var, "Study duration")) {
        sprintf(paste0("x <- gsSurv(\n%s,\n",
                       "  lambdaC = %s, hr = %s, hr0 = %s, eta = %s,\n",
                       "  %s, T = NULL, minfup = %s, ratio = %s)"),
                common, format(lambdaC, digits = 6), p$hr, p$hr0, p$eta,
                gamma_R(), p$minfup, p$ratio)
      } else {
        sprintf(paste0("x <- gsSurv(\n%s,\n",
                       "  lambdaC = %s, hr = %s, hr0 = %s, eta = %s,\n",
                       "  %s, T = %s, minfup = %s, ratio = %s)"),
                common, format(lambdaC, digits = 6), p$hr, p$hr0, p$eta,
                gamma_R(), p$study_duration, p$minfup, p$ratio)
      }
    },
    "Binomial" = {
      bf <- binomial_fixed(p$p1, p$p2, p$alpha, 1 - p$power, p$ratio,
                           if (test.type == 1L) 1L else 2L,
                           p$bin_scale, p$delta0)
      sprintf(paste0("x <- gsDesign(\n%s,\n",
                     "  n.fix = %s, delta = %s, delta0 = 0, delta1 = %s)\n",
                     "# fixed-design N from nBinomial(p1 = %s, p2 = %s, ",
                     "scale = \"%s\")"),
              common, format(bf$n, digits = 6), format(bf$delta, digits = 6),
              format(bf$delta, digits = 6), p$p1, p$p2, p$bin_scale)
    },
    "Normal" = {
      nf <- normal_fixed(p$effect, p$sd_c, p$sd_e, p$alpha, 1 - p$power,
                         p$ratio, if (test.type == 1L) 1L else 2L, p$delta0)
      sprintf(paste0("x <- gsDesign(\n%s,\n",
                     "  n.fix = %s, delta = %s, delta0 = 0, delta1 = %s)\n",
                     "# fixed-design N from nNormal(delta1 = %s, sd = %s, ",
                     "sd2 = %s)"),
              common, format(nf$n, digits = 6), format(nf$delta, digits = 6),
              format(nf$delta, digits = 6),
              p$effect, p$sd_c, p$sd_e)
    },
    sprintf(paste0("x <- gsDesign(\n%s,\n",
                   "  n.fix = %s, delta = %s, delta0 = 0, delta1 = %s)"),
            common, p$n_fixed, p$delta_std, p$delta_std)
  )
  paste0(
    "library(gsDesign)\n\n",
    "# Group sequential design reproduced from the Group Sequential Design Explorer\n",
    "# Endpoint type: ", p$endpoint, "\n\n",
    body, "\n\n",
    "summary(x)\n",
    "gsBoundSummary(x)\n"
  )
}

# R code that reproduces a bound update (used by the Update Bounds tab).
# `d` is the base design, `fracs` the information fractions of the updated
# design, `spend_all_final` the final-analysis spending strategy.
update_code <- function(d, fracs, spend_all_final) {
  k_new <- length(fracs)
  paste0(
    "library(gsDesign)\n\n",
    "# Update bounds for the design created on the Create Design tab\n",
    "# (design parameters reproduced from the saved design object)\n\n",
    sprintf("x <- gsDesign(k = %d, test.type = %d, alpha = %s, beta = %s,\n",
            d$k, d$test.type, d$alpha, d$beta),
    sprintf("  timing = c(%s), sfu = %s, sfupar = c(%s),\n",
            paste(d$timing, collapse = ", "),
            sf_fn_name(d$upper$sf),
            paste(d$upper$param, collapse = ", ")),
    sprintf("  sfl = %s, sflpar = c(%s),\n",
            sf_fn_name(d$lower$sf), paste(d$lower$param, collapse = ", ")),
    sprintf("  n.fix = %s, delta = %s, delta0 = %s, delta1 = %s)\n\n",
            d$n.fix, format(d$delta, digits = 6),
            format(d$delta0, digits = 6), format(d$delta1, digits = 6)),
    sprintf("# Observed / planned information fractions: c(%s)\n",
            paste(fracs, collapse = ", ")),
    sprintf("y <- gsDesign(k = %d, test.type = x$test.type, alpha = x$alpha,\n",
            k_new),
    sprintf("  beta = x$beta, n.fix = x$n.fix,\n"),
    sprintf("  n.I = c(%s) * x$n.I[x$k], maxn.IPlan = x$n.I[x$k],\n",
            paste(fracs, collapse = ", ")),
    sprintf("  nFixSurv = x$nFixSurv,\n"),
    sprintf("  sfu = x$upper$sf, sfupar = x$upper$param,\n"),
    sprintf("  sfl = x$lower$sf, sflpar = x$lower$param,\n"),
    sprintf("  delta = x$delta, delta0 = x$delta0, delta1 = x$delta1,\n"),
    sprintf("  usTime = c(%s))\n\n",
            paste(if (isTRUE(spend_all_final)) c(fracs[-k_new], 1) else fracs,
                  collapse = ", ")),
    "summary(y)\n",
    "gsBoundSummary(y)\n"
  )
}

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
ui <- navbarPage(
  title = "Group Sequential Design Explorer",
  id = "main_tabs",

  # ============================ CREATE DESIGN ==============================
  tabPanel(
    "Create Design",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        h4("Design Parameters"),
        tabsetPanel(
          id = "param_tabs", type = "tabs",

          # ---- Sample Size ----
          tabPanel(
            "Sample Size",
            br(),
            selectInput("endpoint", "Endpoint type",
                        choices = c("Time-to-event", "Binomial", "Normal",
                                    "Fixed design N"),
                        selected = "Time-to-event"),
            numericInput("alpha", "alpha (1-sided)", value = 0.025,
                         min = 0.0001, max = 0.5, step = 0.001),
            numericInput("power", "Power", value = 0.9,
                         min = 0.5, max = 0.999, step = 0.01),
            numericInput("ratio",
                         "Randomization ratio (experimental : control)",
                         value = 1, min = 0.1, max = 10, step = 0.1),
            textInput("time_units", "Time units", value = "Months"),
            conditionalPanel(
              "input.endpoint == 'Time-to-event'",
              selectInput("adjust_var",
                          "Variable to adjust to obtain desired power",
                          choices = c("Sample size", "Number of events",
                                      "Study duration"),
                          selected = "Sample size"),
              numericInput("hr", "Hazard ratio (H1)", value = 0.7,
                           min = 0.01, max = 5, step = 0.01),
              numericInput("hr0", "Hazard ratio (H0)", value = 1,
                           min = 0.01, max = 5, step = 0.01),
              numericInput("eta", "Censoring (dropout) rate", value = 0,
                           min = 0, max = 1, step = 0.01)
            ),
            conditionalPanel(
              "input.endpoint == 'Binomial'",
              numericInput("p1", "Control response rate", value = 0.3,
                           min = 0.001, max = 0.999, step = 0.01),
              numericInput("p2", "Experimental response rate", value = 0.45,
                           min = 0.001, max = 0.999, step = 0.01),
              selectInput("bin_scale", "Treatment effect scale",
                          choices = c("Difference" = "Difference",
                                      "Risk ratio" = "RR",
                                      "Odds ratio" = "OR"),
                          selected = "Difference"),
              numericInput("delta0", "H0 difference", value = 0,
                           step = 0.01),
              helpText("For risk ratio and odds ratio, the H0 difference
                       is on the log scale (0 = no difference).")
            ),
            conditionalPanel(
              "input.endpoint == 'Normal'",
              numericInput("effect",
                           "Treatment effect (alternative hypothesis)",
                           value = 1, step = 0.1),
              numericInput("delta0", "H0 difference", value = 0,
                           step = 0.1),
              numericInput("sd_c", "Standard deviation (control)",
                           value = 1, min = 0.01, step = 0.1),
              numericInput("sd_e", "Standard deviation (experimental)",
                           value = 1, min = 0.01, step = 0.1)
            ),
            conditionalPanel(
              "input.endpoint == 'Fixed design N'",
              numericInput("n_fixed", "N for fixed design", value = 100,
                           min = 2, step = 1),
              numericInput("delta_std", "Standardized treatment effect",
                           value = 0.5, step = 0.05)
            )
          ),

          # ---- Timing ----
          tabPanel(
            "Timing",
            br(),
            selectInput("timing_type", "Timing type",
                        choices = c("Information-based", "Calendar-based"),
                        selected = "Information-based"),
            conditionalPanel(
              "input.endpoint == 'Time-to-event' && input.timing_type == 'Calendar-based'",
              textInput("calendar_times",
                        "Calendar times of analyses (comma-separated)",
                        value = "12, 18"),
              helpText("Analyses occur at these calendar times; spending
                       is by information fraction at those times
                       (gsSurvCalendar). The number of interim analyses
                       above is ignored in calendar mode.")
            ),
            numericInput("n_interim", "Number of interim analyses",
                         value = 1, min = 0, max = 9, step = 1),
            checkboxInput("equal_spacing", "Equal information spacing",
                          value = TRUE),
            conditionalPanel(
              "!input.equal_spacing",
              textInput("info_frac",
                        "Information fraction at interim analyses
                        (comma-separated)",
                        value = "0.5"),
              helpText("Enter one value per interim analysis, increasing,
                       between 0 and 1. The final analysis is at 1.")
            ),
            helpText("Interim analyses must occur prior to enrollment
                     duration + minimum follow-up.")
          ),

          # ---- Boundaries ----
          tabPanel(
            "Boundaries",
            br(),
            selectInput("testing", "Testing",
                        choices = c(
                          "One-sided" = "one-sided",
                          "Two-sided symmetric" = "two-sided symmetric",
                          "Two-sided asymmetric, non-binding futility" =
                            "two-sided asymmetric, non-binding",
                          "Two-sided asymmetric, binding futility" =
                            "two-sided asymmetric, binding"),
                        selected = "one-sided"),
            h5("Upper spending"),
            selectInput("sfu_name", "Spending function",
                        choices = SPENDING_CHOICES, selected = "LDOF"),
            conditionalPanel(
              "input.sfu_name == 'HSD'",
              numericInput("sfu_gamma", "Gamma", value = -4, step = 0.5)
            ),
            conditionalPanel(
              "input.sfu_name == 'Power'",
              numericInput("sfu_rho", "Rho", value = 1, min = 0, step = 0.1)
            ),
            conditionalPanel(
              "input.sfu_name == 'Logistic'",
              numericInput("sfu_nu", "Nu", value = 0.5, step = 0.1)
            ),
            conditionalPanel(
              "input.sfu_name == 'Points'",
              textInput("sfu_points_a",
                        "Cumulative alpha-spend at interim analyses
                        (comma-separated proportions)",
                        value = "0.4"),
              helpText("One proportion per interim analysis, increasing,
                       between 0 and 1. The final analysis spends the
                       remaining alpha.")
            ),
            h5("Lower spending"),
            selectInput("sfl_name", "Spending function",
                        choices = c("None" = "None", SPENDING_CHOICES),
                        selected = "None"),
            conditionalPanel(
              "input.sfl_name == 'HSD'",
              numericInput("sfl_gamma", "Gamma", value = -4, step = 0.5)
            ),
            conditionalPanel(
              "input.sfl_name == 'Power'",
              numericInput("sfl_rho", "Rho", value = 1, min = 0, step = 0.1)
            ),
            conditionalPanel(
              "input.sfl_name == 'Logistic'",
              numericInput("sfl_nu", "Nu", value = 0.5, step = 0.1)
            ),
            conditionalPanel(
              "input.sfl_name == 'Points'",
              textInput("sfl_points_a",
                        "Cumulative beta-spend at interim analyses
                        (comma-separated proportions)",
                        value = "0.4"),
              helpText("One proportion per interim analysis, increasing,
                       between 0 and 1.")
            ),
            hr(),
            h5("Spending summary"),
            verbatimTextOutput("spending_summary")
          ),

          # ---- Hazard Ratios (TTE) ----
          tabPanel(
            "Hazard Ratios",
            br(),
            conditionalPanel(
              "input.endpoint == 'Time-to-event'",
              numericInput("hr2", "Hazard ratio (H1)", value = 0.7,
                           min = 0.01, max = 5, step = 0.01),
              numericInput("hr02", "Hazard ratio (H0)", value = 1,
                           min = 0.01, max = 5, step = 0.01),
              helpText("These mirror the hazard ratios on the Sample Size
                       tab.")
            ),
            conditionalPanel(
              "input.endpoint != 'Time-to-event'",
              helpText("Hazard ratios apply to time-to-event designs only.")
            )
          ),

          # ---- Failure Rates (TTE) ----
          tabPanel(
            "Failure Rates",
            br(),
            conditionalPanel(
              "input.endpoint == 'Time-to-event'",
              selectInput("control_spec",
                          "Control failure rate specification",
                          choices = c("Control median" = "median",
                                      "Control failure rate" = "failrate",
                                      "Control survival rate" = "survrate"),
                          selected = "median"),
              conditionalPanel(
                "input.control_spec == 'median'",
                numericInput("control_median", "Control median",
                             value = 12, min = 0.1, step = 0.5)
              ),
              conditionalPanel(
                "input.control_spec == 'failrate'",
                numericInput("control_fail_rate", "Control failure rate",
                             value = 0.2, min = 0.001, max = 0.999,
                             step = 0.01),
                numericInput("control_fail_time", "Time", value = 12,
                             min = 0.1, step = 0.5)
              ),
              conditionalPanel(
                "input.control_spec == 'survrate'",
                numericInput("control_surv_rate", "Control survival rate",
                             value = 0.5, min = 0.001, max = 0.999,
                             step = 0.01),
                numericInput("control_surv_time", "Time", value = 12,
                             min = 0.1, step = 0.5)
              ),
              numericInput("eta2", "Censoring (dropout) rate", value = 0,
                           min = 0, max = 1, step = 0.01),
              helpText("Failure rates are assumed exponential. The control
                       hazard is derived from the specification above.")
            ),
            conditionalPanel(
              "input.endpoint != 'Time-to-event'",
              helpText("Failure rates apply to time-to-event designs only.")
            )
          ),

          # ---- Enrollment (TTE) ----
          tabPanel(
            "Enrollment",
            br(),
            conditionalPanel(
              "input.endpoint == 'Time-to-event'",
              numericInput("enroll_rate", "Steady-state enrollment rate",
                           value = 10, min = 0.01, step = 1),
              numericInput("enroll_duration", "Enrollment duration",
                           value = 12, min = 0.1, step = 0.5),
              numericInput("minfup", "Minimum follow-up duration",
                           value = 12, min = 0, step = 0.5),
              numericInput("rampup", "Enrollment ramp-up duration",
                           value = 0, min = 0, step = 0.5),
              numericInput("study_duration", "Maximum time (study duration)",
                           value = 24, min = 0.1, step = 0.5),
              helpText("When 'Sample size' is the variable to adjust,
                       gsSurv scales the enrollment rate to obtain the
                       desired power with the given enrollment duration.
                       When 'Study duration' is the variable to adjust,
                       the study duration is solved for with the
                       enrollment rate held fixed. During ramp-up,
                       enrollment runs at half the steady-state rate.")
            ),
            conditionalPanel(
              "input.endpoint != 'Time-to-event'",
              helpText("Enrollment applies to time-to-event designs only.")
            )
          )
        )
      ),

      mainPanel(
        width = 8,
        h4("Design Output"),
        tabsetPanel(
          id = "output_tabs", type = "tabs",
          tabPanel(
            "Plot",
            br(),
            selectInput("plot_type", "Plot type",
                        choices = c(
                          "Efficacy and futility bounds (Z-scale)" = "bounds_z",
                          "Nominal p-value bounds" = "bounds_p",
                          "Cumulative power" = "power",
                          "Enrollment and events over time" = "enroll"),
                        selected = "bounds_z"),
            plotOutput("design_plot", height = "520px")
          ),
          tabPanel(
            "Tabular",
            br(),
            DT::dataTableOutput("design_table"),
            br(),
            downloadButton("download_csv", "Download table (CSV)")
          ),
          tabPanel(
            "Text",
            br(),
            verbatimTextOutput("design_text")
          ),
          tabPanel(
            "Report",
            br(),
            p("The R code below reproduces this design with the gsDesign
              package. Download it as an R Markdown file to re-run later."),
            verbatimTextOutput("design_code"),
            br(),
            downloadButton("download_rmd", "Download R Markdown (.Rmd)"),
            downloadButton("download_html", "Download HTML report")
          )
        )
      )
    )
  ),

  # ============================ UPDATE BOUNDS ==============================
  tabPanel(
    "Update Bounds",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        h4("Bound Update Parameters"),
        checkboxInput("use_current_design",
                      "Use the design from the Create Design tab",
                      value = TRUE),
        conditionalPanel(
          "!input.use_current_design",
          fileInput("update_design_file",
                    "Upload a saved design (.rds)",
                    accept = ".rds")
        ),
        numericInput("ub_n_analyses", "Total number of analyses",
                     value = 3, min = 1, max = 10, step = 1),
        numericInput("ub_done", "Analyses completed so far",
                     value = 1, min = 1, max = 9, step = 1),
        textInput("ub_info_frac",
                  "Observed information fractions at completed analyses
                  (comma-separated)",
                  value = "0.5"),
        selectInput("ub_final_strategy",
                    "Final analysis alpha spending strategy",
                    choices = c(
                      "Use spending function at final analysis" = "spending",
                      "Spend all remaining alpha at final analysis" =
                        "spend-all"),
                    selected = "spending"),
        helpText("Bounds are updated by re-running gsDesign() with the
                 observed information (n.I) and the planned maximum
                 information (maxn.IPlan). Future interim analyses, if
                 any, are placed at equal information spacing between the
                 last observed fraction and 1. The final analysis is
                 always at information fraction 1.")
      ),
      mainPanel(
        width = 8,
        h4("Bound Update Output"),
        tabsetPanel(
          id = "ub_output_tabs", type = "tabs",
          tabPanel(
            "Tabular",
            br(),
            DT::dataTableOutput("ub_table"),
            br(),
            downloadButton("ub_download_csv", "Download table (CSV)")
          ),
          tabPanel(
            "Report",
            br(),
            verbatimTextOutput("ub_code"),
            br(),
            downloadButton("ub_download_rmd", "Download R Markdown (.Rmd)")
          )
        )
      )
    )
  ),

  # ============================ SAVE / RESTORE =============================
  tabPanel(
    "Save / Restore",
    fluidRow(
      column(
        6,
        h4("Save design"),
        p("Download the current design (parameters and the computed
          gsDesign object) to restore it later."),
        downloadButton("download_rds", "Download design (.rds)"),
        br(), br(),
        downloadButton("download_json",
                       "Download parameters (.json)")
      ),
      column(
        6,
        h4("Restore design"),
        p("Upload a design file saved from this app. All design
          parameters are restored into the Create Design tab."),
        fileInput("restore_file", "Choose a design file (.rds)",
                  accept = ".rds"),
        verbatimTextOutput("restore_status")
      )
    )
  ),

  # ================================ ABOUT ==================================
  tabPanel(
    "About",
    fluidRow(
      column(
        8, offset = 1,
        h3("About this app"),
        p("This is a web interface to the open-source R package gsDesign
          for designing group sequential clinical trials. It was built
          independently with Shiny and calls the gsDesign package through
          its public API; it contains no code copied from any other web
          application."),
        h4("Features"),
        tags$ul(
          tags$li(strong("Create design: "), "designs for time-to-event,
                  binomial, and normal endpoints, as well as extending a
                  fixed-design sample size to a group sequential design."),
          tags$li("Information-based design with flexible timing and any
                  number of interim analyses."),
          tags$li("Spending-function boundaries: O'Brien-Fleming, Pocock,
                  Lan-DeMets, Hwang-Shih-DeCani, power family, logistic,
                  and custom spending points, for efficacy and futility."),
          tags$li("For time-to-event (proportional hazards) designs,
                  enrollment rates, dropout rates, and event rates can be
                  specified; expected calendar times of analyses are
                  shown."),
          tags$li("Tabular output that can be copied into a word
                  processor, a brief textual summary, and several plots."),
          tags$li("An R Markdown report generator that saves gsDesign
                  code which can be re-run to reproduce the design."),
          tags$li(strong("Update bounds: "), "update bounds at the time
                  of analysis, with a flexible number of total analyses
                  and alpha-spending strategies for the final analysis."),
          tags$li(strong("Save and restore design: "), "save and restore
                  all parameters when creating a design or updating
                  bounds.")
        ),
        h4("Contributors"),
        p("gsDesign was created by Keaven Anderson, PhD, and is
          maintained by Nan Xiao, PhD. This web interface is an
          independent implementation built on the gsDesign package."),
        h4("License"),
        p("The gsDesign R package is licensed under GPL (>= 3). This
          application is provided for designing group sequential trials
          and carries no warranty; all calculations should be verified
          before use in a real trial."),
        h4("References"),
        tags$ul(
          tags$li(a("gsDesign on CRAN",
                    href = "https://cran.r-project.org/package=gsDesign")),
          tags$li(a("gsDesign documentation",
                    href = "https://keaven.github.io/gsDesign/")),
          tags$li("Jennison, C. and Turnbull, B. W. (2000). Group
                  Sequential Methods with Applications to Clinical Trials.
                  Chapman & Hall/CRC.")
        )
      )
    )
  )
)

# ---------------------------------------------------------------------------
# Input registry (for save / restore)
# ---------------------------------------------------------------------------
INPUT_TYPES <- list(
  endpoint = "select", alpha = "numeric", power = "numeric",
  ratio = "numeric", time_units = "text", testing = "select",
  n_interim = "numeric", equal_spacing = "checkbox",
  info_frac = "text", timing_type = "select", calendar_times = "text",
  sfu_name = "select", sfu_gamma = "numeric", sfu_rho = "numeric",
  sfu_nu = "numeric", sfu_points_a = "text",
  sfl_name = "select", sfl_gamma = "numeric", sfl_rho = "numeric",
  sfl_nu = "numeric", sfl_points_a = "text",
  adjust_var = "select", hr = "numeric", hr0 = "numeric", eta = "numeric",
  control_spec = "select", control_median = "numeric",
  control_fail_rate = "numeric", control_fail_time = "numeric",
  control_surv_rate = "numeric", control_surv_time = "numeric",
  study_duration = "numeric", minfup = "numeric",
  enroll_rate = "numeric", enroll_duration = "numeric", rampup = "numeric",
  p1 = "numeric", p2 = "numeric", bin_scale = "select", delta0 = "numeric",
  effect = "numeric", sd_c = "numeric", sd_e = "numeric",
  n_fixed = "numeric", delta_std = "numeric"
)

collect_params <- function(input) {
  p <- lapply(names(INPUT_TYPES), function(id) input[[id]])
  names(p) <- names(INPUT_TYPES)
  p
}

# Information fractions for the updated design: observed fractions for the
# completed analyses, then (if the total number of analyses exceeds the
# completed ones plus the final) equally spaced future interim fractions
# between the last observed fraction and 1, then the final analysis at 1.
update_fracs <- function(obs, m, k_new) {
  if (k_new < m + 1L) return(NULL)
  if (k_new == m + 1L) return(c(obs, 1))
  n_future <- k_new - m - 1L
  future <- obs[m] + (1 - obs[m]) * seq_len(n_future) / (n_future + 1L)
  c(obs, future, 1)
}

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  # --- keep mirrored hazard-ratio / dropout inputs in sync ---------------
  observeEvent(input$hr, {
    if (!is.null(input$hr2) && !isTRUE(all.equal(input$hr2, input$hr)))
      updateNumericInput(session, "hr2", value = input$hr)
  })
  observeEvent(input$hr2, {
    if (!is.null(input$hr) && !isTRUE(all.equal(input$hr, input$hr2)))
      updateNumericInput(session, "hr", value = input$hr2)
  })
  observeEvent(input$hr0, {
    if (!is.null(input$hr02) && !isTRUE(all.equal(input$hr02, input$hr0)))
      updateNumericInput(session, "hr02", value = input$hr0)
  })
  observeEvent(input$hr02, {
    if (!is.null(input$hr0) && !isTRUE(all.equal(input$hr0, input$hr02)))
      updateNumericInput(session, "hr0", value = input$hr02)
  })
  observeEvent(input$eta, {
    if (!is.null(input$eta2) && !isTRUE(all.equal(input$eta2, input$eta)))
      updateNumericInput(session, "eta2", value = input$eta)
  })
  observeEvent(input$eta2, {
    if (!is.null(input$eta) && !isTRUE(all.equal(input$eta, input$eta2)))
      updateNumericInput(session, "eta", value = input$eta2)
  })

  # --- current design -----------------------------------------------------
  params <- reactive(collect_params(input))

  built <- reactive({
    p <- params()
    req(p$alpha, p$power, p$endpoint)
    build_design(p)
  })

  design_obj <- reactive({
    b <- built()
    if (!is.null(b$error)) {
      showNotification(b$error, type = "error", duration = 8)
      return(NULL)
    }
    b$design
  })

  design_info <- reactive({
    b <- built()
    b$info
  })

  # --- spending summary ---------------------------------------------------
  output$spending_summary <- renderText({
    d <- design_obj()
    if (is.null(d)) return("No design yet.")
    paste0(
      "Total alpha: ", format(d$alpha, digits = 4), "\n",
      "Alpha spent at each analysis: ",
      paste(round(d$upper$spend, 4), collapse = ", "), "\n",
      "Cumulative alpha spent: ",
      paste(round(cumsum(d$upper$spend), 4), collapse = ", "), "\n",
      if (!is.null(d$lower$spend) && any(d$lower$spend > 0))
        paste0("Beta spent at each analysis: ",
               paste(round(d$lower$spend, 4), collapse = ", "), "\n",
               "Cumulative beta spent: ",
               paste(round(cumsum(d$lower$spend), 4), collapse = ", "))
      else "Lower spending: none"
    )
  })

  # --- plots --------------------------------------------------------------
  output$design_plot <- renderPlot({
    d <- design_obj()
    info <- design_info()
    req(d, info)
    k <- d$k
    timing <- d$timing
    if (input$plot_type == "bounds_z") {
      df <- data.frame(
        Analysis = rep(1:k, 2),
        Information = rep(timing, 2),
        Z = c(d$upper$bound,
              if (!is.null(d$lower$bound)) d$lower$bound else rep(NA, k)),
        Bound = rep(c("Efficacy (upper)", "Futility (lower)"),
                    each = k)
      )
      df <- df[is.finite(df$Z), ]
      ggplot(df, aes(x = Information, y = Z, color = Bound)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        geom_line(linewidth = 1) + geom_point(size = 3) +
        geom_text(aes(label = round(Z, 2)), vjust = -0.9, size = 3.5,
                  show.legend = FALSE) +
        scale_x_continuous(breaks = timing, labels = round(timing, 2)) +
        labs(title = "Group sequential boundaries (Z-scale)",
             x = "Information fraction", y = "Critical value (Z)") +
        theme_minimal(base_size = 13) + theme(legend.position = "bottom")
    } else if (input$plot_type == "bounds_p") {
      df <- data.frame(
        Analysis = rep(1:k, 2),
        Information = rep(timing, 2),
        p = c(1 - pnorm(d$upper$bound),
              if (!is.null(d$lower$bound)) pnorm(d$lower$bound)
              else rep(NA, k)),
        Bound = rep(c("Efficacy (upper)", "Futility (lower)"), each = k)
      )
      df <- df[is.finite(df$p), ]
      ggplot(df, aes(x = Information, y = p, color = Bound)) +
        geom_line(linewidth = 1) + geom_point(size = 3) +
        geom_text(aes(label = signif(p, 3)), vjust = -0.9, size = 3.5,
                  show.legend = FALSE) +
        scale_x_continuous(breaks = timing, labels = round(timing, 2)) +
        scale_y_log10() +
        labs(title = "Nominal p-value bounds (one-sided)",
             x = "Information fraction", y = "Nominal p-value (log scale)") +
        theme_minimal(base_size = 13) + theme(legend.position = "bottom")
    } else if (input$plot_type == "power") {
      # gsProbability gives per-analysis crossing probabilities; the
      # cumulative curve is their running sum.
      cum1 <- cumsum(crossing_probs(d, d$delta))
      cum0 <- cumsum(crossing_probs(d, 0))
      df <- data.frame(
        Analysis = rep(1:k, 2),
        Information = rep(timing, 2),
        Prob = c(cum1, cum0),
        Scenario = rep(c("Under H1 (power)", "Under H0 (type I error)"),
                       each = k)
      )
      ggplot(df, aes(x = Information, y = Prob, color = Scenario)) +
        geom_line(linewidth = 1) + geom_point(size = 3) +
        geom_text(aes(label = round(Prob, 3)), vjust = -0.9, size = 3.5,
                  show.legend = FALSE) +
        scale_x_continuous(breaks = timing, labels = round(timing, 2)) +
        ylim(0, 1) +
        labs(title = "Cumulative probability of crossing the efficacy bound",
             x = "Information fraction", y = "Cumulative probability") +
        theme_minimal(base_size = 13) + theme(legend.position = "bottom")
    } else { # enrollment / events over time (TTE)
      if (!inherits(d, "gsSurv")) {
        plot.new()
        title("Enrollment and events plot is available for time-to-event designs only.")
        return(invisible())
      }
      tt <- d$T
      ev <- d$eDC + d$eDE
      df <- data.frame(Time = tt, Events = ev)
      ggplot(df, aes(x = Time, y = Events)) +
        geom_line(linewidth = 1, color = "steelblue") +
        geom_point(size = 3, color = "steelblue") +
        geom_text(aes(label = round(Events, 1)), vjust = -0.9, size = 3.5) +
        labs(title = "Expected cumulative events at each analysis",
             x = paste0("Calendar time (", info$time_units, ")"),
             y = "Expected events") +
        theme_minimal(base_size = 13)
    }
  })

  # --- tabular output -----------------------------------------------------
  table_df <- reactive({
    d <- design_obj()
    info <- design_info()
    req(d, info)
    as.data.frame(design_table(d, info))
  })

  output$design_table <- DT::renderDataTable({
    DT::datatable(table_df(), rownames = FALSE,
                  options = list(pageLength = 25, dom = "t",
                                 scrollX = TRUE))
  })

  output$download_csv <- downloadHandler(
    filename = function() "gsdesign_summary.csv",
    content = function(file) write.csv(table_df(), file, row.names = FALSE)
  )

  # --- text output --------------------------------------------------------
  output$design_text <- renderText({
    d <- design_obj()
    info <- design_info()
    req(d, info)
    design_text(d, info)
  })

  # --- report -------------------------------------------------------------
  code_text <- reactive({
    p <- params()
    req(p$endpoint)
    design_code(p)
  })

  output$design_code <- renderText(code_text())

  output$download_rmd <- downloadHandler(
    filename = function() "gsdesign_report.Rmd",
    content = function(file) {
      rmd <- paste0(
        "---\ntitle: \"Group Sequential Design Report\"\n",
        "output: html_document\n---\n\n",
        "```{r setup, include=FALSE}\n",
        "knitr::opts_chunk$set(echo = TRUE)\n```\n\n",
        "## Design\n\n",
        "```{r}\n", code_text(), "```\n"
      )
      writeLines(rmd, file)
    }
  )

  output$download_html <- downloadHandler(
    filename = function() "gsdesign_report.html",
    content = function(file) {
      if (!requireNamespace("rmarkdown", quietly = TRUE)) {
        stop("The rmarkdown package is not installed on the server; ",
             "download the .Rmd file instead.")
      }
      rmd <- paste0(
        "---\ntitle: \"Group Sequential Design Report\"\n",
        "output: html_document\n---\n\n",
        "```{r setup, include=FALSE}\n",
        "knitr::opts_chunk$set(echo = TRUE)\n```\n\n",
        "## Design\n\n",
        "```{r}\n", code_text(), "```\n"
      )
      tmp <- tempfile(fileext = ".Rmd")
      writeLines(rmd, tmp)
      out <- rmarkdown::render(tmp, quiet = TRUE)
      file.copy(out, file, overwrite = TRUE)
    }
  )

  # --- update bounds ------------------------------------------------------
  base_design <- reactive({
    if (isTRUE(input$use_current_design)) {
      list(design = design_obj(), info = design_info())
    } else {
      req(input$update_design_file)
      obj <- readRDS(input$update_design_file$datapath)
      list(design = obj$design, info = obj$info)
    }
  })

  # Parsed, validated update inputs: list(fracs, spend_all) or NULL.
  update_inputs <- reactive({
    obs <- parse_num_vector(input$ub_info_frac)
    obs <- obs[!is.na(obs)]
    m <- as.integer(input$ub_done)
    k_new <- as.integer(input$ub_n_analyses)
    if (length(obs) != m || any(obs <= 0) || any(obs >= 1) ||
        is.unsorted(obs)) {
      showNotification(paste0("Observed information fractions must be ", m,
                              " increasing values between 0 and 1."),
                       type = "error", duration = 8)
      return(NULL)
    }
    fracs <- update_fracs(obs, m, k_new)
    if (is.null(fracs)) {
      showNotification("Total number of analyses must be at least the
                       number completed plus one (the final analysis).",
                       type = "error", duration = 8)
      return(NULL)
    }
    list(fracs = fracs,
         spend_all = identical(input$ub_final_strategy, "spend-all"))
  })

  updated_bounds <- reactive({
    bd <- base_design()
    req(bd$design)
    ui <- update_inputs()
    req(ui)
    tryCatch(
      update_design_bounds(bd$design, ui$fracs, ui$spend_all),
      error = function(e) {
        showNotification(paste("Bound update failed:",
                               conditionMessage(e)),
                         type = "error", duration = 8)
        NULL
      }
    )
  })

  ub_table_df <- reactive({
    u <- updated_bounds()
    req(u)
    bd <- base_design()
    logdelta <- inherits(bd$design, "gsSurv") ||
      isTRUE(bd$info$events_mode)
    as.data.frame(gsBoundSummary(u, logdelta = logdelta,
                                 ratio = bd$info$ratio,
                                 timename = bd$info$time_units, tdigits = 1))
  })

  output$ub_table <- DT::renderDataTable({
    DT::datatable(ub_table_df(), rownames = FALSE,
                  options = list(pageLength = 25, dom = "t",
                                 scrollX = TRUE))
  })

  output$ub_download_csv <- downloadHandler(
    filename = function() "gsdesign_updated_bounds.csv",
    content = function(file) write.csv(ub_table_df(), file, row.names = FALSE)
  )

  ub_code_text <- reactive({
    bd <- base_design()
    req(bd$design)
    ui <- update_inputs()
    req(ui)
    update_code(bd$design, ui$fracs, ui$spend_all)
  })

  output$ub_code <- renderText(ub_code_text())

  output$ub_download_rmd <- downloadHandler(
    filename = function() "gsdesign_bound_update.Rmd",
    content = function(file) {
      rmd <- paste0(
        "---\ntitle: \"Group Sequential Bound Update\"\n",
        "output: html_document\n---\n\n",
        "```{r setup, include=FALSE}\n",
        "knitr::opts_chunk$set(echo = TRUE)\n```\n\n",
        "## Bound update\n\n",
        "```{r}\n",
        isolate(ub_code_text()),
        "```\n"
      )
      writeLines(rmd, file)
    }
  )

  # --- save / restore -----------------------------------------------------
  output$download_rds <- downloadHandler(
    filename = function() "gsdesign_design.rds",
    content = function(file) {
      saveRDS(list(params = params(), design = design_obj(),
                   info = design_info()), file)
    }
  )

  output$download_json <- downloadHandler(
    filename = function() "gsdesign_params.json",
    content = function(file) {
      writeLines(jsonlite::toJSON(params(), auto_unbox = TRUE, pretty = TRUE),
                 file)
    }
  )

  restore_msg <- reactiveVal("No file restored yet.")
  output$restore_status <- renderText(restore_msg())

  observeEvent(input$restore_file, {
    req(input$restore_file)
    obj <- tryCatch(readRDS(input$restore_file$datapath),
                    error = function(e) e)
    if (inherits(obj, "error")) {
      restore_msg(paste("Could not read file:", conditionMessage(obj)))
      return()
    }
    if (is.null(obj$params)) {
      restore_msg("This file does not contain saved design parameters.")
      return()
    }
    for (id in names(obj$params)) {
      type <- INPUT_TYPES[[id]]
      val <- obj$params[[id]]
      if (is.null(type) || is.null(val)) next
      # Only update inputs that exist in the current session.
      if (is.null(input[[id]])) next
      tryCatch({
        if (type == "numeric") {
          updateNumericInput(session, id, value = val)
        } else if (type == "text") {
          updateTextInput(session, id, value = val)
        } else if (type == "select") {
          updateSelectInput(session, id, selected = val)
        } else if (type == "checkbox") {
          updateCheckboxInput(session, id, value = isTRUE(val))
        }
      }, error = function(e) NULL)
    }
    # Keep the mirrored hazard-ratio / dropout inputs in sync as well.
    if (!is.null(obj$params$hr)) updateNumericInput(session, "hr2", value = obj$params$hr)
    if (!is.null(obj$params$hr0)) updateNumericInput(session, "hr02", value = obj$params$hr0)
    if (!is.null(obj$params$eta)) updateNumericInput(session, "eta2", value = obj$params$eta)
    restore_msg(paste0("Restored design from ", input$restore_file$name,
                       " (endpoint: ", obj$params$endpoint, ")."))
  })
}

# ---------------------------------------------------------------------------
shinyApp(ui = ui, server = server)
