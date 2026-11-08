# Validation script for the Group Sequential Design Explorer (v3).
# Run: Rscript test_validation.R
# Validates gsDesign 3.11.0 API usage and the app's own helper functions
# (build_design, design_power, design_table, design_text, design_code,
# update_design_bounds) by sourcing app.R (shinyApp() is not launched).
suppressPackageStartupMessages(library(gsDesign))

ok <- function(name, cond) cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", name))

cat("=== 1. Two-look O'Brien-Fleming (Lan-DeMets) design ===\n")
# Classic 2-look OF bounds at information fractions 0.5 / 1.
d2 <- gsDesign(k = 2, test.type = 1, alpha = 0.025, beta = 0.1,
               sfu = sfLDOF, timing = c(0.5, 1))
cat("LDOF upper bounds:", d2$upper$bound, "\n")
ok("LDOF 2-look bounds = 2.9626 / 1.9686",
   all(abs(d2$upper$bound - c(2.9626, 1.9686)) < 0.001))
# HSD(gamma = -4) is often quoted as "approximating" OF but is NOT the same.
d2h <- gsDesign(k = 2, test.type = 1, alpha = 0.025, beta = 0.1,
                sfu = sfHSD, sfupar = -4, timing = c(0.5, 1))
cat("HSD(-4) upper bounds:", d2h$upper$bound, " (not OF)\n")

cat("\n=== 2. gsSurv documented example ===\n")
xgs <- gsSurv(lambdaC = 0.2, hr = 0.5, eta = 0.1, T = 2, minfup = 1.5)
cat("N =", xgs$N, " events (n.I) =", xgs$n.I, "\n")
ok("gsSurv example N ~ 460, events ~ 97",
   abs(max(xgs$N) - 460) < 5 && abs(xgs$n.I[xgs$k] - 97) < 5)

cat("\n=== 3. gsSurv gamma scaling / T=NULL / piecewise enrollment ===\n")
a1 <- gsSurv(lambdaC = log(2) / 12, hr = 0.7, eta = 0, T = 24, minfup = 12,
             gamma = 10, R = 12)
cat("gamma=10,R=12: gamma derived =", a1$gamma, " N =", max(a1$N),
    " events =", a1$n.I[a1$k], "\n")
ok("gsSurv scales gamma up to reach power", a1$gamma > 10)
a4 <- gsSurv(lambdaC = log(2) / 12, hr = 0.7, eta = 0, T = NULL, minfup = 12,
             gamma = 10, R = 12)
cat("T=NULL: gamma =", a4$gamma, " R =", a4$R, " T =", a4$T, "\n")
ok("T=NULL keeps gamma fixed at 10", abs(a4$gamma - 10) < 1e-9)
a5 <- gsSurv(lambdaC = log(2) / 12, hr = 0.7, eta = 0, T = 24, minfup = 12,
             gamma = c(5, 10), R = c(2, 10))
cat("piecewise gamma=c(5,10), R=c(2,10): N =", max(a5$N),
    " events =", a5$n.I[a5$k], "\n")
ok("piecewise enrollment works", is.finite(max(a5$N)))

cat("\n=== 4. nNormal / nBinomial return numeric scalars ===\n")
nn <- nNormal(delta1 = 1, sd = 1, alpha = 0.025, beta = 0.1, ratio = 1, sided = 1)
cat("nNormal =", nn, " (hand calc total 42.04)\n")
ok("nNormal total N ~ 42", abs(as.numeric(nn) - 42) < 2)
nb <- nBinomial(p1 = 0.3, p2 = 0.45, alpha = 0.025, beta = 0.1,
                ratio = 1, sided = 1, scale = "Difference")
cat("nBinomial(Difference) =", nb, "\n")
ok("nBinomial total N ~ 434", abs(as.numeric(nb) - 434) < 5)

cat("\n=== 5. Source app.R and test build_design across endpoints ===\n")
source("app.R")  # defines helpers + ui/server; does not launch the app

# Default parameter list mirroring the UI defaults.
base_p <- list(
  endpoint = "Time-to-event", alpha = 0.025, power = 0.9, ratio = 1,
  time_units = "Months", testing = "one-sided",
  n_interim = 1L, equal_spacing = TRUE, info_frac = "0.5",
  timing_type = "Information-based", calendar_times = "12, 18",
  sfu_name = "LDOF", sfu_gamma = -4, sfu_rho = 1, sfu_nu = 0.5,
  sfu_points_a = "0.4",
  sfl_name = "None", sfl_gamma = -4, sfl_rho = 1, sfl_nu = 0.5,
  sfl_points_a = "0.4",
  adjust_var = "Sample size", hr = 0.7, hr0 = 1, eta = 0,
  control_spec = "median", control_median = 12,
  control_fail_rate = 0.2, control_fail_time = 12,
  control_surv_rate = 0.5, control_surv_time = 12,
  study_duration = 24, minfup = 12,
  enroll_rate = 10, enroll_duration = 12, rampup = 0,
  p1 = 0.3, p2 = 0.45, bin_scale = "Difference", delta0 = 0,
  effect = 1, sd_c = 1, sd_e = 1,
  n_fixed = 100, delta_std = 0.5
)

check_design <- function(label, p, target_power = 0.9) {
  b <- build_design(p)
  if (!is.null(b$error)) {
    cat(label, "ERROR:", b$error, "\n")
    ok(paste(label, "builds without error"), FALSE)
    return(NULL)
  }
  d <- b$design
  pwr <- design_power(d)
  cat(label, "power =", round(pwr, 4),
      " N =", round(if (inherits(d, "gsSurv")) max(d$N) else d$n.I[d$k], 1),
      " bounds =", paste(round(d$upper$bound, 3), collapse = ", "), "\n")
  ok(paste(label, "power within 0.01 of", target_power),
     abs(pwr - target_power) < 0.01)
  # Tabular and text output must run without error.
  tb <- tryCatch(as.data.frame(design_table(d, b$info)), error = function(e) e)
  ok(paste(label, "design_table runs"), !inherits(tb, "error"))
  tx <- tryCatch(design_text(d, b$info), error = function(e) e)
  ok(paste(label, "design_text runs"), !inherits(tx, "error"))
  # Generated R code must parse and reproduce the design.
  code <- design_code(p)
  env <- new.env()
  ev <- tryCatch({ eval(parse(text = code), envir = env); env$x },
                 error = function(e) e)
  if (inherits(ev, "error")) {
    cat(label, "generated code ERROR:", conditionMessage(ev), "\n")
    ok(paste(label, "generated code runs"), FALSE)
  } else {
    ok(paste(label, "generated code reproduces bounds"),
       all(abs(ev$upper$bound - d$upper$bound) < 1e-6))
  }
  b
}

b_tte <- check_design("TTE sample-size (default)", base_p)

p_bin <- base_p; p_bin$endpoint <- "Binomial"
b_bin <- check_design("Binomial difference", p_bin)
p_bin$bin_scale <- "RR"
invisible(check_design("Binomial risk ratio", p_bin))
p_bin$bin_scale <- "OR"
invisible(check_design("Binomial odds ratio", p_bin))

p_norm <- base_p; p_norm$endpoint <- "Normal"
invisible(check_design("Normal", p_norm))

p_fix <- base_p; p_fix$endpoint <- "Fixed design N"
invisible(check_design("Fixed design N", p_fix))

cat("\n=== 6. TTE adjust modes and calendar timing ===\n")
p_ev <- base_p; p_ev$adjust_var <- "Number of events"
b_ev <- check_design("TTE number of events", p_ev)
if (!is.null(b_ev)) {
  cat("  events =", round(b_ev$info$n_events, 1),
      " N =", round(b_ev$info$n_total, 1), "\n")
}
p_dur <- base_p; p_dur$adjust_var <- "Study duration"
b_dur <- check_design("TTE study duration", p_dur)
if (!is.null(b_dur)) {
  cat("  solved study duration =", round(b_dur$info$study_duration_solved, 2),
      " months\n")
}
p_cal <- base_p; p_cal$timing_type <- "Calendar-based"
b_cal <- check_design("TTE calendar timing", p_cal)
if (!is.null(b_cal)) {
  cat("  calendar times =", paste(round(b_cal$info$calendar_times, 1),
                                  collapse = ", "), "\n")
}
p_ramp <- base_p; p_ramp$rampup <- 2; p_ramp$enroll_duration <- 10
invisible(check_design("TTE ramp-up 2 months", p_ramp))

cat("\n=== 7. Spending functions through build_design ===\n")
p_hsd <- base_p; p_hsd$sfu_name <- "HSD"; p_hsd$sfu_gamma <- -4
invisible(check_design("HSD(-4) upper", p_hsd))
p_poc <- base_p; p_poc$sfu_name <- "LDPocock"
invisible(check_design("Pocock (Lan-DeMets) upper", p_poc))
p_pts <- base_p; p_pts$sfu_name <- "Points"; p_pts$sfu_points_a <- "0.4"
invisible(check_design("Custom points upper", p_pts))
p_beta <- base_p; p_beta$testing <- "two-sided asymmetric, non-binding"
p_beta$sfl_name <- "LDOF"
invisible(check_design("Two-sided, beta-spending LDOF", p_beta))

cat("\n=== 8. Bound update via update_design_bounds ===\n")
x <- gsDesign(k = 3, test.type = 1, alpha = 0.025, beta = 0.1,
              sfu = sfLDOF, timing = c(0.5, 0.75, 1), n.fix = 100)
cat("original upper:", x$upper$bound, "\n")
y <- update_design_bounds(x, c(0.5, 0.6, 1), spend_all_final = FALSE)
cat("updated upper (spending fn):", y$upper$bound, "\n")
ok("bound update keeps first bound, shifts later bounds",
   abs(y$upper$bound[1] - x$upper$bound[1]) < 1e-6 &&
   abs(y$upper$bound[2] - x$upper$bound[2]) > 0.01)
y2 <- update_design_bounds(x, c(0.5, 0.6, 0.9), spend_all_final = TRUE)
y3 <- update_design_bounds(x, c(0.5, 0.6, 0.9), spend_all_final = FALSE)
cat("spend-all final bound:", y2$upper$bound[3],
    " spending-fn final bound:", y3$upper$bound[3], "\n")
ok("spend-all strategy changes final bound",
   abs(y2$upper$bound[3] - y3$upper$bound[3]) > 0.01)
# Extra interim analysis (k_new = 4) with observed 0.5.
y4 <- update_design_bounds(x, c(0.5, 0.667, 0.833, 1), spend_all_final = FALSE)
ok("update to 4 analyses works", length(y4$upper$bound) == 4)
# Generated update code must parse and run.
uc <- update_code(x, c(0.5, 0.6, 1), FALSE)
env2 <- new.env()
ev2 <- tryCatch({ eval(parse(text = uc), envir = env2); env2$y },
                error = function(e) e)
ok("update code runs and reproduces bounds",
   !inherits(ev2, "error") &&
   all(abs(ev2$upper$bound - y$upper$bound) < 1e-6))
# gsBoundSummary on the updated design.
ut <- tryCatch(as.data.frame(gsBoundSummary(y, ratio = 1)), error = function(e) e)
ok("gsBoundSummary on updated design runs", !inherits(ut, "error"))

cat("\nDone.\n")
