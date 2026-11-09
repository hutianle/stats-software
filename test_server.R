# Server-level tests for the Group Sequential Design Explorer (v1).
# Run: Rscript test_server.R
# Uses shiny::testServer() to exercise the server reactives and outputs
# without a browser: default design, endpoint switching, bound update,
# and the save/restore parameter round-trip.
suppressPackageStartupMessages(library(shiny))

ok <- function(name, cond) cat(sprintf("[%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", name))

source("app.R")  # defines ui/server; does not launch the app

# Full default input set, mirroring the UI defaults.
defaults <- list(
  endpoint = "Time-to-event", alpha = 0.025, power = 0.9, ratio = 1,
  time_units = "Months", testing = "one-sided",
  n_interim = 1L, equal_spacing = TRUE, info_frac = "0.5",
  timing_type = "Information-based", calendar_times = "12, 18",
  sfu_name = "LDOF", sfu_gamma = -4, sfu_rho = 1, sfu_nu = 0.5,
  sfu_points_a = "0.4",
  sfl_name = "None", sfl_gamma = -4, sfl_rho = 1, sfl_nu = 0.5,
  sfl_points_a = "0.4",
  adjust_var = "Sample size", hr = 0.7, hr0 = 1, eta = 0,
  hr2 = 0.7, hr02 = 1, eta2 = 0,
  control_spec = "median", control_median = 12,
  control_fail_rate = 0.2, control_fail_time = 12,
  control_surv_rate = 0.5, control_surv_time = 12,
  study_duration = 24, minfup = 12,
  enroll_rate = 10, enroll_duration = 12, rampup = 0,
  p1 = 0.3, p2 = 0.45, bin_scale = "Difference", delta0 = 0,
  effect = 1, sd_c = 1, sd_e = 1,
  n_fixed = 100, delta_std = 0.5,
  use_current_design = TRUE,
  ub_n_analyses = 3L, ub_done = 1L, ub_info_frac = "0.5",
  ub_final_frac = 1, ub_final_strategy = "spending",
  plot_type = "bounds_z"
)

testServer(server, {
  do.call(session$setInputs, defaults)

  cat("=== 1. Default design (TTE, LDOF) ===\n")
  d <- design_obj()
  ok("default design builds", !is.null(d))
  ok("default design power ~ 0.9",
     !is.null(d) && abs(design_power(d) - 0.9) < 0.01)
  txt <- output$design_text
  ok("text output mentions sample size",
     is.character(txt) && grepl("Group-sequential sample size", txt))
  ss <- output$spending_summary
  ok("spending summary shows total alpha",
     is.character(ss) && grepl("Total alpha: 0.025", ss))
  code <- output$design_code
  ok("report code contains gsSurv",
     is.character(code) && grepl("gsSurv", code))
  tb <- tryCatch(output$design_table, error = function(e) e)
  ok("tabular output renders", !inherits(tb, "error") && !is.null(tb))
  pl <- tryCatch(output$design_plot, error = function(e) e)
  ok("bounds plot renders", !inherits(pl, "error") && !is.null(pl))

  cat("\n=== 2. Endpoint switching ===\n")
  session$setInputs(endpoint = "Binomial")
  d_bin <- design_obj()
  ok("binomial design builds",
     !is.null(d_bin) && abs(design_power(d_bin) - 0.9) < 0.01)
  session$setInputs(endpoint = "Normal")
  d_norm <- design_obj()
  ok("normal design builds",
     !is.null(d_norm) && abs(design_power(d_norm) - 0.9) < 0.01)
  session$setInputs(endpoint = "Fixed design N")
  d_fix <- design_obj()
  ok("fixed-design N builds",
     !is.null(d_fix) && abs(design_power(d_fix) - 0.9) < 0.01)

  cat("\n=== 3. Spending function switching ===\n")
  session$setInputs(endpoint = "Time-to-event", sfu_name = "LDPocock")
  d_poc <- design_obj()
  ok("Pocock design builds",
     !is.null(d_poc) && abs(design_power(d_poc) - 0.9) < 0.01)
  session$setInputs(sfu_name = "Points", sfu_points_a = "0.4")
  d_pts <- design_obj()
  ok("custom-points design builds",
     !is.null(d_pts) && abs(design_power(d_pts) - 0.9) < 0.01)

  cat("\n=== 4. Bound update ===\n")
  session$setInputs(sfu_name = "LDOF")
  # Final analysis at 90% of planned information: the spend-all strategy
  # must then differ from the spending-function strategy.
  session$setInputs(ub_n_analyses = 3L, ub_done = 1L,
                    ub_info_frac = "0.5", ub_final_frac = 0.9,
                    ub_final_strategy = "spending")
  u <- updated_bounds()
  ok("bound update returns 3 bounds",
     !is.null(u) && length(u$upper$bound) == 3)
  session$setInputs(ub_final_strategy = "spend-all")
  u2 <- updated_bounds()
  ok("spend-all changes final bound (final at 0.9)",
     !is.null(u) && !is.null(u2) &&
     abs(u2$upper$bound[3] - u$upper$bound[3]) > 0.01)
  # At full information (final fraction 1) both strategies coincide.
  session$setInputs(ub_final_frac = 1, ub_final_strategy = "spending")
  u_full <- updated_bounds()
  session$setInputs(ub_final_strategy = "spend-all")
  u_full2 <- updated_bounds()
  ok("strategies coincide at full information",
     !is.null(u_full) && !is.null(u_full2) &&
     abs(u_full2$upper$bound[3] - u_full$upper$bound[3]) < 1e-9)
  ubt <- tryCatch(ub_table_df(), error = function(e) e)
  ok("updated bounds table renders",
     !inherits(ubt, "error") && !is.null(ubt) && nrow(ubt) > 0)
  ubc <- output$ub_code
  ok("update report code contains n.I",
     is.character(ubc) && grepl("n\\.I", ubc))

  cat("\n=== 5. Save / restore parameter round-trip ===\n")
  p <- collect_params(input)
  ok("collect_params covers all registered inputs",
     all(names(INPUT_TYPES) %in% names(p)) &&
     all(!vapply(p[names(INPUT_TYPES)], is.null, logical(1))))
  # Restoring the collected parameters must reproduce the same design.
  session$setInputs(endpoint = "Normal")
  for (id in names(p)) {
    val <- p[[id]]
    if (is.null(val)) next
    args <- list(val)
    names(args) <- id
    do.call(session$setInputs, args)
  }
  d_restored <- design_obj()
  ok("restored parameters reproduce the TTE design",
     !is.null(d_restored) && inherits(d_restored, "gsSurv") &&
     abs(design_power(d_restored) - 0.9) < 0.01)
})

cat("\nDone.\n")
