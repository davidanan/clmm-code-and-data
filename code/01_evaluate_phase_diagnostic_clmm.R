#!/usr/bin/env Rscript

# Evaluate the phase-diagnostic cumulative link mixed model used in the manuscript.
#
# Default usage:
#   Rscript code/01_evaluate_phase_diagnostic_clmm.R
#
# Optional environment variables:
#   CLMM_MODEL_RESULTS_ROOT   Root directory that contains the saved fit and prep objects.
#   CLMM_OUTPUT_RESULTS_ROOT  Root directory for evaluation outputs.
#   CLMM_DO_LAYER1            TRUE/FALSE, run convergence diagnostics.
#   CLMM_DO_PPC               TRUE/FALSE, run posterior predictive checks.
#   CLMM_DO_LAYER3            TRUE/FALSE, run predictive summaries.
#   CLMM_DO_LOO               TRUE/FALSE, run leave-one-out diagnostics.
#   CLMM_EVAL_NDRAWS          Posterior draws used for predictive summaries.
#   CLMM_PPC_NDRAWS_PLOT      Posterior draws used in pp_check plots.
#   CLMM_PPC_NDRAWS_TABLE     Posterior draws used in PPC summary tables.

script_path <- {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)

  if (length(file_arg) > 0) {
    normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)
  } else if (!is.null(sys.frames()[[1]]$ofile)) {
    normalizePath(sys.frames()[[1]]$ofile, mustWork = FALSE)
  } else {
    normalizePath(file.path(getwd(), "code", "01_evaluate_phase_diagnostic_clmm.R"), mustWork = FALSE)
  }
}

source(file.path(dirname(script_path), "utils_clmm_release.R"))

load_required_packages(c(
  "dplyr", "readr", "tibble", "tidyr",
  "brms", "posterior", "bayesplot", "ggplot2", "openxlsx"
))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(brms)
  library(posterior)
  library(bayesplot)
  library(ggplot2)
  library(openxlsx)
})

ensure_eval_df <- function(prep_object) {
  if (is.null(prep_object$df) || nrow(prep_object$df) == 0) {
    stop("prep$df is empty.")
  }

  required_columns <- c(
    "source_row_id", "SiteID", "EventID", "SeqID", "SeqSegID",
    "t", "Y_obs_num", "Y_sim_num", "Y_obs_ord", "Phase"
  )

  missing_columns <- setdiff(required_columns, names(prep_object$df))
  if (length(missing_columns) > 0) {
    stop("prep$df is missing columns: ", paste(missing_columns, collapse = ", "))
  }

  prep_object$df
}

compute_in_sample_ordinal_metrics <- function(fit, eval_data, ndraws = 120) {
  prob_array <- posterior_epred(
    fit,
    newdata = eval_data,
    ndraws = ndraws,
    re_formula = NULL,
    allow_new_levels = FALSE,
    sort = FALSE
  )

  if (length(dim(prob_array)) != 3) {
    stop("posterior_epred returned unexpected dimensions.")
  }
  if (dim(prob_array)[2] != nrow(eval_data)) {
    stop("posterior_epred row count does not match eval_data.")
  }

  prob_mean <- apply(prob_array, c(2, 3), mean)
  if (is.vector(prob_mean)) {
    prob_mean <- matrix(prob_mean, nrow = nrow(eval_data), byrow = TRUE)
  }

  y_levels <- 0:4
  obs <- eval_data$Y_obs_num
  pred_class <- y_levels[max.col(prob_mean, ties.method = "first")]
  max_prob <- apply(prob_mean, 1, max)
  expected_y <- as.vector(prob_mean %*% y_levels)
  abs_error <- abs(expected_y - obs)
  squared_error <- (expected_y - obs)^2
  rps <- compute_rps(prob_mean, obs_values = obs, levels_num = y_levels)

  category_calibration <- tibble(
    Y_obs = factor(y_levels, levels = y_levels, ordered = TRUE),
    obs_n = as.numeric(table(factor(obs, levels = y_levels))),
    obs_prop = as.numeric(prop.table(table(factor(obs, levels = y_levels)))),
    pred_prop = colMeans(prob_mean)
  )

  confusion_matrix <- as.data.frame.matrix(
    table(
      obs = factor(obs, levels = y_levels),
      pred = factor(pred_class, levels = y_levels)
    )
  ) |>
    tibble::rownames_to_column("obs")

  row_predictions <- tibble(
    row_id = seq_len(nrow(eval_data)),
    source_row_id = eval_data$source_row_id,
    SiteID = eval_data$SiteID,
    EventID = eval_data$EventID,
    SeqSegID = eval_data$SeqSegID,
    t = eval_data$t,
    Phase = eval_data$Phase,
    Y_sim = eval_data$Y_sim_num,
    Y_obs = obs,
    Y_pred_class = pred_class,
    Y_pred_E = expected_y,
    pred_max_prob = max_prob,
    abs_error_E = abs_error,
    rps = rps
  )

  metrics <- tibble(
    n_obs = nrow(eval_data),
    ndraws_eval = ndraws,
    exact_accuracy = mean(pred_class == obs),
    adjacent_accuracy = mean(abs(pred_class - obs) <= 1),
    quadratic_weighted_kappa = compute_quadratic_weighted_kappa(obs, pred_class, y_levels),
    mae_expected_y = mean(abs_error),
    rmse_expected_y = sqrt(mean(squared_error)),
    mean_rps = mean(rps)
  )

  list(
    metrics = metrics,
    category_calibration = category_calibration,
    confusion = confusion_matrix,
    row_predictions = row_predictions
  )
}

compute_raw_stage_summaries <- function(eval_data, phase_levels) {
  summary_data <- eval_data |>
    mutate(
      under = Y_obs_num > Y_sim_num,
      match = Y_obs_num == Y_sim_num,
      over = Y_obs_num < Y_sim_num,
      Delta_num = Y_obs_num - Y_sim_num
    )

  overall <- summary_data |>
    group_by(Phase) |>
    summarise(
      n = n(),
      P_under = mean(under),
      P_match = mean(match),
      P_over = mean(over),
      E_delta = mean(Delta_num),
      .groups = "drop"
    ) |>
    mutate(Phase = factor(Phase, levels = phase_levels)) |>
    arrange(Phase)

  by_event <- summary_data |>
    group_by(EventID, Phase) |>
    summarise(
      n = n(),
      P_under = mean(under),
      P_match = mean(match),
      P_over = mean(over),
      E_delta = mean(Delta_num),
      .groups = "drop"
    ) |>
    mutate(Phase = factor(Phase, levels = phase_levels)) |>
    arrange(EventID, Phase)

  list(overall = overall, by_event = by_event)
}

summarise_draw_metric <- function(draw_summary, phase_levels) {
  draw_summary |>
    pivot_longer(
      cols = c("P_under", "P_match", "P_over", "E_delta"),
      names_to = "metric",
      values_to = "value"
    ) |>
    group_by(Phase, metric) |>
    summarise(
      mean = mean(value),
      q2.5 = quantile(value, 0.025),
      q50 = quantile(value, 0.5),
      q97.5 = quantile(value, 0.975),
      .groups = "drop"
    ) |>
    mutate(Phase = factor(Phase, levels = phase_levels)) |>
    arrange(metric, Phase)
}

make_stage_contrasts <- function(draw_summary) {
  phase_pairs <- list(
    c("peak", "rising"),
    c("peak", "recession"),
    c("rising", "recession")
  )

  long_summary <- draw_summary |>
    pivot_longer(
      cols = c("P_under", "P_match", "P_over", "E_delta"),
      names_to = "metric",
      values_to = "value"
    )

  bind_rows(lapply(unique(long_summary$metric), function(metric_name) {
    wide_summary <- long_summary |>
      filter(metric == metric_name) |>
      select(draw, Phase, value) |>
      pivot_wider(names_from = Phase, values_from = value)

    bind_rows(lapply(phase_pairs, function(pair) {
      diff_value <- wide_summary[[pair[[1]]]] - wide_summary[[pair[[2]]]]
      tibble(
        metric = metric_name,
        contrast = paste(pair[[1]], "-", pair[[2]]),
        mean = mean(diff_value),
        q2.5 = quantile(diff_value, 0.025),
        q50 = quantile(diff_value, 0.5),
        q97.5 = quantile(diff_value, 0.975),
        p_gt_0 = mean(diff_value > 0),
        p_lt_0 = mean(diff_value < 0)
      )
    }))
  }))
}

compute_stage_model_summaries <- function(fit, eval_data, phase_levels, ndraws = 120) {
  y_levels <- 0:4
  n_rows <- nrow(eval_data)
  n_levels <- length(y_levels)
  event_ids <- levels(eval_data$EventID)
  event_indices <- split(seq_len(n_rows), eval_data$EventID)

  category_matrix <- matrix(rep(y_levels, each = n_rows), nrow = n_rows, ncol = n_levels)
  ysim_matrix <- matrix(eval_data$Y_sim_num, nrow = n_rows, ncol = n_levels)
  under_mask <- category_matrix > ysim_matrix
  match_mask <- category_matrix == ysim_matrix
  over_mask <- category_matrix < ysim_matrix
  delta_mask <- category_matrix - ysim_matrix

  draw_tables <- list()
  event_tables <- list()

  for (phase_name in phase_levels) {
    new_data <- eval_data
    new_data$Phase <- factor(phase_name, levels = phase_levels)

    prob_array <- posterior_epred(
      fit,
      newdata = new_data,
      ndraws = ndraws,
      re_formula = NULL,
      allow_new_levels = FALSE,
      sort = FALSE
    )

    if (length(dim(prob_array)) != 3) {
      stop("posterior_epred returned unexpected dimensions in stage summary.")
    }
    if (dim(prob_array)[2] != n_rows) {
      stop("posterior_epred row count mismatch in stage summary.")
    }

    prob_mean <- apply(prob_array, c(2, 3), mean)
    if (is.vector(prob_mean)) {
      prob_mean <- matrix(prob_mean, nrow = n_rows, byrow = TRUE)
    }

    event_tables[[phase_name]] <- tibble(
      EventID = eval_data$EventID,
      Phase = factor(phase_name, levels = phase_levels),
      P_under = rowSums(prob_mean * under_mask),
      P_match = rowSums(prob_mean * match_mask),
      P_over = rowSums(prob_mean * over_mask),
      E_delta = rowSums(prob_mean * delta_mask)
    ) |>
      group_by(EventID, Phase) |>
      summarise(
        P_under = mean(P_under),
        P_match = mean(P_match),
        P_over = mean(P_over),
        E_delta = mean(E_delta),
        .groups = "drop"
      )

    draw_metrics <- matrix(NA_real_, nrow = dim(prob_array)[1], ncol = 4)
    colnames(draw_metrics) <- c("P_under", "P_match", "P_over", "E_delta")

    for (draw_index in seq_len(dim(prob_array)[1])) {
      prob_draw <- prob_array[draw_index, , ]
      if (is.null(dim(prob_draw))) {
        prob_draw <- matrix(prob_draw, nrow = n_rows, byrow = TRUE)
      }

      under_row <- rowSums(prob_draw * under_mask)
      match_row <- rowSums(prob_draw * match_mask)
      over_row <- rowSums(prob_draw * over_mask)
      delta_row <- rowSums(prob_draw * delta_mask)

      event_means <- do.call(
        rbind,
        lapply(event_ids, function(event_id) {
          idx <- event_indices[[event_id]]
          c(
            mean(under_row[idx]),
            mean(match_row[idx]),
            mean(over_row[idx]),
            mean(delta_row[idx])
          )
        })
      )

      draw_metrics[draw_index, ] <- colMeans(event_means)
    }

    draw_tables[[phase_name]] <- tibble(
      draw = seq_len(nrow(draw_metrics)),
      Phase = factor(phase_name, levels = phase_levels),
      P_under = draw_metrics[, "P_under"],
      P_match = draw_metrics[, "P_match"],
      P_over = draw_metrics[, "P_over"],
      E_delta = draw_metrics[, "E_delta"]
    )
  }

  draw_summary <- bind_rows(draw_tables)
  event_summary <- bind_rows(event_tables) |>
    arrange(EventID, Phase)

  list(
    overall = summarise_draw_metric(draw_summary, phase_levels = phase_levels),
    by_event = event_summary,
    contrasts = make_stage_contrasts(draw_summary)
  )
}

model_tag <- "clmm_phase_0410"
paths <- build_eval_paths(model_tag)

do_layer1 <- read_env_flag("CLMM_DO_LAYER1", default = TRUE)
do_ppc <- read_env_flag("CLMM_DO_PPC", default = TRUE)
do_layer3 <- read_env_flag("CLMM_DO_LAYER3", default = TRUE)
do_loo <- read_env_flag("CLMM_DO_LOO", default = FALSE)
reuse_existing_loo <- read_env_flag("CLMM_REUSE_EXISTING_LOO", default = TRUE)
moment_match_loo <- read_env_flag("CLMM_MOMENT_MATCH_LOO", default = FALSE)
eval_ndraws <- read_env_int("CLMM_EVAL_NDRAWS", 120L)
pp_ndraws_plot <- read_env_int("CLMM_PPC_NDRAWS_PLOT", 80L)
pp_ndraws_table <- read_env_int("CLMM_PPC_NDRAWS_TABLE", 200L)
fig_width <- read_env_num("CLMM_FIG_WIDTH", 7.5)
fig_height <- read_env_num("CLMM_FIG_HEIGHT", 5.2)
fig_dpi <- read_env_int("CLMM_FIG_DPI", 300L)
pp_check_type <- Sys.getenv("CLMM_PPC_TYPE", "bars")

ensure_dir(paths$output_results_root)
ensure_dir(paths$output_run_dir)
ensure_dir(paths$output_tab_dir)
ensure_dir(paths$output_diag_dir)

fit_path <- file.path(paths$input_run_dir, paste0(model_tag, "_fit.rds"))
prep_path <- file.path(paths$input_run_dir, paste0(model_tag, "_prep_objects.rds"))
xlsx_path <- file.path(paths$output_tab_dir, paste0(model_tag, "_results.xlsx"))
loo_path <- file.path(paths$output_diag_dir, paste0(model_tag, "_loo.rds"))
eval_runlog_path <- file.path(paths$output_diag_dir, paste0(model_tag, "_eval_runlog.txt"))
eval_session_path <- file.path(paths$output_diag_dir, paste0(model_tag, "_eval_sessionInfo.txt"))
all_parameters_path <- file.path(paths$output_diag_dir, paste0(model_tag, "_summary_all_parameters.csv"))
convergence_path <- file.path(paths$output_diag_dir, paste0(model_tag, "_convergence_summary.csv"))
sampler_path <- file.path(paths$output_diag_dir, paste0(model_tag, "_sampler_diagnostics.csv"))
ppc_category_path <- file.path(paths$output_diag_dir, paste0(model_tag, "_ppc_category_proportions.csv"))
ppc_scalar_path <- file.path(paths$output_diag_dir, paste0(model_tag, "_ppc_scalar_summary.csv"))
ppc_plot_base <- file.path(paths$output_diag_dir, paste0(model_tag, "_ppc_bars"))
ppc_category_plot_base <- file.path(paths$output_diag_dir, paste0(model_tag, "_ppc_category_proportions"))
in_sample_path <- file.path(paths$output_diag_dir, paste0(model_tag, "_in_sample_metrics.csv"))
calibration_path <- file.path(paths$output_diag_dir, paste0(model_tag, "_category_calibration.csv"))
confusion_path <- file.path(paths$output_diag_dir, paste0(model_tag, "_confusion_matrix.csv"))
row_predictions_path <- file.path(paths$output_diag_dir, paste0(model_tag, "_row_predictions.csv"))
raw_stage_overall_path <- file.path(paths$output_diag_dir, paste0(model_tag, "_raw_stage_overall.csv"))
raw_stage_event_path <- file.path(paths$output_diag_dir, paste0(model_tag, "_raw_stage_by_event.csv"))
model_stage_overall_path <- file.path(paths$output_diag_dir, paste0(model_tag, "_model_stage_overall.csv"))
model_stage_event_path <- file.path(paths$output_diag_dir, paste0(model_tag, "_model_stage_by_event.csv"))
stage_contrast_path <- file.path(paths$output_diag_dir, paste0(model_tag, "_stage_contrasts.csv"))
loo_summary_path <- file.path(paths$output_diag_dir, paste0(model_tag, "_loo_summary.csv"))
pareto_table_path <- file.path(paths$output_diag_dir, paste0(model_tag, "_pareto_k_table.csv"))
pareto_plot_base <- file.path(paths$output_diag_dir, paste0(model_tag, "_pareto_k"))

if (!file.exists(fit_path)) {
  stop("Cannot find fit file: ", fit_path)
}
if (!file.exists(prep_path)) {
  stop("Cannot find prep file: ", prep_path)
}

fit <- readRDS(fit_path)
prep <- readRDS(prep_path)
eval_data <- ensure_eval_df(prep)

fit_n <- count_fit_rows(fit)
if (!is.na(fit_n) && fit_n != nrow(eval_data)) {
  stop("Model row count does not match the prep rows.")
}

phase_levels <- prep$config$phase_levels
max_treedepth <- fit$call$control$max_treedepth %||% NA_integer_
workbook <- load_or_create_workbook(xlsx_path)

all_parameters <- NULL
convergence_summary <- NULL
sampler_summary <- NULL

if (do_layer1) {
  message("Running convergence diagnostics ...")

  all_parameters <- summarise_all_draws(fit) |>
    rename(param = variable) |>
    arrange(desc(rhat))
  readr::write_csv(all_parameters, all_parameters_path)
  write_sheet(workbook, "all_parameters", all_parameters)

  convergence_summary <- tibble(
    max_rhat = max(all_parameters$rhat, na.rm = TRUE),
    min_ess_bulk = min(all_parameters$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(all_parameters$ess_tail, na.rm = TRUE),
    n_rhat_gt_1_01 = sum(all_parameters$rhat > 1.01, na.rm = TRUE),
    n_rhat_gt_1_05 = sum(all_parameters$rhat > 1.05, na.rm = TRUE)
  )
  sampler_summary <- summarise_sampler_diagnostics(fit, max_treedepth = max_treedepth)

  readr::write_csv(convergence_summary, convergence_path)
  readr::write_csv(sampler_summary, sampler_path)
  write_sheet(workbook, "convergence_summary", convergence_summary)
  write_sheet(workbook, "sampler_diagnostics", sampler_summary)
}

if (do_ppc) {
  message("Running posterior predictive checks ...")

  ppc_plot <- pp_check(fit, type = pp_check_type, ndraws = pp_ndraws_plot) +
    labs(
      title = "Posterior predictive check",
      subtitle = paste0("Type = ", pp_check_type, ", ndraws = ", pp_ndraws_plot)
    ) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
  save_plot(ppc_plot, ppc_plot_base, width = fig_width, height = fig_height, dpi = fig_dpi)

  y_rep <- posterior_predict(
    fit,
    newdata = eval_data,
    ndraws = pp_ndraws_table,
    re_formula = NULL,
    allow_new_levels = FALSE,
    sort = FALSE
  )
  y_rep_chr <- coerce_yrep_to_levels(y_rep, level_order = as.character(0:4))

  ppc_category <- make_ppc_category_table(
    y_obs = eval_data$Y_obs_ord,
    y_rep = y_rep_chr,
    level_order = as.character(0:4)
  )
  readr::write_csv(ppc_category, ppc_category_path)
  write_sheet(workbook, "ppc_category_props", ppc_category)

  y_rep_num <- matrix(as.numeric(y_rep_chr), nrow = nrow(y_rep_chr), ncol = ncol(y_rep_chr))
  ppc_scalar <- make_ppc_scalar_summary(y_obs_num = eval_data$Y_obs_num, y_rep_num = y_rep_num)
  readr::write_csv(ppc_scalar, ppc_scalar_path)
  write_sheet(workbook, "ppc_scalar_summary", ppc_scalar)

  ppc_category_plot <- make_ppc_category_plot(ppc_category)
  save_plot(
    ppc_category_plot,
    ppc_category_plot_base,
    width = fig_width,
    height = fig_height,
    dpi = fig_dpi
  )
}

if (do_layer3) {
  message("Running predictive summaries ...")

  prediction_summary <- compute_in_sample_ordinal_metrics(fit, eval_data = eval_data, ndraws = eval_ndraws)
  readr::write_csv(prediction_summary$metrics, in_sample_path)
  readr::write_csv(prediction_summary$category_calibration, calibration_path)
  readr::write_csv(prediction_summary$confusion, confusion_path)
  readr::write_csv(prediction_summary$row_predictions, row_predictions_path)
  write_sheet(workbook, "in_sample_metrics", prediction_summary$metrics)
  write_sheet(workbook, "category_calibration", prediction_summary$category_calibration)
  write_sheet(workbook, "confusion_matrix", prediction_summary$confusion)

  raw_stage <- compute_raw_stage_summaries(eval_data, phase_levels = phase_levels)
  readr::write_csv(raw_stage$overall, raw_stage_overall_path)
  readr::write_csv(raw_stage$by_event, raw_stage_event_path)
  write_sheet(workbook, "raw_stage_overall", raw_stage$overall)
  write_sheet(workbook, "raw_stage_by_event", raw_stage$by_event)

  model_stage <- compute_stage_model_summaries(fit, eval_data, phase_levels = phase_levels, ndraws = eval_ndraws)
  readr::write_csv(model_stage$overall, model_stage_overall_path)
  readr::write_csv(model_stage$by_event, model_stage_event_path)
  readr::write_csv(model_stage$contrasts, stage_contrast_path)
  write_sheet(workbook, "model_stage_overall", model_stage$overall)
  write_sheet(workbook, "model_stage_by_event", model_stage$by_event)
  write_sheet(workbook, "stage_contrasts", model_stage$contrasts)

  if (do_loo) {
    loo_result <- if (reuse_existing_loo && file.exists(loo_path)) {
      readRDS(loo_path)
    } else {
      safe_loo(fit, moment_match = moment_match_loo)
    }

    if (inherits(loo_result, "error")) {
      stop("LOO failed: ", conditionMessage(loo_result))
    }

    saveRDS(loo_result, loo_path)

    pareto_k <- loo_result$diagnostics$pareto_k
    loo_summary <- tibble(
      status = "ok",
      elpd_loo = loo_result$estimates["elpd_loo", "Estimate"],
      se_elpd_loo = loo_result$estimates["elpd_loo", "SE"],
      looic = loo_result$estimates["looic", "Estimate"],
      se_looic = loo_result$estimates["looic", "SE"],
      p_loo = loo_result$estimates["p_loo", "Estimate"],
      n_pareto_k_gt_0_7 = sum(pareto_k > 0.7, na.rm = TRUE),
      n_pareto_k_gt_1 = sum(pareto_k > 1, na.rm = TRUE)
    )
    pareto_table <- make_pareto_table(pareto_k)

    readr::write_csv(loo_summary, loo_summary_path)
    readr::write_csv(pareto_table, pareto_table_path)
    write_sheet(workbook, "loo_summary", loo_summary)
    write_sheet(workbook, "pareto_k_table", pareto_table)

    pareto_plot <- make_pareto_plot(pareto_table)
    save_plot(pareto_plot, pareto_plot_base, width = 6.4, height = 4.5, dpi = fig_dpi)
  }
}

saveWorkbook(workbook, xlsx_path, overwrite = TRUE)

writeLines(
  c(
    paste("Model tag:", model_tag),
    paste("Fit path:", fit_path),
    paste("Prep path:", prep_path),
    paste("Output diagnostics dir:", paths$output_diag_dir),
    paste("Layer 1:", do_layer1),
    paste("PPC:", do_ppc),
    paste("Layer 3:", do_layer3),
    paste("LOO enabled:", do_loo),
    if (!is.null(convergence_summary)) paste("Max Rhat:", convergence_summary$max_rhat) else NULL,
    if (!is.null(sampler_summary)) paste("N divergent:", sampler_summary$n_divergent) else NULL
  ),
  con = eval_runlog_path
)
capture.output(sessionInfo(), file = eval_session_path)

message("Done. Evaluation outputs saved under: ", paths$output_diag_dir)
