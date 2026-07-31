#!/usr/bin/env Rscript

# Shared utilities for the CLMM release package.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)

  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
  }

  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile, mustWork = FALSE))
  }

  normalizePath(getwd(), mustWork = TRUE)
}

get_release_root <- function() {
  normalizePath(file.path(dirname(get_script_path()), ".."), mustWork = TRUE)
}

resolve_path <- function(path, base_dir) {
  if (is.null(path) || !nzchar(path)) {
    stop("A non-empty path is required.")
  }

  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    return(normalizePath(path, mustWork = FALSE))
  }

  normalizePath(file.path(base_dir, path), mustWork = FALSE)
}

read_env_flag <- function(name, default = FALSE) {
  raw_value <- tolower(Sys.getenv(name, if (default) "true" else "false"))

  if (raw_value %in% c("true", "t", "1", "yes", "y")) {
    return(TRUE)
  }

  if (raw_value %in% c("false", "f", "0", "no", "n")) {
    return(FALSE)
  }

  stop(sprintf("Environment variable %s must be TRUE or FALSE.", name))
}

read_env_int <- function(name, default) {
  raw_value <- Sys.getenv(name, as.character(default))
  value <- suppressWarnings(as.integer(raw_value))

  if (is.na(value)) {
    stop(sprintf("Environment variable %s must be an integer.", name))
  }

  value
}

read_env_num <- function(name, default) {
  raw_value <- Sys.getenv(name, as.character(default))
  value <- suppressWarnings(as.numeric(raw_value))

  if (is.na(value)) {
    stop(sprintf("Environment variable %s must be numeric.", name))
  }

  value
}

detect_default_cores <- function(max_cores = 4L) {
  detected_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))

  if (length(detected_cores) == 0 || is.na(detected_cores) || detected_cores < 1) {
    return(1L)
  }

  as.integer(min(max_cores, detected_cores))
}

load_required_packages <- function(packages) {
  missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing_packages) > 0) {
    stop(
      "Missing required packages: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(packages)
}

source_release_utils <- function() {
  utils_path <- file.path(dirname(get_script_path()), "utils_clmm_release.R")
  source(utils_path, local = FALSE)
}

build_fit_paths <- function(model_tag, default_data_file) {
  release_root <- get_release_root()
  data_path <- resolve_path(
    Sys.getenv("CLMM_DATA_FILE", file.path("data", default_data_file)),
    base_dir = release_root
  )
  results_root <- resolve_path(
    Sys.getenv("CLMM_RESULTS_ROOT", "results"),
    base_dir = release_root
  )

  run_dir <- file.path(results_root, model_tag)
  tab_dir <- file.path(run_dir, "tables")
  diag_dir <- file.path(run_dir, "diagnostics")

  list(
    release_root = release_root,
    data_path = data_path,
    results_root = results_root,
    run_dir = run_dir,
    tab_dir = tab_dir,
    diag_dir = diag_dir
  )
}

build_eval_paths <- function(model_tag) {
  release_root <- get_release_root()
  model_results_root <- resolve_path(
    Sys.getenv("CLMM_MODEL_RESULTS_ROOT", "results"),
    base_dir = release_root
  )
  output_results_root <- resolve_path(
    Sys.getenv("CLMM_OUTPUT_RESULTS_ROOT", Sys.getenv("CLMM_MODEL_RESULTS_ROOT", "results")),
    base_dir = release_root
  )

  input_run_dir <- file.path(model_results_root, model_tag)
  output_run_dir <- file.path(output_results_root, model_tag)
  output_tab_dir <- file.path(output_run_dir, "tables")
  output_diag_dir <- file.path(output_run_dir, "diagnostics")

  list(
    release_root = release_root,
    model_results_root = model_results_root,
    output_results_root = output_results_root,
    input_run_dir = input_run_dir,
    output_run_dir = output_run_dir,
    output_tab_dir = output_tab_dir,
    output_diag_dir = output_diag_dir
  )
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

write_sheet <- function(workbook, sheet_name, data_frame) {
  if (sheet_name %in% names(workbook)) {
    openxlsx::removeWorksheet(workbook, sheet_name)
  }
  openxlsx::addWorksheet(workbook, sheet_name)
  openxlsx::writeData(workbook, sheet_name, data_frame)
}

load_or_create_workbook <- function(path) {
  if (file.exists(path)) {
    openxlsx::loadWorkbook(path)
  } else {
    openxlsx::createWorkbook()
  }
}

save_plot <- function(plot_object, path_base, width = 7, height = 4.5, dpi = 400) {
  ggplot2::ggsave(
    filename = paste0(path_base, ".png"),
    plot = plot_object,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white",
    limitsize = FALSE
  )
  ggplot2::ggsave(
    filename = paste0(path_base, ".pdf"),
    plot = plot_object,
    width = width,
    height = height,
    bg = "white",
    limitsize = FALSE
  )
}

summarise_all_draws <- function(fit) {
  draws <- posterior::as_draws_df(fit)

  five_number_quantiles <- function(x) {
    q <- stats::quantile(x, probs = c(0.025, 0.05, 0.5, 0.95, 0.975), na.rm = TRUE)
    names(q) <- c("q2.5", "q5", "q50", "q95", "q97.5")
    q
  }

  summary_table <- withCallingHandlers(
    posterior::summarise_draws(
      draws,
      "mean",
      "sd",
      five_number_quantiles,
      "rhat",
      "ess_bulk",
      "ess_tail"
    ) |>
      as.data.frame(),
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("The ESS has been capped to avoid unstable estimates", msg, fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )

  names(summary_table) <- sub("^five_number_quantiles\\.", "", names(summary_table))
  summary_table
}

extract_ranef_table <- function(fit, group_name) {
  random_effects <- tryCatch(brms::ranef(fit)[[group_name]][, , "Intercept"], error = function(e) NULL)

  if (is.null(random_effects)) {
    return(tibble::tibble())
  }

  as.data.frame(random_effects) |>
    dplyr::mutate(level = rownames(random_effects)) |>
    dplyr::rename(mu = Estimate, lo = Q2.5, hi = Q97.5)
}

count_fit_rows <- function(fit) {
  fit_rows <- tryCatch(stats::nobs(fit), error = function(e) NA_integer_)

  if (!is.na(fit_rows)) {
    return(as.integer(fit_rows))
  }

  fit_data <- tryCatch(as.data.frame(fit$data), error = function(e) NULL)
  if (is.null(fit_data) || nrow(fit_data) == 0) {
    return(NA_integer_)
  }

  as.integer(nrow(fit_data))
}

summarise_sampler_diagnostics <- function(fit, max_treedepth) {
  sampler_params <- tryCatch(bayesplot::nuts_params(fit), error = function(e) NULL)

  if (is.null(sampler_params) || nrow(sampler_params) == 0) {
    return(
      tibble::tibble(
        n_divergent = NA_integer_,
        n_max_treedepth = NA_integer_,
        max_treedepth_setting = max_treedepth,
        mean_stepsize = NA_real_,
        median_accept_stat = NA_real_
      )
    )
  }

  tibble::tibble(
    n_divergent = sum(sampler_params$Parameter == "divergent__" & sampler_params$Value > 0, na.rm = TRUE),
    n_max_treedepth = sum(
      sampler_params$Parameter == "treedepth__" & sampler_params$Value >= max_treedepth,
      na.rm = TRUE
    ),
    max_treedepth_setting = max_treedepth,
    mean_stepsize = mean(sampler_params$Value[sampler_params$Parameter == "stepsize__"], na.rm = TRUE),
    median_accept_stat = stats::median(
      sampler_params$Value[sampler_params$Parameter == "accept_stat__"],
      na.rm = TRUE
    )
  )
}

safe_loo <- function(fit, moment_match = FALSE) {
  tryCatch(
    brms::loo(fit, pointwise = TRUE, moment_match = moment_match),
    error = function(e) e
  )
}

make_pareto_table <- function(pareto_k) {
  tibble::tibble(
    k_band = c("<=0.5", "0.5-0.7", "0.7-1.0", ">1.0"),
    n = c(
      sum(pareto_k <= 0.5, na.rm = TRUE),
      sum(pareto_k > 0.5 & pareto_k <= 0.7, na.rm = TRUE),
      sum(pareto_k > 0.7 & pareto_k <= 1.0, na.rm = TRUE),
      sum(pareto_k > 1.0, na.rm = TRUE)
    )
  ) |>
    dplyr::mutate(prop = n / sum(n))
}

make_pareto_plot <- function(pareto_table) {
  ggplot2::ggplot(pareto_table, ggplot2::aes(x = k_band, y = n)) +
    ggplot2::geom_col(fill = "#4C78A8", width = 0.7) +
    ggplot2::labs(title = "Pareto k diagnostic bins for LOO", x = "Pareto k range", y = "Count") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

compute_quadratic_weighted_kappa <- function(obs, pred, levels_num) {
  confusion <- table(
    factor(obs, levels = levels_num),
    factor(pred, levels = levels_num)
  )
  observed <- prop.table(confusion)
  row_marginals <- rowSums(observed)
  col_marginals <- colSums(observed)
  expected <- outer(row_marginals, col_marginals)
  n_levels <- length(levels_num)
  weights <- outer(
    seq_len(n_levels),
    seq_len(n_levels),
    function(i, j) ((i - j)^2) / ((n_levels - 1)^2)
  )
  denominator <- sum(weights * expected)

  if (denominator == 0) {
    return(NA_real_)
  }

  1 - sum(weights * observed) / denominator
}

compute_rps <- function(prob_mat, obs_values, levels_num) {
  n_levels <- length(levels_num)
  output <- numeric(length(obs_values))

  for (i in seq_along(obs_values)) {
    cumulative_prob <- cumsum(prob_mat[i, ])[seq_len(n_levels - 1)]
    cumulative_obs <- as.numeric(obs_values[i] <= levels_num[seq_len(n_levels - 1)])
    output[i] <- sum((cumulative_prob - cumulative_obs)^2)
  }

  output
}

coerce_yrep_to_levels <- function(y_rep, level_order) {
  if (is.vector(y_rep)) {
    y_rep <- matrix(y_rep, nrow = 1)
  }

  if (is.numeric(y_rep) || is.integer(y_rep)) {
    y_rep_int <- as.integer(y_rep)
    if (all(y_rep_int %in% seq_along(level_order), na.rm = TRUE)) {
      return(matrix(level_order[y_rep_int], nrow = nrow(y_rep), ncol = ncol(y_rep)))
    }
  }

  output <- apply(y_rep, c(1, 2), as.character)
  if (is.vector(output)) {
    output <- matrix(output, nrow = nrow(y_rep))
  }
  output
}

make_ppc_category_table <- function(y_obs, y_rep, level_order) {
  observed_prop <- prop.table(table(factor(as.character(y_obs), levels = level_order)))
  replicated_prop <- t(apply(y_rep, 1, function(row_i) {
    as.numeric(prop.table(table(factor(as.character(row_i), levels = level_order))))
  }))

  tibble::tibble(
    Y_obs = factor(level_order, levels = level_order, ordered = TRUE),
    obs_prop = as.numeric(observed_prop),
    rep_mean = colMeans(replicated_prop),
    rep_lo95 = apply(replicated_prop, 2, stats::quantile, probs = 0.025, na.rm = TRUE),
    rep_hi95 = apply(replicated_prop, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
  )
}

make_ppc_scalar_summary <- function(y_obs_num, y_rep_num) {
  replicated_mean <- apply(y_rep_num, 1, mean, na.rm = TRUE)
  replicated_sd <- apply(y_rep_num, 1, stats::sd, na.rm = TRUE)

  tibble::tibble(
    stat = c("mean_y_obs", "sd_y_obs"),
    observed = c(mean(y_obs_num, na.rm = TRUE), stats::sd(y_obs_num, na.rm = TRUE)),
    rep_mean = c(mean(replicated_mean, na.rm = TRUE), mean(replicated_sd, na.rm = TRUE)),
    rep_lo95 = c(
      unname(stats::quantile(replicated_mean, 0.025, na.rm = TRUE)),
      unname(stats::quantile(replicated_sd, 0.025, na.rm = TRUE))
    ),
    rep_hi95 = c(
      unname(stats::quantile(replicated_mean, 0.975, na.rm = TRUE)),
      unname(stats::quantile(replicated_sd, 0.975, na.rm = TRUE))
    )
  )
}

make_ppc_category_plot <- function(ppc_category) {
  ggplot2::ggplot(ppc_category, ggplot2::aes(x = Y_obs)) +
    ggplot2::geom_linerange(
      ggplot2::aes(ymin = rep_lo95, ymax = rep_hi95),
      linewidth = 0.6,
      color = "#444444"
    ) +
    ggplot2::geom_point(ggplot2::aes(y = rep_mean), size = 2.1, color = "#444444") +
    ggplot2::geom_point(ggplot2::aes(y = obs_prop), size = 2.2, color = "#C0392B") +
    ggplot2::labs(
      title = "Posterior predictive check: category proportions",
      subtitle = "Black = replicated mean and 95% interval; red = observed",
      x = "Observed category",
      y = "Proportion"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )
}
