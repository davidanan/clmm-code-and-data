#!/usr/bin/env Rscript

# Fit the phase-diagnostic cumulative link mixed model used in the manuscript.
#
# Default usage:
#   Rscript code/01_fit_phase_diagnostic_clmm.R
#
# Optional environment variables:
#   CLMM_DATA_FILE         Relative or absolute path to the model-input CSV.
#   CLMM_RESULTS_ROOT      Relative or absolute path to the output root directory.
#   CLMM_DO_FIT            TRUE/FALSE, fit the model or load an existing fit.
#   CLMM_PREP_ONLY         TRUE/FALSE, stop after preprocessing.
#   CLMM_CHAINS            Number of MCMC chains.
#   CLMM_ITER              Total iterations per chain.
#   CLMM_WARMUP            Warmup iterations per chain.
#   CLMM_CORES             Number of CPU cores to use.
#   CLMM_ADAPT_DELTA       Stan adapt_delta setting.
#   CLMM_MAX_TREEDEPTH     Stan max_treedepth setting.
#   CLMM_REFRESH           brms refresh interval.

script_path <- {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)

  if (length(file_arg) > 0) {
    normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)
  } else if (!is.null(sys.frames()[[1]]$ofile)) {
    normalizePath(sys.frames()[[1]]$ofile, mustWork = FALSE)
  } else {
    normalizePath(file.path(getwd(), "code", "01_fit_phase_diagnostic_clmm.R"), mustWork = FALSE)
  }
}

source(file.path(dirname(script_path), "utils_clmm_release.R"))

load_required_packages(c(
  "dplyr", "purrr", "readr", "stringr", "tibble",
  "brms", "cmdstanr", "posterior", "openxlsx"
))

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(brms)
  library(cmdstanr)
  library(posterior)
  library(openxlsx)
})

assert_cmdstan_ready <- function() {
  cmdstan_version <- tryCatch(cmdstanr::cmdstan_version(error_on_NA = FALSE), error = function(e) NA)
  if (is.na(cmdstan_version)) {
    stop("cmdstanr is not configured. Please install and configure CmdStan before fitting.")
  }
}

make_missingness_report <- function(data_frame, required_columns) {
  tibble(
    variable = required_columns,
    n_missing = map_int(required_columns, ~ sum(is.na(data_frame[[.x]]))),
    prop_missing = n_missing / nrow(data_frame)
  ) |>
    arrange(desc(n_missing), variable)
}

make_prefit_validation <- function(model_data, sequence_overview, required_columns, phase_levels, response_levels) {
  tibble(
    check = c(
      "required_columns_present",
      "n_rows_after_complete_case",
      "y_obs_in_0_4",
      "y_sim_in_0_4",
      "y_obs_has_all_5_levels",
      "y_sim_has_all_5_levels",
      "phase_levels_valid",
      "t_positive_integer",
      "no_duplicate_seqseg_t",
      "n_segments_len_1"
    ),
    status = c(
      TRUE,
      nrow(model_data) > 0,
      all(model_data$Y_obs_num %in% 0:4),
      all(model_data$Y_sim_num %in% 0:4),
      identical(levels(model_data$Y_obs_ord), response_levels),
      setequal(unique(model_data$Y_sim_num), 0:4),
      !any(is.na(model_data$Phase)) && identical(levels(model_data$Phase), phase_levels),
      all(!is.na(model_data$t)) &&
        all(model_data$t > 0) &&
        all(abs(model_data$t - round(model_data$t)) < .Machine$double.eps^0.5),
      nrow(model_data |> count(SeqSegID, t, name = "n") |> filter(n > 1)) == 0,
      TRUE
    ),
    value = c(
      paste(required_columns, collapse = "; "),
      as.character(nrow(model_data)),
      as.character(all(model_data$Y_obs_num %in% 0:4)),
      as.character(all(model_data$Y_sim_num %in% 0:4)),
      paste(levels(model_data$Y_obs_ord), collapse = ","),
      paste(sort(unique(model_data$Y_sim_num)), collapse = ","),
      paste(levels(model_data$Phase), collapse = ","),
      as.character(all(!is.na(model_data$t)) && all(model_data$t > 0)),
      as.character(nrow(model_data |> count(SeqSegID, t, name = "n") |> filter(n > 1)) == 0),
      as.character(sum(sequence_overview$n_time == 1))
    ),
    note = c(
      "Required columns are available in the source table.",
      "Number of rows after complete-case filtering.",
      "Observed levels remain within the 0-4 ordinal scale.",
      "Simulated levels remain within the 0-4 ordinal scale.",
      "Observed response levels span all five ordinal categories.",
      "Simulated anchor levels span all five ordinal categories.",
      "Flood_Phase contains only rising, peak, and recession.",
      "Within-segment time index t is a positive integer.",
      "Each (SeqSegID, t) pair is unique.",
      "Length-1 segments are permitted but not expected."
    )
  )
}

model_tag <- "clmm_phase_0410"
phase_levels <- c("rising", "peak", "recession")
response_levels <- as.character(0:4)
event_reference <- "1"
required_columns <- c(
  "source_row_id", "Event_ID", "Site_ID", "SeqID", "SeqSegID", "Actual_Time",
  "t", "Y_obs", "Y_sim", "Flood_Phase"
)

paths <- build_fit_paths(
  model_tag = model_tag,
  default_data_file = "clmm_model_input_v0410.csv"
)

do_fit <- read_env_flag("CLMM_DO_FIT", default = TRUE)
prep_only <- read_env_flag("CLMM_PREP_ONLY", default = FALSE)
seed <- read_env_int("CLMM_SEED", 20260410L)
chains <- read_env_int("CLMM_CHAINS", 4L)
iter <- read_env_int("CLMM_ITER", 3000L)
warmup <- read_env_int("CLMM_WARMUP", 1500L)
cores <- read_env_int("CLMM_CORES", detect_default_cores())
adapt_delta <- read_env_num("CLMM_ADAPT_DELTA", 0.98)
max_treedepth <- read_env_int("CLMM_MAX_TREEDEPTH", 12L)
refresh <- read_env_int("CLMM_REFRESH", 40L)

control <- list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)

ensure_dir(paths$results_root)
ensure_dir(paths$run_dir)
ensure_dir(paths$tab_dir)
ensure_dir(paths$diag_dir)

fit_path <- file.path(paths$run_dir, paste0(model_tag, "_fit.rds"))
prep_path <- file.path(paths$run_dir, paste0(model_tag, "_prep_objects.rds"))
xlsx_path <- file.path(paths$tab_dir, paste0(model_tag, "_results.xlsx"))
model_info_path <- file.path(paths$tab_dir, paste0(model_tag, "_model_info.csv"))
missingness_path <- file.path(paths$tab_dir, paste0(model_tag, "_missingness_report.csv"))
used_rows_path <- file.path(paths$tab_dir, paste0(model_tag, "_used_rows.csv"))
prefit_path <- file.path(paths$tab_dir, paste0(model_tag, "_prefit_validation.csv"))
all_parameters_path <- file.path(paths$tab_dir, paste0(model_tag, "_summary_all_parameters.csv"))
fixed_effects_path <- file.path(paths$tab_dir, paste0(model_tag, "_fixed_effects.csv"))
thresholds_path <- file.path(paths$tab_dir, paste0(model_tag, "_thresholds.csv"))
hyperparameters_path <- file.path(paths$tab_dir, paste0(model_tag, "_hyperparameters.csv"))
fixed_or_path <- file.path(paths$tab_dir, paste0(model_tag, "_fixed_effects_with_OR.csv"))
site_re_path <- file.path(paths$tab_dir, paste0(model_tag, "_site_random_intercepts.csv"))
sequence_re_path <- file.path(paths$tab_dir, paste0(model_tag, "_sequence_random_intercepts.csv"))
data_overview_path <- file.path(paths$tab_dir, paste0(model_tag, "_data_overview.csv"))
sequence_overview_path <- file.path(paths$tab_dir, paste0(model_tag, "_sequence_overview.csv"))
fit_runlog_path <- file.path(paths$diag_dir, paste0(model_tag, "_fit_runlog.txt"))
fit_session_path <- file.path(paths$diag_dir, paste0(model_tag, "_fit_sessionInfo.txt"))

assert_cmdstan_ready()
options(mc.cores = cores)

raw_data <- readr::read_csv(paths$data_path, show_col_types = FALSE)

missingness_report <- make_missingness_report(raw_data, required_columns)
readr::write_csv(missingness_report, missingness_path)

model_data <- raw_data |>
  filter(complete.cases(across(all_of(required_columns)))) |>
  mutate(
    SiteID = factor(Site_ID),
    EventID = stats::relevel(factor(Event_ID), ref = event_reference),
    SeqID = factor(SeqID),
    SeqSegID = factor(SeqSegID),
    Y_obs_num = as.integer(Y_obs),
    Y_sim_num = as.integer(Y_sim),
    Y_obs_ord = factor(Y_obs_num, levels = 0:4, ordered = TRUE),
    Y_sim_mo = Y_sim_num + 1L,
    Phase = factor(Flood_Phase, levels = phase_levels),
    Delta_num = Y_obs_num - Y_sim_num
  )

if (nrow(model_data) == 0) {
  stop("No complete-case rows remain after filtering.")
}

if (anyNA(model_data$Phase)) {
  stop("Flood_Phase contains values outside the expected phase levels.")
}

dropped_rows <- nrow(raw_data) - nrow(model_data)

sequence_overview <- model_data |>
  count(SeqSegID, SeqID, SiteID, EventID, name = "n_time") |>
  arrange(desc(n_time))

data_overview <- tibble(
  n_obs_raw = nrow(raw_data),
  n_obs_used = nrow(model_data),
  n_rows_dropped_missing = dropped_rows,
  n_sites = n_distinct(model_data$SiteID),
  n_events = n_distinct(model_data$EventID),
  n_base_sequences = n_distinct(model_data$SeqID),
  n_segments = n_distinct(model_data$SeqSegID),
  n_segments_len_1 = sum(sequence_overview$n_time == 1),
  min_segment_len = min(sequence_overview$n_time),
  median_segment_len = median(sequence_overview$n_time),
  mean_segment_len = mean(sequence_overview$n_time),
  max_segment_len = max(sequence_overview$n_time)
)

prefit_validation <- make_prefit_validation(
  model_data = model_data,
  sequence_overview = sequence_overview,
  required_columns = required_columns,
  phase_levels = phase_levels,
  response_levels = response_levels
)

readr::write_csv(prefit_validation, prefit_path)
readr::write_csv(model_data, used_rows_path)

prep_object <- list(
  df = model_data,
  scalers = list(),
  config = list(
    data_path = paths$data_path,
    results_root = paths$results_root,
    run_dir = paths$run_dir,
    tab_dir = paths$tab_dir,
    diag_dir = paths$diag_dir,
    model_tag = model_tag,
    phase_levels = phase_levels,
    response_levels = response_levels,
    formula = "Y_obs_ord ~ mo(Y_sim_mo) + Phase + EventID + (1 | SiteID) + (1 | SeqSegID)"
  ),
  data_overview = data_overview,
  seq_overview = sequence_overview,
  prefit_validation = prefit_validation
)
saveRDS(prep_object, prep_path)

if (prep_only) {
  message("Preprocessing completed. Prep objects were saved and model fitting was skipped.")
  quit(save = "no", status = 0)
}

formula_text <- "Y_obs_ord ~ mo(Y_sim_mo) + Phase + EventID + (1 | SiteID) + (1 | SeqSegID)"
model_formula <- bf(Y_obs_ord ~ mo(Y_sim_mo) + Phase + EventID + (1 | SiteID) + (1 | SeqSegID))

model_info <- tibble(
  item = c(
    "model_tag", "formula", "family", "seed", "chains", "iter", "warmup",
    "cores", "adapt_delta", "max_treedepth", "do_fit", "prep_only"
  ),
  value = c(
    model_tag,
    formula_text,
    "cumulative(link = 'logit')",
    as.character(seed),
    as.character(chains),
    as.character(iter),
    as.character(warmup),
    as.character(cores),
    as.character(adapt_delta),
    as.character(max_treedepth),
    as.character(do_fit),
    as.character(prep_only)
  )
)
readr::write_csv(model_info, model_info_path)

if (do_fit) {
  priors <- c(
    prior(normal(0, 1), class = "b"),
    prior(student_t(3, 0, 2.5), class = "Intercept"),
    prior(student_t(3, 0, 1), class = "sd")
  )

  message("Fitting formula: ", formula_text)

  fit <- brm(
    formula = model_formula,
    data = model_data,
    family = cumulative(link = "logit"),
    prior = priors,
    chains = chains,
    iter = iter,
    warmup = warmup,
    cores = cores,
    backend = "cmdstanr",
    seed = seed,
    control = control,
    save_pars = save_pars(all = TRUE),
    refresh = refresh
  )

  saveRDS(fit, fit_path)
} else {
  if (!file.exists(fit_path)) {
    stop("No saved fit was found at: ", fit_path)
  }
  fit <- readRDS(fit_path)
  message("Loaded saved fit: ", fit_path)
}

workbook <- createWorkbook()

write_sheet(workbook, "missingness_report", missingness_report)
write_sheet(workbook, "prefit_validation", prefit_validation)

all_parameters <- summarise_all_draws(fit) |>
  rename(param = variable) |>
  arrange(desc(rhat))
readr::write_csv(all_parameters, all_parameters_path)
write_sheet(workbook, "all_parameters", all_parameters)

fixed_effects <- all_parameters |> filter(str_detect(param, "^b_"))
thresholds <- all_parameters |> filter(str_detect(param, "^Intercept\\["))
hyperparameters <- all_parameters |>
  filter(!str_detect(param, "^b_"), !str_detect(param, "^Intercept\\["), param != "lp__")

readr::write_csv(fixed_effects, fixed_effects_path)
readr::write_csv(thresholds, thresholds_path)
readr::write_csv(hyperparameters, hyperparameters_path)

write_sheet(workbook, "fixed_effects", fixed_effects)
write_sheet(workbook, "thresholds", thresholds)
write_sheet(workbook, "hyperparameters", hyperparameters)

fixed_effects_or <- fixed_effects |>
  mutate(
    OR = exp(mean),
    OR_lo = exp(q2.5),
    OR_hi = exp(q97.5)
  )
readr::write_csv(fixed_effects_or, fixed_or_path)
write_sheet(workbook, "fixed_effects_OR", fixed_effects_or)

site_random_intercepts <- extract_ranef_table(fit, "SiteID") |>
  rename(SiteID = level) |>
  arrange(desc(mu))
sequence_random_intercepts <- extract_ranef_table(fit, "SeqSegID") |>
  rename(SeqSegID = level) |>
  arrange(desc(mu))

readr::write_csv(site_random_intercepts, site_re_path)
readr::write_csv(sequence_random_intercepts, sequence_re_path)
write_sheet(workbook, "site_random_intercepts", site_random_intercepts)
write_sheet(workbook, "sequence_random_intercepts", sequence_random_intercepts)

readr::write_csv(data_overview, data_overview_path)
readr::write_csv(sequence_overview, sequence_overview_path)
write_sheet(workbook, "data_overview", data_overview)
write_sheet(workbook, "sequence_overview", sequence_overview)
write_sheet(workbook, "model_info", model_info)

saveWorkbook(workbook, xlsx_path, overwrite = TRUE)

writeLines(
  c(
    paste("Model tag:", model_tag),
    paste("Data path:", paths$data_path),
    paste("Fit path:", fit_path),
    paste("Prep path:", prep_path),
    paste("Rows used:", nrow(model_data)),
    paste("Segments:", n_distinct(model_data$SeqSegID))
  ),
  con = fit_runlog_path
)
capture.output(sessionInfo(), file = fit_session_path)

message("Done. Outputs saved under: ", paths$run_dir)
