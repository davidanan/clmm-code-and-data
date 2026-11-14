#!/usr/bin/env Rscript
# Flood-phase CLMM: observed inundation during rising, peak and recession.
# Run: Rscript code/fit_flood_phase_clmm.R [output_directory]
# Requires R, cmdstanr and an installed CmdStan toolchain.
# This standalone refit adapts its own sampler; no saved pilot states are required.

prepare_data <- function(csv) {
  records <- read.csv(csv)
  records <- records[order(records$SeqSegID, records$t), ]
  segment <- match(records$SeqSegID, unique(records$SeqSegID))
  site <- match(records$Site_ID, sort(unique(records$Site_ID)))
  event <- factor(records$Event_ID, levels = 1:5)
  stopifnot(!anyNA(event), all(records$Y_obs %in% 0:4), all(records$Y_sim %in% 0:4))
  starts <- which(!duplicated(segment))
  following <- which(duplicated(segment))
  stopifnot(all(records$t[following] - records$t[following - 1L] == 1L))

  # Event 1 (12 August 2023) is the reference storm.
  events <- model.matrix(~ event, contrasts.arg = list(event = "contr.treatment"))[, -1, drop = FALSE]
  phase <- factor(records$Flood_Phase, levels = c("rising", "peak", "recession"))
  stopifnot(!anyNA(phase))
  phases <- model.matrix(~ phase, contrasts.arg = list(phase = "contr.treatment"))[, -1, drop = FALSE]
  design <- cbind(phases, events)
  stopifnot(nrow(design) == nrow(records), all(is.finite(design)))
  # The Student-t threshold prior applies at the centred design, as in the fitted model.
  design <- scale(design, center = TRUE, scale = FALSE)
  list(N = nrow(records), K = ncol(design), X = design,
       observed = as.integer(records$Y_obs + 1L), simulated = as.integer(records$Y_sim + 1L),
       n_sites = max(site), site = site, n_segments = length(starts), segment = segment,
       first = array(starts), next_minute = array(following), previous_minute = array(following - 1L))
}

model_code <- "
data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  array[N] int<lower=1, upper=5> observed;
  array[N] int<lower=1, upper=5> simulated;
  int<lower=1> n_sites;
  array[N] int<lower=1, upper=n_sites> site;
  int<lower=1> n_segments;
  array[N] int<lower=1, upper=n_segments> segment;
  array[n_segments] int<lower=1, upper=N> first;
  array[N - n_segments] int<lower=2, upper=N> next_minute;
  array[N - n_segments] int<lower=1, upper=N> previous_minute;
}
parameters {
  ordered[4] cutpoints;
  vector[K] beta;
  real simulated_slope;
  simplex[4] class_weights;
  vector[n_sites] site_z;
  vector[n_segments] segment_z;
  real<lower=0> site_sd;
  real<lower=0> segment_sd;
  real<lower=-1, upper=1> rho;
  real<lower=0> innovation_sd;
  vector[N] state;  // AR(1) state in units of the innovation SD
}
model {
  vector[5] simulated_effect = append_row(0.0, 4 * simulated_slope * cumulative_sum(class_weights));
  vector[N] eta = X * beta + simulated_effect[simulated] + site_sd * site_z[site]
                  + segment_sd * segment_z[segment] + innovation_sd * state;

  cutpoints ~ student_t(3, 0, 2.5);
  beta ~ std_normal();
  simulated_slope ~ std_normal();
  class_weights ~ dirichlet(rep_vector(1.0, 4));
  site_z ~ std_normal();
  segment_z ~ std_normal();
  site_sd ~ student_t(3, 0, 1);
  segment_sd ~ student_t(3, 0, 1);
  innovation_sd ~ student_t(3, 0, 2.5);

  // Restart at stationarity for the first minute of each continuous video segment.
  state[first] ~ normal(0, inv_sqrt(1 - square(rho)));
  state[next_minute] ~ normal(rho * state[previous_minute], 1);
  observed ~ ordered_logistic(eta, cutpoints);
}
"

fit_model <- function(csv, output_dir, chains = 4L, warmup = 1500L, draws = 2000L) {
  data <- prepare_data(csv)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  stan_file <- cmdstanr::write_stan_file(model_code, dir = output_dir, basename = "flood_phase")
  model <- cmdstanr::cmdstan_model(stan_file)
  fit <- model$sample(data = data, seed = 202612051, chains = chains, parallel_chains = chains,
                      iter_warmup = warmup, iter_sampling = draws, adapt_delta = 0.95,
                      max_treedepth = 12, output_dir = output_dir)
  writeLines(colnames(data$X), file.path(output_dir, "beta_column_order.txt"))
  print(fit$diagnostic_summary())
  invisible(fit)
}

if (sys.nframe() == 0L) {
  script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  root <- dirname(dirname(normalizePath(script)))
  args <- commandArgs(trailingOnly = TRUE)
  output <- if (length(args)) args[[1]] else file.path(root, "results", "flood_phase")
  fit_model(file.path(root, "data", "paired_observed_simulated_inundation.csv"), output)
}
