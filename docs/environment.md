# Running the R scripts

The scripts were checked on Apple silicon with R 4.5.2, cmdstanr 0.8.0 and CmdStan 2.38.0. They use base R to read and prepare the data and cmdstanr to compile and sample the embedded Stan models.

## Setup

Install R and a C++ toolchain for CmdStan. In an R session, run:

```r
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", "https://cloud.r-project.org"))
cmdstanr::check_cmdstan_toolchain()
cmdstanr::install_cmdstan(version = "2.38.0")
```

If CmdStan is already installed, check the installation with `cmdstanr::cmdstan_version()`. For a non-default location, set the `CMDSTAN` environment variable to the installation directory before starting R.

## Sampling

Each script runs four parallel chains with 1,500 adaptive warmup iterations and 2,000 retained draws per chain. The target acceptance probability is 0.95 and the maximum tree depth is 12. These scripts use the manuscript's likelihood and priors and tune the sampler at each new fit. The manuscript runs used chain-specific settings from pilot fits and 200 further warmup iterations at fixed settings. Individual posterior draws will therefore differ.

The standalone entry points were checked by compilation and short sampling runs. For a full refit, inspect sampler diagnostics, R-hat and effective sample sizes before summarising the posterior.
