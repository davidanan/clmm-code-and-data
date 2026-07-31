# Software Environment

The release scripts were checked in the environment below.

## Tested Platform

- Operating system: macOS Sequoia 15.6.1
- Architecture: Apple silicon (`aarch64-apple-darwin20`)
- R version: 4.5.2
- CmdStan version: 2.38.0

Other recent operating systems and R versions may work but were not tested for this
release.

## R Packages

| Package | Tested version |
| --- | --- |
| `dplyr` | 1.2.0 |
| `tidyr` | 1.3.2 |
| `stringr` | 1.6.0 |
| `purrr` | 1.2.1 |
| `readr` | 2.1.6 |
| `tibble` | 3.3.1 |
| `brms` | 2.23.0 |
| `cmdstanr` | 0.8.0 |
| `posterior` | 1.6.1 |
| `openxlsx` | 4.2.8.1 |
| `bayesplot` | 1.15.0 |
| `ggplot2` | 4.0.2 |

## Installation

Install the required R packages from CRAN and the Stan package repository:

```r
install.packages(c(
  "dplyr", "tidyr", "stringr", "purrr", "readr", "tibble",
  "brms", "posterior", "openxlsx", "bayesplot", "ggplot2"
))

install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", "https://cloud.r-project.org")
)
```

Install and verify the CmdStan version used for this release:

```r
cmdstanr::install_cmdstan(version = "2.38.0")
cmdstanr::cmdstan_version()
cmdstanr::check_cmdstan_toolchain()
```

The installation of CmdStan requires a working C++ toolchain. See the official
`cmdstanr` installation documentation if toolchain configuration is required for the
operating system.

## Computational Requirements

The released fits use four Markov chain Monte Carlo chains with 3,000 iterations per
chain, including 1,500 warmup iterations. Full fitting is computationally intensive,
and wall-clock time depends strongly on the processor, compiler, and available cores.
No portable runtime benchmark is reported. Four available CPU cores and 16 GB of RAM
are a practical starting configuration.

Before a full fit, run the preprocessing-only checks:

```bash
CLMM_PREP_ONLY=TRUE Rscript code/01_fit_phase_diagnostic_clmm.R
CLMM_PREP_ONLY=TRUE Rscript code/01_fit_conditional_discrepancy_clmm.R
```

## Notes

- fitting uses `brms` with the `cmdstanr` backend;
- evaluation scripts can read model objects from one results directory and write
  evaluation outputs to another;
- all default file references use paths relative to the repository root;
- numerical results may vary slightly across platforms or package versions.
