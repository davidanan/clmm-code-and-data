# Urban flood model evaluation in Dalian

This repository accompanies *Sparse validation can mask process errors in urban flood modelling: evidence from large multimodal model-assisted video observations*, by Shunan Zhou, Heng Lyu, Anette Eltner, Pedro Zamboni, Nan Wang, Jiaxing Zhao and Chi Zhang.

Two Bayesian cumulative link mixed models (CLMMs) relate observed street inundation to simulated inundation during five storms in Dalian, China. The flood-phase model compares rising, peak and recession. The rainfall and site model examines associations with recent rainfall and local site characteristics. Both models account for differences among events, sites and video segments, and for correlation between consecutive minutes within a segment.

## Repository contents

| File | Purpose |
| --- | --- |
| [Flood-phase CLMM](code/fit_flood_phase_clmm.R) | Fit the phase model. The rising phase is the reference. |
| [Rainfall and site CLMM](code/fit_rainfall_site_clmm.R) | Fit the model with two rainfall descriptors and seven site attributes. |
| [Paired observations and simulations](data/paired_observed_simulated_inundation.csv) | Input to both scripts. Contains 20,825 paired minutes from 60 sites and 164 continuous video segments. |
| [Verified video observations](data/verified_video_inundation_levels.csv) | Full observation archive. Contains 36,490 verified minute records from 222 camera records at 124 sites. |

The [data dictionary](docs/data_dictionary.md) describes the columns, units, event dates and inundation scale.

## Fit the models

With R, cmdstanr and CmdStan installed, run the following commands from the repository root. See the [environment guide](docs/environment.md) for installation and sampler settings.

```sh
Rscript code/fit_flood_phase_clmm.R
Rscript code/fit_rainfall_site_clmm.R
```

Each script reads the paired CSV and creates its own directory under `results/`: `flood_phase/` or `rainfall_site/`. The directory contains the Stan source, compiled model, posterior chain CSVs and `beta_column_order.txt`, which identifies the predictor for each regression coefficient. The script prints sampler diagnostics after fitting. To choose another output directory, pass its path as the first argument.

## Citation and reuse

Please cite the manuscript and identify the code version used. The two CSVs retain the contents of the [original v1.0.0 data archive](https://doi.org/10.5281/zenodo.21728327) under descriptive filenames. The R scripts implement the manuscript models as standalone refits with fresh sampler adaptation.

The code is distributed under the [MIT License](LICENSE) and the derived data under [CC BY 4.0](DATA_LICENSE.md).
