# bedcaster

**Bed occupancy forecasting for epidemics**

`bedcaster` fits a Bayesian Stan model to daily surveillance data and produces nowcasts and short-term projections of cases, deaths, ETU occupancy, community alerts, and isolation bed use. It is designed for outbreak response when reporting delays and partial observation of healthcare indicators need to be accounted for explicitly.

## Features

- Joint modelling of cases, deaths, ETU beds, alerts, and isolation occupancy
- Reporting-delay adjustment for case and death series
- Spline-based time-varying growth rates with optional extrapolation over the projection window
- Informative priors on delays, case fatality, and alert processes
- `bedcast` objects with S3 methods for printing, summarising, extracting, and plotting
- ggplot2 visualisations for fit, parameters, and growth rate

## Installation

```r
# install.packages("remotes")
remotes::install_github("finlaycampbell/bedcaster")
```

To install from a local clone:

```r
remotes::install_local("path/to/bedcaster")
```

**Requirements:** R (>= 4.1.0), a C++ toolchain, and [rstan](https://mc-stan.org/rstan/) with StanHeaders. On Windows, install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) before building.

## Input data

Prepare a data frame with one row per day and these columns:

| Column   | Description                                      |
|----------|--------------------------------------------------|
| `date`   | Date of observation                              |
| `cases`  | Reported cases (use `NA` if not available)       |
| `deaths` | Reported deaths                                  |
| `etu`    | ETU occupancy                                    |
| `alerts` | New community alerts                             |
| `iso`    | Isolation bed occupancy                          |

Series that are not recorded every day can contain `NA`; the model treats missing values accordingly.

## Quick start

```r
library(bedcaster)

# df: daily data frame as described above
as_of <- max(df$date, na.rm = TRUE)

bedcast <- fit_bedcaster(
  df,
  as_of = as_of,
  n_proj = 28,
  n_iter = 1000,
  n_chains = 4,
  n_cores = 4
)

print(bedcast)

# Summarise a Stan quantity (e.g. nowcasted cases)
summary(bedcast, "cases_nowcast_sim", alpha = 0.95)

# Plots
plot(bedcast, type = "fit")
plot(bedcast, type = "parameters")
plot(bedcast, type = "growthrate")

# Or call visualisation functions directly
vis_bedcast_fit(bedcast)
vis_bedcast_parameters(bedcast)
vis_bedcast_growthrate(bedcast)

# Extract MCMC draws for custom analysis
extract(bedcast, "cfr")
```

## Priors

Delay priors are length-4 named vectors (`meanlog.mean`, `meanlog.sd`, `sdlog.mean`, `sdlog.sd`) passed to the Stan model. The onset-to-reporting delay used for nowcasting uses `meanlog.mean` and `sdlog.mean` with [distcrete](https://cran.r-project.org/package=distcrete).

```r
delay_prior <- c(
  meanlog.mean = log(5), meanlog.sd = 0.5,
  sdlog.mean = log(2), sdlog.sd = 0.5
)

bedcast <- fit_bedcaster(
  df,
  as_of = as_of,
  prior_onset_to_etu = delay_prior,
  prior_cfr = c(qlogis(0.3), 0.25),
  extrapolate_growthrate = FALSE
)
```

See `?fit_bedcaster` for all prior and sampling arguments.

## Workflow

1. **Prepare data** — Build a daily time series from linelists, bed counts, and alert systems.
2. **Fit** — Call `fit_bedcaster()` with an `as_of` date matching your surveillance snapshot.
3. **Review** — Use `print()`, `summary()`, and `plot()` to assess fit and parameters.
4. **Export** — Use `extract()` for draws or `summary()` for quantiles by date.

A longer walkthrough is in the package vignette:

```r
vignette("bedcaster-workflow", package = "bedcaster")
```

## Main functions

| Function | Purpose |
|----------|---------|
| `fit_bedcaster()` | Fit the Stan model and return a `bedcast` object |
| `summary()` | Posterior quantiles for a Stan parameter |
| `extract()` | MCMC draws in long format |
| `plot()` | Dispatch to fit, parameter, or growth-rate plots |
| `vis_bedcast_fit()` | Cases, deaths, ETU, alerts, and isolation fit plot |
| `vis_bedcast_parameters()` | Prior vs posterior parameter densities |
| `vis_bedcast_growthrate()` | Time-varying growth rate with intervals |
| `get_auc()` | Normalisation helper for prior density plots |

## License

MIT — see [LICENSE](LICENSE).

## Citation

```
Campbell F (2026). bedcaster: Bed Occupancy Forecasting for Epidemics.
R package version 0.1.0. https://github.com/finlaycampbell/bedcaster
```
