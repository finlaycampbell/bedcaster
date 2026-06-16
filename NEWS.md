# bedcaster (development version)

## Bug fixes

- Death likelihood is now skipped for days where `deaths` is `NA`, matching ETU,
  alerts, and isolation handling. Previously `NA` deaths were passed as `-1000`
  and still included in the vectorised death likelihood, causing sampler failures.

# bedcaster 0.1.0

Initial release.

## Features

- `fit_bedcaster()` — Bayesian Stan model for cases, deaths, ETU occupancy, alerts, and isolation beds with reporting-delay adjustment and configurable projection horizon.
- `bedcast` S3 class with `print()`, `summary()`, `extract()`, and `plot()` methods.
- Visualisation helpers: `vis_bedcast_fit()`, `vis_bedcast_parameters()`, `vis_bedcast_growthrate()`.
- `get_auc()` for prior density normalisation in parameter plots.
- Vignette: *Bed Occupancy Forecasting Workflow*.

## Notes

- Compiled Stan model ships with the package via rstantools.
- Requires R >= 4.1.0, rstan, and a C++ toolchain for installation from source.
