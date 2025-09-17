# bedcast

Bed Occupancy Forecasting for Epidemics

## Overview

The `bedcast` package provides tools for forecasting bed occupancy during epidemics using Bayesian nowcasting methods. It is designed to help healthcare systems prepare for and respond to epidemic outbreaks by providing accurate forecasts of ETU (Ebola Treatment Unit) and isolation bed demand.

## Features

- **Data Import and Processing**: Import and clean epidemiological data from multiple sources
- **Delay Estimation**: Calculate delay distributions between key epidemiological events
- **Bayesian Nowcasting**: Account for reporting delays in case data
- **Stan Integration**: Leverage Stan for efficient MCMC sampling
- **Comprehensive Visualization**: Create publication-ready plots of results
- **Uncertainty Quantification**: Provide credible intervals for all estimates

## Installation

```r
# Install from source
devtools::install_github("your-username/bedcast")

# Or install locally
devtools::install("path/to/bedcast")
```

## Quick Start

```r
library(bedcast)

# Import data
ll <- get_real_linelist()
alerts <- get_alerts()
etu_iso <- get_etu_iso()
delays <- get_delays(ll)

# Merge data
data <- merge_data(ll, alerts, etu_iso, delays)

# Compile and fit model
model <- stan_model("inst/stan/occupancy_model_nowcast.stan")
results <- fit_stan(data, model)

# Visualize results
vis_stan_fit(results) %>% save_plot("model_fit.png")
vis_parameters(results) %>% save_plot("parameters.png")
vis_growth_rate(results) %>% save_plot("growth_rate.png")
```

## Workflow

1. **Data Import**: Import linelist, alerts, and occupancy data
2. **Delay Estimation**: Calculate delay distributions from linelist data
3. **Data Merging**: Combine all data sources into time series format
4. **Model Fitting**: Fit Bayesian Stan model for nowcasting and forecasting
5. **Visualization**: Create comprehensive plots of results

## Documentation

- See the [workflow vignette](vignettes/bedcast-workflow.Rmd) for a complete example
- Function documentation is available via `?function_name`
- The package follows the tidyverse style guide

## Dependencies

The package depends on:
- `tidyverse` for data manipulation
- `rstan` for Bayesian modeling
- `tsibble` for time series data
- `EpiNow2` and `simulacr` for epidemiological functions
- `ggplot2` for visualization

## License

MIT License

## Contributing

Contributions are welcome! Please feel free to submit issues and pull requests.

## Citation

If you use this package in your research, please cite it as:

```
Bed Occupancy Forecasting for Epidemics (bedcast) R Package
https://github.com/your-username/bedcast
```
