#' Visualize parameter distributions and priors
#'
#' Creates visualizations comparing posterior parameter distributions with their
#' prior distributions. This helps assess how much the data informed the parameter
#' estimates and whether the priors were appropriate.
#'
#' @param results A \code{bedcast} object from \code{\link{fit_bedcaster}}.
#' @param base_size Numeric base font size for the plot (default: 12).
#'
#' @return A ggplot object showing parameter distributions with priors overlaid.
#'
#' @importFrom dplyr group_by summarise mutate select bind_rows rename left_join
#' @importFrom purrr map_dfr pmap set_names
#' @importFrom tidyr unnest
#' @importFrom ggplot2 ggplot aes geom_area geom_density facet_wrap labs
#' @importFrom ggplot2 scale_y_continuous scale_x_continuous theme_minimal theme
#' @importFrom ggplot2 element_rect labeller
#' @importFrom tibble tibble
#' @importFrom stats qnorm dnorm plogis qlogis density rnorm
#' @importFrom stringr str_remove
#' @export
#'
vis_bedcast_parameters <- function(results, base_size = 12) {

  # define plot labels
  labels <- c(
    onset_to_reporting_logmean = "Delay onset to reporting",
    onset_to_etu_logmean = "Delay onset to hospitalisation",
    etu_to_survival_logmean = "Delay hospitalisation to survival",
    etu_to_death_logmean = "Delay hospitalisation to death",
    onset_to_iso_logmean = "Delay onset to contact isolation",
    iso_to_release_logmean = "Lab turnaround time",
    cfr = "In-hospital case fatality rate",
    prop_iso = "Proportion of alerts isolated",
    alerts_per_case = "Number of alerts per case",
    alerts_background = "Background number of alerts"
  )

  delay_prefixes <- c(
    "onset_to_reporting",
    "onset_to_etu",
    "etu_to_survival",
    "etu_to_death",
    "onset_to_iso",
    "iso_to_release"
  )

  # delay posteriors on the lognormal mean: exp(meanlog + sdlog^2 / 2)
  delay_samples <- map_dfr(delay_prefixes, function(prefix) {
    ml <- extract(results, paste0(prefix, "_logmean"))
    sl <- extract(results, paste0(prefix, "_sd"))
    ml |>
      rename(meanlog = value) |>
      left_join(
        sl |> rename(sdlog = value),
        by = c("index", "iter")
      ) |>
      mutate(
        var = paste0(prefix, "_logmean"),
        value = exp(meanlog + sdlog^2 / 2)
      ) |>
      select(var, index, iter, value)
  })

  other_vars <- setdiff(names(labels), unique(delay_samples$var))
  other_samples <- map_dfr(
    set_names(other_vars),
    ~ extract(results, .x),
    .id = "var"
  )

  samples <- bind_rows(delay_samples, other_samples)

  # generate priors
  prior <- samples |>
    group_by(var) |>
    summarise(minval = min(value), maxval = max(value), .groups = "drop") |>
    mutate(
      prior = pmap(
        list(var, minval, maxval),
        function(v, mi, ma) make_prior(results, v, mi, ma)
      )
    ) |>
    select(var, prior) |>
    unnest(prior) |>
    mutate(var = factor(var, names(labels)), group = "AA")

  # generate plot
  ggplot() +
    geom_area(
      data = prior,
      aes(value, y = density),
      fill = "grey20", alpha = 0.75, color = "black"
    ) +
    geom_density(
      data = samples,
      aes(value, group = factor(index)),
      fill = "darkgreen", alpha = 0.5,
    ) +
    facet_wrap(~ var, scales = "free", labeller = labeller(var = labels)) +
    scale_y_continuous(expand = c(0.01, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(
      x = "Value",
      y = "Density",
      caption = "* Green represents the posterior estimate informed by data, grey the prior distribution"
    ) +
    theme_minimal(base_size) +
    theme(plot.background = element_rect(fill = "white", color = NA))
}


#' Build prior density curves for parameter plots
#'
#' @keywords internal
#' @noRd
make_prior <- function(results, var, minval, maxval) {
  varname <- paste0("prior_", str_remove(var, "_logmean"))
  prior <- results$data[[varname]]
  if (varname %in% c("prior_cfr", "prior_prop_iso")) {
    rnge <- qnorm(c(0.01, 0.99), prior[1], prior[2])
    sq <- seq(
      min(rnge[1], qlogis(minval)),
      max(rnge[2], qlogis(maxval)),
      length = 100
    )
    tibble(
      value = plogis(sq),
      density = dnorm(sq, prior[1], prior[2]) / get_auc(prior, "plogis")
    )
  } else if (varname == "prior_growth_rate") {
    rnge <- qnorm(c(0.01, 0.99), prior[1], prior[2])
    sq <- seq(min(rnge[1], minval), max(rnge[2], maxval), length = 100)
    tibble(
      value = sq,
      density = dnorm(sq, prior[1], prior[2]) / get_auc(prior, "identity")
    )
  } else if (grepl("alerts", varname)) {
    rnge <- qnorm(c(0.01, 0.99), prior[1], prior[2])
    sq <- seq(
      min(rnge[1], log(minval)),
      max(rnge[2], log(maxval)),
      length = 100
    )
    tibble(
      value = exp(sq),
      density = dnorm(sq, prior[1], prior[2]) / get_auc(prior, "exp")
    )
  } else if (length(prior) >= 4L && grepl("_logmean$", var)) {
    # Pathway delay: normal priors on meanlog and sdlog -> lognormal mean
    n <- 2000L
    ml <- rnorm(n, prior[1], prior[2])
    sl <- pmax(rnorm(n, prior[3], prior[4]), 1e-6)
    mean_nat <- exp(ml + sl^2 / 2)
    from <- min(minval, stats::quantile(mean_nat, 0.001))
    to <- max(maxval, stats::quantile(mean_nat, 0.999))
    dens <- density(mean_nat, from = from, to = to)
    tibble(value = dens$x, density = dens$y)
  } else {
    rnge <- qnorm(c(0.01, 0.99), prior[1], prior[2])
    sq <- seq(
      min(rnge[1], minval),
      max(rnge[2], maxval),
      length = 100
    )
    tibble(
      value = exp(sq),
      density = dnorm(sq, prior[1], prior[2]) / get_auc(prior, "exp")
    )
  }
}
