#' Visualize Stan model fit results
#'
#' Creates comprehensive visualizations of the fitted Stan model showing observed
#' data, nowcasts, and forecasts for cases, ETU occupancy, alerts, and isolation
#' occupancy. The plot includes credible intervals and distinguishes between
#' observed and projected data.
#'
#' @param results A list containing the fitted Stan model results from fit_stan().
#' @param ylim_factor Numeric factor for controlling y-axis limits to prevent
#'   extreme values from dominating the plot (default: 4).
#' @param base_size Numeric base font size for the plot (default: 12).
#' @param quantiles Numeric vector specifying the credible interval levels to
#'   display (default: c(0.25, 0.5, 0.75, 0.95)).
#'
#' @return A ggplot object showing the model fit visualization.
#'
#' @details The visualization includes:
#' \itemize{
#'   \item Observed data as bar charts
#'   \item Nowcasts and forecasts as ribbon plots with credible intervals
#'   \item Different alpha levels for different credible interval widths
#'   \item Vertical line separating observed from projected data
#'   \item Faceted plots for cases, ETU occupancy, alerts, and isolation occupancy
#' }
#'
#' The function automatically handles extreme values by capping them at a
#' reasonable multiple of the observed maximum to ensure readable plots.
#'
#' @examples
#' \dontrun{
#' # Create visualization with default settings
#' p <- vis_stan_fit(results)
#'
#' # Customize visualization
#' p <- vis_stan_fit(results, ylim_factor = 2, base_size = 14)
#'
#' # Save the plot
#' p %>% save_plot("model_fit.png", height = 20)
#' }
#'
#' @importFrom dplyr map_dfr bind_rows group_by summarise mutate filter
#' @importFrom purrr map2_dfr imap_dfr
#' @importFrom ggplot2 ggplot aes geom_col geom_ribbon geom_line geom_vline
#' @importFrom ggplot2 scale_linetype scale_alpha_discrete facet_wrap labs
#' @importFrom ggplot2 scale_x_date scale_y_continuous theme_minimal theme
#' @importFrom ggplot2 element_rect guide_legend
#' @importFrom glue glue
#' @importFrom tibble tibble
#' @importFrom forcats fct_rev
#' @importFrom scales percent
#' @export
vis_stan_fit <- function(results,
                         ylim_factor = 4,
                         base_size = 12,
                         quantiles = c(0.25, 0.5, 0.75, 0.95)) {
  ## function to manually remove too large values
  get_max <- function(upper, observed, mid) {
    ## if no reported deaths, set max val per million pop
    if (all(is.na(observed))) {
      max_val <- max(mid)
    } else {
      max_val <- max(observed, na.rm = TRUE)
    }
    if (max_val == 0) max_val <- 1e7

    replace(
      upper,
      ifelse(is.na(upper), FALSE, upper > ylim_factor * max_val),
      if (any(upper < ylim_factor * max_val, na.rm = TRUE)) {
        max(upper[upper < ylim_factor * max_val], na.rm = TRUE)
      } else {
        ylim_factor * max_val
      }
    )
  }

  df <- map_dfr(
    c("etu", "iso", "alerts", "cases"),
    ~ tibble(
      day = seq_len(results$data$n_days),
      observed = pluck(results, "data", glue("{.x}_observed")),
      what = .x,
      inflated = FALSE
    )
  ) %>%
    full_join(
      bind_rows(
        map2_dfr(
          c(
            "log_cases_inflated", "log_cases_projected", "etu_modelled",
            "alerts_modelled", "iso_modelled"
          ),
          c("cases", "cases", "etu", "alerts", "iso"),
          function(varname, var) {
            quantiles %<>% sort()
            breaks <- sort(as_vector(map(quantiles, ~ c(0.5 - .x / 2, 0.5 + .x / 2))))
            values <- summarise_stan(varname, results, breaks)
            mid <- summarise_stan(varname, results, 0.5)
            imap_dfr(
              c(rev(quantiles), tail(quantiles, -1)),
              ~ tibble(
                quantile = .x,
                group = letters[.y],
                day = values$day,
                lower = values[[.y + 1]],
                mid = mid[[2]],
                upper = values[[.y + 2]]
              )
            ) %>%
              mutate(what = var, inflated = TRUE)
          }
        ),
        summarise_stan("log_cases_fitted", results) %>%
          mutate(what = "cases", inflated = FALSE, lower = NA, upper = NA, quantile = 1)
      ),
      by = c("day", "what", "inflated")
    ) %>%
    group_by(what) %>%
    mutate(
      upper = get_max(upper, observed, mid),
      mid = replace(mid, mid > max(upper, na.rm = TRUE), NA),
      lower = replace(lower, lower > max(upper, na.rm = TRUE), NA),
      what = factor(what, c("cases", "etu", "alerts", "iso")),
      date = as.Date(day, origin = min(results$data$date) - 1)
    )

  df %>%
    ggplot(aes(date)) +
    geom_col(aes(y = observed)) +
    geom_ribbon(
      aes(
        ymin = lower, ymax = upper,
        alpha = fct_rev(factor(quantile)),
        group = group
      ),
      fill = "firebrick",
      color = NA,
      data = filter(df, !is.na(group))
    ) +
    geom_line(
      aes(y = mid, linetype = fct_rev(factor(as.numeric(inflated)))),
      data = filter(df, quantile == quantiles[1] | is.na(group)),
      color = "black",
      size = 1
    ) +
    geom_vline(xintercept = max(results$data$date) + 0.5, linetype = 3) +
    scale_linetype(
      labels = c("1" = "Nowcast", "0" = "Reported")
    ) +
    scale_alpha_discrete(
      range = c(0.25, 0.8),
      labels = ~ scales::percent(as.numeric(as.character(.x)), 1)
    ) +
    facet_wrap(
      ~what,
      ncol = 1,
      strip.position = "left",
      scales = "free_y",
      labeller = labeller(what = c(
        cases = "Daily number of cases",
        etu = "Daily ETU occupancy",
        alerts = "Daily number of alerts",
        iso = "Daily isolation occupancy"
      ))
    ) +
    labs(
      x = NULL,
      y = NULL,
      alpha = "Credible interval",
      linetype = "Type"
    ) +
    guides(alpha = guide_legend(reverse = TRUE)) +
    scale_x_date(
      expand = c(0, 0),
      date_labels = "%b %d",
      date_breaks = "2 weeks"
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_minimal(base_size) +
    theme(
      strip.placement = "outside",
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "bottom"
    )
}
