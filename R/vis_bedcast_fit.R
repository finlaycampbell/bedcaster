#' Visualize Stan model fit results
#'
#' Creates comprehensive visualizations of the fitted Stan model showing reported
#' data, nowcasts, and forecasts for cases, ETU occupancy, alerts, and isolation
#' occupancy. The plot includes credible intervals and distinguishes between
#' reported and projected data.
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
#' The function automatically handles extreme values by capping them at a
#' reasonable multiple of the reported maximum to ensure readable plots.
#'
#' @examples
#' \dontrun{
#' # Create visualization with default settings
#' p <- vis_bedcast_fit(results)
#'
#' # Customize visualization
#' p <- vis_bedcast_fit(results, ylim_factor = 2, base_size = 14)
#'
#' # Save the plot
#' p %>% save_plot("model_fit.png", height = 20)
#' }
#'
#' @importFrom dplyr bind_rows group_by summarise mutate filter
#' @importFrom purrr map_dfr
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
vis_bedcast_fit <- function(results,
                            ylim_factor = 4,
                            base_size = 12,
                            quantiles = c(0.25, 0.5, 0.75, 0.95)) {

  # function to manually remove too large values
  set_max <- function(vals, max_val) {
    replace(
      vals,
      ifelse(is.na(vals), FALSE, vals > ylim_factor * max_val),
      if (any(vals < ylim_factor * max_val, na.rm = TRUE)) {
        max(vals[vals < ylim_factor * max_val], na.rm = TRUE)
      } else {
        ylim_factor * max_val
      }
    )
  }

  # get maximum reported for a given variable
  get_max <- function(var) {
    max(filter(data, what == var)$reported, na.rm = TRUE)
  }

  # define fitting variables
  fittings_vars <- c("cases_nowcast_sim", "cases_projected_sim",
                     "deaths_nowcast_sim", "deaths_projected_sim",
                     "etu_nowcast_sim", "alerts_nowcast_sim",
                     "iso_nowcast_sim")

  # define data variables
  data_vars <- c(
    cases = "Daily number of cases",
    deaths = "Daily number of deaths",
    etu = "Daily ETU occupancy",
    iso = "Daily isolation occupancy",
    alerts = "Daily number of alerts"
  )

  # extract data
  data <- map_dfr(
    purrr::set_names(names(data_vars)),
    ~ tibble(
      index = results$data$date_fit,
      reported = pluck(results, "data", glue("{.x}_reported"))
    ),
    .id = "what"
  ) |>
    mutate(what = fct_inorder(what))

  # extract quantiles for ribbons
  ribbons <- map_dfr(
    vars,
    function(varname) {
      quantiles <- sort(quantiles)
      breaks <- sort(
        as_vector(map(quantiles, ~ c(0.5 - .x / 2, 0.5 + .x / 2)))
      )
      values <- summary(results, varname, probs = breaks)
      var <- str_split(varname, "_")[[1]][1]
      imap_dfr(
        c(rev(quantiles), tail(quantiles, -1)),
        ~ tibble(
          quantile = .x,
          group = letters[.y],
          index = values$index,
          lower = values[[.y + 2]],
          upper = values[[.y + 3]]
        )
      ) %>%
        mutate(
          upper = set_max(upper, get_max(var)),
          lower = replace(lower, lower > max(upper, na.rm = TRUE), NA),
          what = var,
          inflated = str_split(varname, "_")[[1]][2] != "truncated"
        )
    }
  ) |>
    mutate(what = factor(what, names(data_vars)))

  # extract median estimates for lines
  mids <- map_dfr(
    c(vars, "cases_truncated_sim", "deaths_truncated_sim"),
    function(varname) {
      var <- str_split(varname, "_")[[1]][1]
      summary(results, varname, probs = 0.5) |>
        transmute(
          index,
          mid = set_max(q_0.5, get_max(var)),
          what = var,
          inflated = str_split(varname, "_")[[1]][2] != "truncated"
        )
    }
  ) |>
    mutate(what = factor(what, names(data_vars)))

  # plot
  ggplot() +
    geom_col(
      data = data,
      aes(index, reported)
    ) +
    geom_ribbon(
      data = ribbons,
      aes(
        index, ymin = lower, ymax = upper,
        alpha = fct_rev(factor(quantile)),
        group = group
      ),
      fill = "firebrick",
      color = NA
    ) +
    geom_line(
      data = mids,
      aes(index, mid, linetype = fct_rev(factor(as.numeric(inflated)))),
      color = "black",
      size = 1
    ) +
    geom_vline(
      xintercept = max(results$data$date_fit) + 0.5,
      linetype = 3
    ) +
    scale_linetype(
      labels = c("1" = "Nowcast", "0" = "Reported")
    ) +
    scale_alpha_discrete(
      range = c(0.25, 0.8),
      labels = ~ scales::percent(as.numeric(as.character(.x)), 1)
    ) +
    facet_wrap(
      ~ what,
      ncol = 1,
      strip.position = "left",
      scales = "free_y",
      labeller = labeller(what = data_vars)
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
