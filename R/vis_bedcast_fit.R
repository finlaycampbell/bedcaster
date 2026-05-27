#' Visualize Stan model fit results
#'
#' Creates comprehensive visualizations of the fitted Stan model showing reported
#' data, nowcasts, and forecasts for cases, ETU occupancy, alerts, and isolation
#' occupancy. The plot includes credible intervals and distinguishes between
#' reported and projected data.
#'
#' @param results A \code{bedcast} object from \code{\link{fit_bedcaster}}.
#' @param ylim_factor Numeric factor for controlling y-axis limits to prevent
#'   extreme values from dominating the plot (default: 4).
#' @param base_size Numeric base font size for the plot (default: 12).
#' @param quantiles Numeric vector specifying the credible interval levels to
#'   display (default: c(0.25, 0.5, 0.75, 0.95)).
#'
#' @return A ggplot object showing the model fit visualization.
#'
#' @importFrom dplyr bind_rows group_by summarise mutate filter transmute
#' @importFrom purrr map_dfr map imap_dfr set_names pluck as_vector
#' @importFrom ggplot2 ggplot aes geom_col geom_ribbon geom_line geom_vline
#' @importFrom ggplot2 scale_linetype scale_alpha_continuous facet_wrap labs
#' @importFrom ggplot2 scale_x_date scale_y_continuous theme_minimal theme
#' @importFrom ggplot2 element_rect guide_legend guides labeller
#' @importFrom glue glue
#' @importFrom tibble tibble
#' @importFrom forcats fct_inorder fct_rev
#' @importFrom tidyr drop_na
#' @importFrom stringr str_split
#' @importFrom stats quantile
#' @importFrom utils tail
#' @export
#'
vis_bedcast_fit <- function(results,
                            ylim_factor = 4,
                            base_size = 12,
                            quantiles = c(0.25, 0.5, 0.75, 0.95)) {

  # function to manually remove too large values
  set_max <- function(vals, max_val) {
    if (!is.finite(max_val) || max_val <= 0) {
      return(vals)
    }
    limit <- max_val * ylim_factor
    replace(vals, !is.na(vals) & vals > limit, limit)
  }

  # anchor for ylim_factor capping
  get_max <- function(var) {
    reported <- dplyr::filter(data, what == var)$reported
    reported_max <- max(reported, na.rm = TRUE)
    if (is.finite(reported_max) && reported_max > 0) {
      return(reported_max)
    }

    # no reported data: anchor on latest central estimate in fit window
    sim_par <- paste0(var, "_nowcast_sim")
    sim_med <- summary(results, sim_par, probs = 0.5)
    fit_dates <- results$data$date_fit
    sim_fit <- sim_med[sim_med$index %in% fit_dates, , drop = FALSE]

    if (nrow(sim_fit) == 0) {
      return(max(sim_med$q_0.5, na.rm = TRUE))
    }

    latest_idx <- which.max(sim_fit$index)
    anchor <- sim_fit$q_0.5[latest_idx]

    if (!is.finite(anchor) || anchor <= 0) {
      anchor <- max(sim_fit$q_0.5, na.rm = TRUE)
    }

    anchor
  }

  # define fitting variables
  fitting_vars <- c("cases_nowcast_sim", "cases_projected_sim",
                     "deaths_nowcast_sim", "deaths_projected_sim",
                     "etu_nowcast_sim", "alerts_nowcast_sim",
                     "iso_nowcast_sim")

  # define data variables
  data_vars <- c(
    cases = "Daily reported cases",
    deaths = "Daily reported deaths",
    etu = "ETU occupancy",
    iso = "Isolation occupancy",
    alerts = "Daily alerts"
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
    mutate(what = forcats::fct_inorder(what))

  # extract quantiles for ribbons
  ribbons <- map_dfr(
    fitting_vars,
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
          upper = values[[.y + 3]],
          inflated = str_split(varname, "_")[[1]][2] != "truncated"
        )
      ) |>
        mutate(
          upper = set_max(upper, get_max(var)),
          lower = replace(
            lower,
            lower > max(upper, na.rm = TRUE),
            max(upper, na.rm = TRUE)
          ),
          what = var
        )
    }
  ) |>
    mutate(what = factor(what, names(data_vars)))

  # extract median estimates for lines
  mids <- map_dfr(
    c(fitting_vars, "cases_truncated_sim", "deaths_truncated_sim"),
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
      data = drop_na(data, reported),
      aes(index, reported)
    ) +
    geom_ribbon(
      data = drop_na(ribbons, lower, upper),
      aes(
        index, ymin = lower, ymax = upper,
        alpha = quantile,
        ## alpha = fct_rev(factor(quantile)),
        group = group
      ),
      fill = "firebrick",
      color = NA
    ) +
    geom_line(
      data = mids,
      aes(index, mid, linetype = forcats::fct_rev(factor(as.numeric(inflated)))),
      color = "black",
      linewidth = 1
    ) +
    geom_vline(
      xintercept = max(results$data$date_fit) + 0.5,
      linetype = 3
    ) +
    scale_linetype(
      labels = c("1" = "Nowcast", "0" = "Reported")
    ) +
    scale_alpha_continuous(
      range = c(0.9, 0.25),
      labels = ~ scales::percent(as.numeric(as.character(.x)), 1),
      breaks = quantiles
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
    guides(alpha = guide_legend(reverse = FALSE)) +
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
