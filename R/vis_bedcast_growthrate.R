#' Visualize growth rate over time
#'
#' Creates a visualization showing the estimated daily growth rate of cases over
#' time, including uncertainty bounds. This helps identify periods of increasing
#' or decreasing transmission.
#'
#' @param bedcast A \code{bedcast} object from \code{\link{fit_bedcaster}}.
#' @param alpha Coverage of the credible interval (default: 0.95).
#' @param base_size Numeric base font size for the plot (default: 12).
#'
#' @return A ggplot object showing the growth rate over time.
#'
#' @importFrom dplyr bind_rows mutate transmute .data
#' @importFrom ggplot2 ggplot aes geom_line geom_hline geom_ribbon
#' @importFrom ggplot2 scale_y_continuous scale_x_date scale_linetype_discrete
#' @importFrom ggplot2 labs theme_minimal theme element_rect
#' @importFrom stringr str_to_title
#' @export
#'
vis_bedcast_growthrate <- function(bedcast,
                                   alpha = 0.95,
                                   base_size = 12) {

  q_lo <- paste0("q_", format(0.5 - alpha / 2, scientific = FALSE))
  q_md <- paste0("q_", format(0.5, scientific = FALSE))
  q_hi <- paste0("q_", format(0.5 + alpha / 2, scientific = FALSE))

  bind_rows(
    summary(bedcast, "growthrate_reported", alpha = alpha) |>
      mutate(type = "observed"),
    summary(bedcast, "growthrate_projected", alpha = alpha) |>
      mutate(type = "projected")
  ) |>
    mutate(
      lower = .data[[q_lo]],
      mid = .data[[q_md]],
      upper = .data[[q_hi]]
    ) |>
    ggplot(aes(index, mid, ymin = lower, ymax = upper, linetype = type)) +
    geom_ribbon(alpha = 0.5) +
    geom_line(color = "firebrick", size = 1) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_date(
      expand = c(0, 0),
      date_labels = "%b %d",
      date_breaks = "2 weeks"
    ) +
    scale_linetype_discrete(labels = str_to_title) +
    labs(
      x = NULL,
      y = "Daily growth rate",
      linetype = NULL
    ) +
    theme_minimal(base_size) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "bottom"
    )

}
