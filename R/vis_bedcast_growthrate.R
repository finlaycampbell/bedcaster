#' Visualize growth rate over time
#'
#' Creates a visualization showing the estimated daily growth rate of cases over
#' time, including uncertainty bounds. This helps identify periods of increasing
#' or decreasing transmission.
#'
#' @param bedcast A list containing the output from \code{fit_bedcast}.
#' @param base_size Numeric base font size for the plot (default: 12).
#'
#' @return A ggplot object showing the growth rate over time.
#'
#' @importFrom dplyr group_by summarise
#' @importFrom ggplot2 ggplot aes geom_line geom_hline scale_y_continuous
#' @importFrom ggplot2 scale_x_date labs theme_minimal theme element_rect
#' @importFrom stringr str_to_title
#' @export
#'
vis_bedcast_growthrate <- function(bedcast,
                                   alpha = 0.95,
                                   base_size = 12) {

  # extract individual estimates
  bind_rows(
    summary(bedcast, "growthrate_reported", alpha = alpha) |>
      mutate(type = "observed"),
    summary(bedcast, "growthrate_projected", alpha = alpha) |>
      mutate(type = "projected")
  ) %>%
    mutate(lower = .[[3]], mid = .[[4]], upper = .[[5]]) |>
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
