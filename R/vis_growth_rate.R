#' Visualize growth rate over time
#'
#' Creates a visualization showing the estimated daily growth rate of cases over
#' time, including uncertainty bounds. This helps identify periods of increasing
#' or decreasing transmission.
#'
#' @param results A list containing the fitted Stan model results from fit_stan().
#' @param base_size Numeric base font size for the plot (default: 12).
#'
#' @return A ggplot object showing the growth rate over time.
#'
#' @details The visualization shows:
#' \itemize{
#'   \item Individual MCMC samples as light lines (showing uncertainty)
#'   \item Median growth rate as a prominent red line
#'   \item Horizontal reference line at zero growth
#'   \item Time series of daily growth rates
#' }
#'
#' Positive growth rates indicate increasing transmission, while negative rates
#' indicate decreasing transmission. The visualization helps identify:
#' \itemize{
#'   \item Periods of epidemic growth or decline
#'   \item Uncertainty in growth rate estimates
#'   \item Changes in transmission dynamics over time
#' }
#'
#' @examples
#' \dontrun{
#' # Create growth rate visualization
#' p <- vis_growth_rate(results)
#'
#' # Customize base size
#' p <- vis_growth_rate(results, base_size = 14)
#'
#' # Save the plot
#' p %>% save_plot("growth_rate.png")
#' }
#'
#' @importFrom dplyr group_by summarise
#' @importFrom ggplot2 ggplot aes geom_line geom_hline scale_y_continuous
#' @importFrom ggplot2 scale_x_date labs theme_minimal theme element_rect
#' @importFrom reshape2 melt
#' @importFrom tibble tibble
#' @export
vis_growth_rate <- function(results, base_size = 12) {
  lines <- extract_stan("growth_rate_vec", results) %>%
    reshape2::melt() %>%
    mutate(date = as.Date(Var2, origin = min(results$data$date) - 1))

  med <- lines %>%
    group_by(date) %>%
    summarise(value = median(value))

  ggplot(mapping = aes(date, value)) +
    geom_line(aes(group = iterations), data = lines, alpha = 0.04) +
    geom_line(data = med, color = "firebrick", size = 1) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_date(
      expand = c(0, 0),
      date_labels = "%b %d",
      date_breaks = "2 weeks"
    ) +
    labs(
      x = NULL,
      y = "Daily growth rate"
    ) +
    theme_minimal(base_size) +
    theme(plot.background = element_rect(fill = "white", color = NA))
}
