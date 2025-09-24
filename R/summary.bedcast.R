#' Summarize Stan model outputs
#'
#' Extracts and summarizes parameter estimates from fitted Stan models with
#' appropriate transformations and formatting for different parameter types.
#'
#' @param x A bedcast object.
#' @param par Character string specifying the parameter name to summarize.
#' @param alpha The coverage of the credible interval.
#' @param probs Numeric vector specifying the quantiles to calculate
#'   (default: c(0.25, 0.5, 0.75)), overrides specification of alpha.
#'
#' @examples
#' \dontrun{
#' # Summarize case estimates
#' cases <- summary(bedcast, "log_cases_inflated")
#'
#' # Get specific quantiles
#' cases <- summary(bedcast, "log_cases_inflated",
#'   probs = c(0.025, 0.5, 0.975)
#' )
#'
#' # Summarize static parameters
#' cfr <- summary(bedcast, "cfr")
#' }
#'
#' @importFrom dplyr transmute mutate select
#' @importFrom tibble as_tibble
#' @importFrom stats plogis
#' @export
#'
summary.bedcast <- function(x, par, alpha = 0.95, probs) {

  if (missing(probs))
    probs <- c(0.5 - alpha / 2, 0.5, 0.5 + alpha / 2)

  out <- rstan::summary(x$fit, pars = par, probs = probs)$summary
  colnames(out)[grepl("%", colnames(out))] <- paste0(
    "q_", format(probs, scientific = FALSE)
  )

  # match variable names to their indices
  mtch <- list(
    date_fit = c("reported", "truncated", "nowcast"),
    date_projection = c("projected"),
    date_total = c("etu", "alerts", "iso")
  )

  index <- names(mtch)[
    map_lgl(mtch, ~ grepl(paste(.x, collapse = "|"), par))
  ]

  # etu/iso/alerts match date_projection and date_total - use total
  index <- if (length(index) == 0) seq_len(nrow(out))
  else x$data[[tail(index, 1)]]

  mutate(as_tibble(out), index = index, .before = mean) |>
    select(-c(se_mean, sd, n_eff, Rhat))

}
