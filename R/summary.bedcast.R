#' Summarize Stan model outputs
#'
#' Extracts and summarizes parameter estimates from fitted Stan models with
#' appropriate transformations and formatting for different parameter types.
#'
#' @param x A bedcast object.
#' @param par Character string specifying the parameter name to summarize.
#' @param probs Numeric vector specifying the quantiles to calculate
#'   (default: c(0.25, 0.5, 0.75)).
#'
#' @return A tibble containing summarized parameter estimates with columns:
#' \itemize{
#'   \item{day}{Day index (for time-varying parameters)}
#'   \item{lower}{Lower quantile estimate}
#'   \item{mid}{Median estimate}
#'   \item{upper}{Upper quantile estimate}
#' }
#'
#' @details The function handles different parameter types:
#' \itemize{
#'   \item{Time-varying parameters}{Cases, growth rates, ETU/ISO occupancy, alerts}
#'   \item{Projected parameters}{Future projections beyond observed data}
#'   \item{Static parameters}{Single values for the entire time period}
#' }
#'
#' Automatic transformations are applied based on parameter type:
#' \itemize{
#'   \item{Log-mean parameters}{Exponentiated to natural scale}
#' \item{Log parameters}{Exponentiated (cases, alerts, etc.)}
#' \item{Logit parameters}{Inverse logit transformation (CFR, proportions)}
#' }
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
summary.bedcast <- function(x, par, probs = c(0.25, 0.5, 0.75), ...) {

  add_days <- function(par, days) {
    out <- rstan::summary(x$fit, pars = par, probs = probs) %$%
      as_tibble(summary)
    if (length(probs) == 3) {
      transmute(
        out,
        day = days,
        lower = out[[4]],
        mid = out[[5]],
        upper = out[[6]]
      )
    } else {
      mutate(
        out[4:(4 + length(probs) - 1)],
        day = days
      ) %>%
        select(any_of("day"), everything())
    }
  }

  if (grepl(paste(c("cases", "deaths", "growth"), collapse = "|"), par) &
        !grepl(paste(c("proj", "slope"), collapse = "|"), par)) {
    out <- add_days(par, x$data$day)
  } else if (
    grepl(paste(c("etu", "iso", "alerts"), collapse = "|"), par) &
      !grepl(paste(c("per", "prop", "to"), collapse = "|"), par) &
      !grepl("background", par)) {
    out <- add_days(par, seq_len(x$data$n_obs + x$data$n_proj))
  } else if (grepl("proj", par)) {
    out <- add_days(par, (x$data$n_obs + 1):(x$data$n_obs + x$data$n_proj))
  } else {
    out <- add_days(par, NULL)
  }

  if (grepl("logmean", par)) out %<>% mutate(across(-any_of("day"), exp))
  if (par %in% c("alerts_per_case", "alerts_background",
                 "log_cases_inflated", "log_cases_fitted",
                 "log_cases_missed", "log_cases_projected")) {
    out %<>% mutate(across(-any_of("day"), exp))
  }

  return(out)

}
