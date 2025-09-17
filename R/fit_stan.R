#' Fit Stan model for bed occupancy forecasting
#'
#' Fits a Bayesian model using Stan to forecast ETU and isolation bed occupancy
#' based on case incidence, alerts, and delay distributions. The model accounts
#' for reporting delays and provides nowcasts and forecasts of healthcare demand.
#'
#' @param data A tsibble containing the merged time series data.
#' @param model A compiled Stan model object.
#' @param prior_onset_to_etu Numeric vector of length 2 specifying the prior
#'   for log-mean and log-sd of onset-to-ETU delay (default: c(log(5), 0.25)).
#' @param prior_etu_to_survival Numeric vector of length 2 specifying the prior
#'   for log-mean and log-sd of ETU-to-survival delay (default: c(log(10), 0.25)).
#' @param prior_etu_to_death Numeric vector of length 2 specifying the prior
#'   for log-mean and log-sd of ETU-to-death delay (default: c(log(6), 0.25)).
#' @param prior_onset_to_iso Numeric vector of length 2 specifying the prior
#'   for log-mean and log-sd of onset-to-isolation delay (default: c(log(3), 0.25)).
#' @param prior_iso_to_release Numeric vector of length 2 specifying the prior
#'   for log-mean and log-sd of isolation-to-release delay (default: c(log(3), 0.25)).
#' @param prior_cfr Numeric vector of length 2 specifying the prior for
#'   case fatality rate (logit scale) (default: c(qlogis(0.25), 0.25)).
#' @param prior_prop_iso Numeric vector of length 2 specifying the prior for
#'   proportion of alerts isolated (logit scale) (default: c(qlogis(0.8), 0.8)).
#' @param prior_alerts_background Numeric vector of length 2 specifying the prior
#'   for background alerts (log scale) (default: c(log(50), 0.25)).
#' @param prior_alerts_per_case Numeric vector of length 2 specifying the prior
#'   for alerts per case (log scale) (default: c(log(10), 0.25)).
#' @param days_ahead Integer specifying number of days to forecast ahead
#'   (default: 28).
#' @param n_knots Integer specifying number of knots for the growth rate spline
#'   (default: 15).
#' @param n_iter Integer specifying number of MCMC iterations (default: 100).
#' @param alerts_background_window Integer specifying window size for background
#'   alerts estimation (default: 14).
#' @param n_chains Integer specifying number of MCMC chains (default: 1).
#'
#' @return A list containing:
#' \itemize{
#'   \item{stan_fit}{The fitted Stan model object}
#'   \item{data}{The processed data used for fitting, including original data}
#' }
#'
#' @details The Stan model implements a hierarchical Bayesian approach that:
#' \itemize{
#'   \item Models case incidence using a spline-based growth rate
#'   \item Accounts for reporting delays using delay distributions
#'   \item Links cases to ETU occupancy through admission and discharge processes
#'   \item Models alerts as a function of cases and background noise
#'   \item Provides nowcasts of unobserved cases and forecasts of future demand
#' }
#'
#' The model automatically saves results to the outputs directory with a timestamp.
#'
#' @examples
#' \dontrun{
#' # Compile model
#' model <- stan_model("inst/stan/occupancy_model_nowcast.stan")
#'
#' # Fit with default priors
#' results <- fit_stan(data, model)
#'
#' # Fit with custom priors
#' results <- fit_stan(data, model,
#'   prior_cfr = c(qlogis(0.3), 0.2),
#'   days_ahead = 14
#' )
#' }
#'
#' @importFrom rstan sampling
#' @importFrom splines bs
#' @importFrom here here
#' @importFrom rio export
#' @importFrom glue glue
#' @importFrom dplyr mutate
#' @importFrom stats qlogis
#' @export
fit_stan <- function(data, model,
                     prior_onset_to_etu = c(log(5), 0.25),
                     prior_etu_to_survival = c(log(10), 0.25),
                     prior_etu_to_death = c(log(6), 0.25),
                     prior_onset_to_iso = c(log(3), 0.25),
                     prior_iso_to_release = c(log(3), 0.25),
                     prior_cfr = c(qlogis(0.25), 0.25),
                     prior_prop_iso = c(qlogis(0.8), 0.8),
                     prior_alerts_background = c(log(50), 0.25),
                     prior_alerts_per_case = c(log(10), 0.25),
                     days_ahead = 28,
                     n_knots = 15,
                     n_iter = 100,
                     alerts_background_window = 14,
                     n_chains = 1) {
  ## generate spline matrix with buffer on either side to prevent extremes
  buffer <- 16
  spline <- t(bs(
    x = seq_len(nrow(data) + buffer),
    knots = seq(1, nrow(data), length = n_knots),
    degree = 3,
    intercept = TRUE
  ))

  sq_keep <- seq(buffer / 2 + 1, nrow(data) + buffer / 2)
  spline <- spline[, sq_keep]

  n_alerts_background <- floor(nrow(data) / alerts_background_window)
  alerts_background_ind <- rep(
    seq_len(n_alerts_background),
    each = alerts_background_window
  )
  alerts_background_ind <- c(
    rep(1, nrow(data) - length(alerts_background_ind)),
    alerts_background_ind
  )

  stan_data <- list(
    n_days = nrow(data),
    max_delay = 50,
    days_ahead = days_ahead,
    cases_observed = replace_na(data$cases, -1000),
    etu_observed = replace_na(data$etu, -1000),
    etu_n = sum(!is.na(data$etu)),
    etu_ind = which(!is.na(data$etu)),
    alerts_observed = replace_na(data$alerts, -1000),
    alerts_n = sum(!is.na(data$alerts)),
    alerts_ind = which(!is.na(data$alerts)),
    iso_observed = replace_na(data$iso, -1000),
    iso_n = sum(!is.na(data$iso)),
    iso_ind = which(!is.na(data$iso)),
    prior_onset_to_etu = prior_onset_to_etu,
    prior_etu_to_survival = prior_etu_to_survival,
    prior_etu_to_death = prior_etu_to_death,
    prior_onset_to_iso = prior_onset_to_iso,
    prior_iso_to_release = prior_iso_to_release,
    prior_cfr = prior_cfr,
    prior_prop_iso = prior_prop_iso,
    prior_alerts_background = prior_alerts_background,
    prior_alerts_per_case = prior_alerts_per_case,
    prop_cases_observed = data$prop,
    log_prop_cases_observed = log(data$prop),
    n_alerts_background = n_alerts_background,
    alerts_background_ind = alerts_background_ind,
    n_spline_param = nrow(spline),
    spline = spline
  )

  init_fun <- function(...) {
    list(
      log_cases_missed =
        log((1 + data$cases) / (stan_data$prop_cases_observed) - 1)
    )
  }

  options(mc.cores = parallel::detectCores() - 1)

  stan_fit <<- sampling(
    model,
    data = stan_data,
    chains = n_chains,
    iter = n_iter,
    open_progress = TRUE,
    verbose = TRUE,
    refresh = TRUE,
    init = init_fun
  )

  stan_data$etu_observed <- data$etu
  stan_data$alerts_observed <- data$alerts
  stan_data$iso_observed <- data$iso
  stan_data$day <- seq_len(stan_data$n_days)
  stan_data$date <- data$date

  out <- list(stan_fit = stan_fit, data = stan_data)

  export(
    out,
    here(
      glue("outputs/results_{format(Sys.time(), '%Y-%m-%d_%H-%m')}.rds")
    )
  )

  return(out)
}
