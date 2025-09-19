#' Fit Stan model for bed occupancy forecasting
#'
#' Fits a Bayesian model using Stan to forecast ETU and isolation bed occupancy
#' based on case incidence, alerts, and delay distributions. The model accounts
#' for reporting delays and provides nowcasts and forecasts of healthcare demand.
#'
#' @param data A tibble containing the merged time series data.
#' @param data_as_of Date of dataset used for calculating reporting delay.
#' @param prior_onset_to_reporting Numeric vector of length 2 specifying the prior
#'   for log-mean and log-sd of onset-to-reporting delay (default: c(log(5), 0.25)).
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
#' @param max_delay The maximum number of days for a delay.
#' @param n_proj Integer specifying number of days to project ahead
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
#'   \item Provides nowcasts of unreported cases and forecasts of future demand
#' }
#'
#' The model automatically saves results to the outputs directory with a timestamp.
#'
#' @examples
#' \dontrun{
#'
#' # Fit with default priors
#' results <- fit_stan(data)
#'
#' # Fit with custom priors
#' results <- fit_stan(data,
#'   prior_cfr = c(qlogis(0.3), 0.2),
#'   n_proj = 14
#' )
#' }
#'
#' @importFrom rstan sampling
#' @importFrom splines bs
#' @importFrom stats qlogis
#' @importFrom distcrete distcrete
#' @export
#'
fit_bedcaster <- function(data, data_as_of,
                          prior_onset_to_reporting = c(0, 5, 0, 5),
                          prior_onset_to_etu = c(0, 5, 0, 5),
                          prior_etu_to_survival = c(0, 5, 0, 5),
                          prior_etu_to_death = c(0, 5, 0, 5),
                          prior_onset_to_iso = c(0, 5, 0, 5),
                          prior_iso_to_release = c(0, 5, 0, 5),
                          prior_cfr = c(qlogis(0.25), 0.25),
                          prior_prop_iso = c(qlogis(0.8), 0.8),
                          prior_alerts_background = c(log(50), 0.25),
                          prior_alerts_per_case = c(log(10), 0.25),
                          growthrate_asymptote_time = 20,
                          growthrate_asymptote_spread = 3,
                          n_proj = 28,
                          n_knots = 15,
                          n_iter = 100,
                          alerts_background_window = 28,
                          max_delay = 50,
                          n_chains = 1,
                          n_cores = 1) {

  # specify discrete lognormal
  onset_to_reporting <- distcrete(
    name = "lnorm", interval = 1,
    meanlog = prior_onset_to_reporting["meanlog.mean"],
    sdlog = prior_onset_to_reporting["sdlog.mean"]
  )

  # define probability of reporting for given delays
  prop_cases_reported <- onset_to_reporting$p(
    as.numeric(data_as_of - data$date)
  )

  ## generate spline matrix with buffer on either side to prevent extremes
  buffer <- 16
  spline <- t(bs(
    x = seq_len(nrow(data) + buffer),
    knots = seq(1, nrow(data), length = n_knots),
    degree = 3,
    intercept = TRUE
  ))

  # define the weight of the growth rate asymptote
  growthrate_asymptote_weight <- plogis(
    seq_len(n_proj),
    growthrate_asymptote_time,
    growthrate_asymptote_spread
  )

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
    n_obs = nrow(data),
    n_proj = n_proj,
    max_delay = max_delay,
    cases_reported = replace_na(data$cases, -1000),
    deaths_reported = replace_na(data$deaths, -1000),
    etu_reported = replace_na(data$etu, -1000),
    etu_n = sum(!is.na(data$etu)),
    etu_ind = which(!is.na(data$etu)),
    alerts_reported = replace_na(data$alerts, -1000),
    alerts_n = sum(!is.na(data$alerts)),
    alerts_ind = which(!is.na(data$alerts)),
    iso_reported = replace_na(data$iso, -1000),
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
    log_prop_cases_reported = log(prop_cases_reported),
    prop_deaths_reported = prop_cases_reported,
    n_alerts_background = n_alerts_background,
    alerts_background_ind = alerts_background_ind,
    n_spline_param = nrow(spline),
    spline = spline,
    growthrate_asymptote_weight = growthrate_asymptote_weight
  )

  init_fun <- function(...) {
    list(
      log_cases_missed =
        log((1 + data$cases) / (stan_data$prop_cases_reported) - 1)
    )
  }

  options(mc.cores = n_cores)

  fit <- rstan::sampling(
    stanmodels$bedcaster,
    data = stan_data,
    chains = n_chains,
    iter = n_iter,
    verbose = TRUE,
    refresh = 10,
    init = init_fun
  )

  stan_data$etu_reported <- data$etu
  stan_data$alerts_reported <- data$alerts
  stan_data$iso_reported <- data$iso
  stan_data$day <- seq_len(stan_data$n_obs)
  stan_data$date <- data$date

  out <- list(fit = fit, data = stan_data)
  class(out) <- "bedcast"

  return(out)

}
