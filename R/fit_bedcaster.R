#' Fit Stan model for bed occupancy forecasting
#'
#' Fits a Bayesian model using Stan to forecast ETU and isolation bed occupancy
#' based on case incidence, alerts, and delay distributions. The model accounts
#' for reporting delays and provides nowcasts and forecasts of healthcare demand.
#'
#' @param df A data frame with columns `date`, `cases`, `deaths`, `etu`,
#'   `alerts`, and `iso` (counts may be `NA` when not observed; `NA` deaths
#'   exclude those days from the death likelihood).
#' @param as_of Date used to compute reporting-delay weights (typically the
#'   latest date in the surveillance system).
#' @param prior_onset_to_reporting Named numeric vector passed to
#'   [distcrete::distcrete()] for the onset-to-reporting delay (`meanlog`,
#'   `sdlog` on the log scale).
#' @param prior_onset_to_etu Named length-4 vector giving Normal priors on
#'   onset-to-ETU log-mean and log-sd (`meanlog.mean`, `meanlog.sd`,
#'   `sdlog.mean`, `sdlog.sd`).
#' @param prior_etu_to_survival Same format as `prior_onset_to_etu` for
#'   hospitalisation-to-survival delay.
#' @param prior_etu_to_death Same format as `prior_onset_to_etu` for
#'   hospitalisation-to-death delay.
#' @param prior_onset_to_iso Same format as `prior_onset_to_etu` for
#'   onset-to-isolation delay.
#' @param prior_iso_to_release Same format as `prior_onset_to_etu` for
#'   isolation-to-release delay.
#' @param prior_cfr Length-2 vector (mean, sd) for in-hospital case fatality
#'   rate on the logit scale.
#' @param prior_prop_iso Length-2 vector (mean, sd) for the proportion of
#'   alerts isolated, on the logit scale.
#' @param prior_alerts_background Length-2 vector (mean, sd) for background
#'   alerts on the log scale.
#' @param prior_alerts_per_case Length-2 vector (mean, sd) for alerts per
#'   case on the log scale.
#' @param extrapolate_growthrate If \code{FALSE} (default), projected growth
#'   rates are held at the last estimated value over the fit period. If
#'   \code{TRUE}, growth rates are extrapolated using recent slope and an
#'   asymptote (see \code{growthrate_asymptote_time}).
#' @param growthrate_asymptote_time Day of projection at which the growth-rate
#'   asymptote reaches half weight (only used when \code{extrapolate_growthrate
#'   = TRUE}).
#' @param growthrate_asymptote_spread Spread of the logistic weighting for the
#'   growth-rate asymptote (only used when \code{extrapolate_growthrate =
#'   TRUE}).
#' @param max_delay The maximum number of days for a delay.
#' @param n_proj Integer specifying number of days to project ahead
#'   (default: 28).
#' @param n_knots Integer specifying number of knots for the growth rate spline
#'   (default: 15).
#' @param n_iter Integer specifying number of MCMC iterations (default: 100).
#' @param alerts_background_window Integer specifying window size for background
#'   alerts estimation (default: 14).
#' @param n_chains Integer specifying number of MCMC chains (default: 1).
#' @param n_cores Number of CPU cores for parallel sampling (default: 1).
#'
#' @return An object of class `"bedcast"` with elements:
#'   \describe{
#'     \item{`fit`}{A `stanfit` object from [rstan::sampling()].}
#'     \item{`data`}{Stan input data plus dates and reporting series used
#'       in plots.}
#'   }
#'
#' @importFrom splines bs
#' @importFrom stats qlogis plogis
#' @importFrom distcrete distcrete
#' @importFrom tidyr replace_na
#' @export
#'
fit_bedcaster <- function(df, as_of,
                          prior_onset_to_reporting = c(
                            meanlog.mean = 0, sdlog.mean = 5
                          ),
                          prior_onset_to_etu = c(
                            meanlog.mean = 0, meanlog.sd = 5,
                            sdlog.mean = 0, sdlog.sd = 5
                          ),
                          prior_etu_to_survival = c(
                            meanlog.mean = 0, meanlog.sd = 5,
                            sdlog.mean = 0, sdlog.sd = 5
                          ),
                          prior_etu_to_death = c(
                            meanlog.mean = 0, meanlog.sd = 5,
                            sdlog.mean = 0, sdlog.sd = 5
                          ),
                          prior_onset_to_iso = c(
                            meanlog.mean = 0, meanlog.sd = 5,
                            sdlog.mean = 0, sdlog.sd = 5
                          ),
                          prior_iso_to_release = c(
                            meanlog.mean = 0, meanlog.sd = 5,
                            sdlog.mean = 0, sdlog.sd = 5
                          ),
                          prior_cfr = c(qlogis(0.25), 0.25),
                          prior_prop_iso = c(qlogis(0.8), 0.8),
                          prior_alerts_background = c(log(50), 0.25),
                          prior_alerts_per_case = c(log(10), 0.25),
                          extrapolate_growthrate = FALSE,
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
    as.numeric(as_of - df$date)
  )

  ## generate spline matrix with buffer on either side to prevent extremes
  buffer <- 16
  spline <- t(bs(
    x = seq_len(nrow(df) + buffer),
    knots = seq(1, nrow(df), length = n_knots),
    degree = 3,
    intercept = TRUE
  ))

  # weight of the growth-rate asymptote (only used when extrapolating)
  growthrate_asymptote_weight <- if (isTRUE(extrapolate_growthrate)) {
    plogis(
      seq_len(n_proj),
      growthrate_asymptote_time,
      growthrate_asymptote_spread
    )
  } else {
    rep(0, n_proj)
  }

  sq_keep <- seq(buffer / 2 + 1, nrow(df) + buffer / 2)
  spline <- spline[, sq_keep]

  n_alerts_background <- floor(nrow(df) / alerts_background_window)
  alerts_background_ind <- rep(
    seq_len(n_alerts_background),
    each = alerts_background_window
  )
  alerts_background_ind <- c(
    rep(1, nrow(df) - length(alerts_background_ind)),
    alerts_background_ind
  )

  deaths_ind <- as.integer(which(!is.na(df$deaths)))

  data <- list(
    n_obs = nrow(df),
    n_proj = n_proj,
    max_delay = max_delay,
    cases_reported = replace_na(df$cases, -1000),
    deaths_reported = replace_na(df$deaths, -1000),
    deaths_n = length(deaths_ind),
    deaths_ind = deaths_ind,
    etu_reported = replace_na(df$etu, -1000),
    etu_n = sum(!is.na(df$etu)),
    etu_ind = which(!is.na(df$etu)),
    alerts_reported = replace_na(df$alerts, -1000),
    alerts_n = sum(!is.na(df$alerts)),
    alerts_ind = which(!is.na(df$alerts)),
    iso_reported = replace_na(df$iso, -1000),
    iso_n = sum(!is.na(df$iso)),
    iso_ind = which(!is.na(df$iso)),
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
    growthrate_asymptote_weight = growthrate_asymptote_weight,
    extrapolate_growthrate = as.integer(isTRUE(extrapolate_growthrate))
  )

  init_fun <- function(...) {
    list(
      log_cases_missed =
        log((1 + df$cases) / (data$prop_cases_reported) - 1)
    )
  }

  options(mc.cores = n_cores)

  fit <- rstan::sampling(
    ## model,
    stanmodels$bedcaster,
    data = data,
    chains = n_chains,
    iter = n_iter,
    verbose = TRUE,
    refresh = 10,
    init = init_fun
  )

  data <- utils::modifyList(
    data,
    list(
      extrapolate_growthrate = isTRUE(extrapolate_growthrate),
      deaths_reported = df$deaths,
      etu_reported = df$etu,
      alerts_reported = df$alerts,
      iso_reported = df$iso,
      date_fit = df$date,
      date_projection = max(df$date) + seq_len(n_proj),
      date_total = c(df$date, max(df$date) + seq_len(n_proj))
    )
  )

  out <- list(fit = fit, data = data)
  class(out) <- "bedcast"

  return(out)

}
