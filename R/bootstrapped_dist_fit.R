#' Bootstrap distribution fitting
#'
#' Fits a log-normal distribution to data using bootstrap resampling to account
#' for uncertainty in parameter estimates. This function is used internally
#' for estimating delay distributions from linelist data.
#'
#' @param x Numeric vector of data to fit distribution to.
#' @param max_value Numeric value specifying the maximum value to consider.
#'   Values larger than this are excluded from the analysis.
#' @param n_bootstrap Integer specifying number of bootstrap samples (default: 1000).
#'
#' @return A list containing:
#' \itemize{
#'   \item{meanlog}{Bootstrap samples of log-mean parameter}
#'   \item{sdlog}{Bootstrap samples of log-standard deviation parameter}
#'   \item{meanlog.mean}{Mean of log-mean bootstrap samples}
#'   \item{sdlog.mean}{Mean of log-standard deviation bootstrap samples}
#' }
#'
#' @details This function:
#' \itemize{
#'   \item Removes missing values and values exceeding max_value
#'   \item Performs bootstrap resampling of the data
#'   \item Fits log-normal distribution to each bootstrap sample
#'   \item Returns both individual bootstrap samples and their means
#' }
#'
#' @examples
#' \dontrun{
#' # Fit distribution to delay data
#' delays <- bootstrapped_dist_fit(c(1, 2, 3, 4, 5), max_value = 10)
#' }
#'
#' @importFrom stats rlnorm
#' @export
bootstrapped_dist_fit <- function(x, max_value = 30, n_bootstrap = 1000) {
  # Remove missing values and values exceeding max_value
  x_clean <- x[!is.na(x) & x <= max_value & x > 0]

  if (length(x_clean) < 2) {
    # Return default values if insufficient data
    return(list(
      meanlog = rep(log(5), n_bootstrap),
      sdlog = rep(0.5, n_bootstrap),
      meanlog.mean = log(5),
      sdlog.mean = 0.5
    ))
  }

  # Bootstrap resampling
  bootstrap_samples <- replicate(n_bootstrap, {
    sample(x_clean, size = length(x_clean), replace = TRUE)
  })

  # Fit log-normal distribution to each bootstrap sample
  bootstrap_params <- apply(bootstrap_samples, 2, function(sample_data) {
    if (length(sample_data) < 2) {
      return(c(meanlog = log(5), sdlog = 0.5))
    }
    log_data <- log(sample_data)
    c(meanlog = mean(log_data), sdlog = sd(log_data))
  })

  # Extract parameters
  meanlog_samples <- bootstrap_params[1, ]
  sdlog_samples <- bootstrap_params[2, ]

  # Return results
  list(
    meanlog = meanlog_samples,
    sdlog = sdlog_samples,
    meanlog.mean = mean(meanlog_samples),
    sdlog.mean = mean(sdlog_samples)
  )
}
