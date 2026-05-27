#' Calculate area under curve for parameter transformations
#'
#' Calculates the area under the curve for prior distributions after applying
#' transformations. This is used to properly normalize prior densities when
#' visualizing parameter distributions.
#'
#' @param prior Numeric vector of length 2 containing the mean and standard
#'   deviation of the prior distribution (on the transformed scale).
#' @param transform Character string specifying the transformation to apply.
#'   Options are "plogis" (logistic), "exp" (exponential), or "identity" (no transformation).
#'
#' @return A numeric value representing the area under the curve after transformation.
#'
#' @details This function is used internally by \code{vis_bedcast_parameters()}
#' to properly
#' normalize prior distributions when they are transformed to different scales.
#' The calculation involves:
#' \itemize{
#'   \item Generating a fine grid of values covering the prior range
#'   \item Applying the specified transformation
#'   \item Calculating the area using numerical integration
#' }
#'
#' @examples
#' \dontrun{
#' # Calculate AUC for logistic transformation
#' auc <- get_auc(c(0, 1), "plogis")
#'
#' # Calculate AUC for exponential transformation
#' auc <- get_auc(c(log(5), 0.5), "exp")
#'
#' # Calculate AUC for identity transformation
#' auc <- get_auc(c(0, 1), "identity")
#' }
#'
#' @importFrom stats qnorm dnorm plogis
#' @importFrom zoo rollmean
#' @export
get_auc <- function(prior, transform = c("plogis", "exp", "identity")) {
  full_range <- qnorm(c(0.01, 0.99), prior[1], prior[2])
  full_sq <- seq(full_range[1], full_range[2], length = 100)
  full_density <- dnorm(full_sq, prior[1], prior[2])
  trans_sq <- if (transform == "plogis") {
    plogis(full_sq)
  } else if (transform == "exp") {
    exp(full_sq)
  } else if (transform == "identity") full_sq

  id <- order(full_sq)
  sum(diff(trans_sq[id]) * zoo::rollmean(full_density[id], 2))
}
