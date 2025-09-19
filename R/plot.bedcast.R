#' Plotting generic for class bedcast
#'
#' @param x Object of class \code{bedcast} returned by
#'   fit_bedcaster.
#' @param type What type of visualisation. "fit" will call
#'   \code{vis_bedcast_fit}, "growthrate" will call
#'   \code{vis_bedcast_growthrate} and "vis_bedcast_parameters" will call
#'   \code{vis_parameters}.
#' @param ... Other arguments passed to plotting functions.
#'
#' @export
#'
plot.bedcast <- function(x,
                         type = c("fit", "growthrate", "parameters"),
                         ...) {

  # match arg
  type <- match.arg(type)

  # get plotting function
  f <- list(
    fit = vis_bedcast_fit,
    growthrate = vis_bedcast_growthrate,
    parameters = vis_bedcast_parameters
  )[[type]]

  # execute
  f(x, ...)

}
