#' Printing generic for class bedcast
#'
#' @param bedcast Object of class \code{bedcast} returned by
#'   fit_bedcaster.
#' @param type What type of visualisation. "fit" will call
#'   \code{vis_bedcast_fit}, "growthrate" will call
#'   \code{vis_bedcast_growthrate} and "vis_bedcast_parameters" will call
#'   \code{vis_parameters}.
#' @param ... Other arguments passed to plotting functions.
#'
#' @export
#'
print.bedcast <- function(x, ...) {

  cat("Bedcast object\n")
  cat("==============\n")

  invisible(x)

}
