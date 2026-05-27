#' Package setup and configuration
#'
#' These directives ensure the NAMESPACE is generated correctly
#' for use with Rcpp, rstan, and rstantools.
#'
#' @import Rcpp
#' @import methods
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib bedcaster, .registration = TRUE

#' @keywords internal
#' @noRd
.onLoad <- function(libname, pkgname) {
  # Load Stan modules when package loads
  Rcpp::loadModule("stan_fit4bedcaster_mod", what = TRUE)
}

# NSE columns used in dplyr and ggplot pipelines
utils::globalVariables(c(
  ".", "Rhat", "density", "group", "index", "inflated", "lower", "maxval",
  "mid", "minval", "n_eff", "q_0.5", "reported", "se_mean", "type", "upper",
  "value", "var", "what"
))
