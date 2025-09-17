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

.onLoad <- function(libname, pkgname) {
  # Load Stan modules when package loads
  Rcpp::loadModule("stan_fit4bedcaster_mod", what = TRUE)
}
