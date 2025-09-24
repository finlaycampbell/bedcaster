#' Extract raw samples from Stan model
#'
#' Extracts raw MCMC samples from fitted Stan models with appropriate
#' transformations applied based on parameter type.
#'
#' @param par Character string specifying the parameter name to extract.
#' @param results A list containing the fitted Stan model results from fit_stan().
#'
#' @return A numeric vector or matrix containing the raw MCMC samples for the
#'   specified parameter.
#'
#' @details The function applies automatic transformations based on parameter type:
#' \itemize{
#'   \item{Log-mean parameters}{Exponentiated to natural scale}
#'   \item{Log parameters}{Exponentiated (alerts, cases, etc.)}
#'   \item{Logit parameters}{Inverse logit transformation (CFR, proportions)}
#' }
#'
#' This function is useful for:
#' \itemize{
#'   \item Custom analysis of parameter distributions
#'   \item Calculating derived quantities
#'   \item Creating custom visualizations
#'   \item Accessing raw samples for further processing
#' }
#'
#' @examples
#' \dontrun{
#' # Extract case fatality rate samples
#' cfr_samples <- extract_stan("cfr", results)
#'
#' # Extract growth rate samples
#' growth_samples <- extract_stan("growth_rate_vec", results)
#'
#' # Extract alerts per case samples
#' alerts_samples <- extract_stan("alerts_per_case", results)
#' }
#'
#' @importFrom stats plogis
#' @importFrom tibble tibble
#' @export
#' @method extract bedcast
#'
#'
extract.bedcast <- function(bedcast, par) {

  # extract value
  out <- rstan::extract(bedcast$fit, pars = par)[[1]]

  if (length(dim(out)) == 2) {

    # match variable names to their indices
    mtch <- list(
      date_fit = c("reported", "truncated", "nowcast"),
      date_proj = c("projected")
    )

    index <- names(mtch)[
      map_lgl(mtch, ~ grepl(paste(.x, collapse = "|"), par))
    ]

    index <- if (length(index) == 0) seq_len(ncol(out)) else bedcast$data[[index]]

    out <- tibble(
      index = rep(index, each = nrow(out)),
      iter = rep(seq_len(nrow(out)), times = ncol(out)),
      value = c(out)
    )

  } else {

    out <- tibble(index = 1, iter = seq_along(out), value = out)

  }

  return(out)

}


#' Extract information from an object
#'
#' A generic function to extract information from objects.
#'
#' @param x An object to extract from.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The result of extraction, depending on the class of \code{x}.
#' @export
#'
extract <- function(x, ...) {
  UseMethod("extract")
}
