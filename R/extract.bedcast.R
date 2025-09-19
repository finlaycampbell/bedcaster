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
#' @importFrom rstan extract
#' @importFrom stats plogis
#' @export
extract.bedcast <- function(bedcast, par, output = c("raw", "shaped")) {

  # match arg
  output <- match.arg(output)

  # extract value
  out <- rstan::extract(bedcast$fit, pars = par)[[1]]

  # transform if need be
  if (grepl("logmean", par)) out <- exp(out)
  if (par %in% c("alerts_per_case", "alerts_background")) out <- exp(out)
  if (grepl("cfr", par) | grepl("prop_iso", par)) out <- plogis(out)

  if (output == "shaped") {
    if (length(dim(out)) == 2) {
      out <- data.frame(
        iter = rep(seq_len(nrow(out)), times = ncol(out)),
        index = rep(seq_len(ncol(out)), each = nrow(out)),
        value = c(out)
      )
    } else {
      out <- data.frame(
        iter = seq_along(out),
        value = out
      )
    }
  }

  return(out)

}
