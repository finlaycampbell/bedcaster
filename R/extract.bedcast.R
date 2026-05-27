#' Extract raw samples from Stan model
#'
#' Extracts raw MCMC samples from fitted Stan models with appropriate
#' transformations applied based on parameter type.
#'
#' @param x A \code{bedcast} object from \code{\link{fit_bedcaster}}.
#' @param par Character string naming a Stan parameter to extract.
#'
#' @return A tibble with columns \code{index}, \code{iter}, and \code{value}.
#'
#' @examples
#' \dontrun{
#' bedcast <- fit_bedcaster(df, as_of = max(df$date))
#' cfr_samples <- extract(bedcast, "cfr")
#' }
#'
#' @importFrom purrr map_lgl
#' @importFrom tibble tibble
#' @export
#' @method extract bedcast
#'
extract.bedcast <- function(x, par) {

  # extract value
  out <- rstan::extract(x$fit, pars = par)[[1]]

  if (length(dim(out)) == 2) {

    # match variable names to their indices
    mtch <- list(
      date_fit = c("reported", "truncated", "nowcast"),
      date_proj = c("projected")
    )

    index <- names(mtch)[
      map_lgl(mtch, ~ grepl(paste(.x, collapse = "|"), par))
    ]

    index <- if (length(index) == 0) seq_len(ncol(out)) else x$data[[index]]

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
