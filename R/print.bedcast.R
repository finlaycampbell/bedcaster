#' Print method for bedcaster objects
#'
#' @param x A bedcaster object
#' @param ... Additional arguments (not used)
#' @export
print.bedcast <- function(x, ...) {

  cat("\n Bedcast object\n\n")

  ## --- 1. Data summary ---
  n_days       <- length(x$data$cases_reported)
  total_cases  <- sum(x$data$cases_reported,  na.rm = TRUE)
  total_deaths <- sum(x$data$deaths_reported, na.rm = TRUE)
  total_etu    <- sum(x$data$etu_reported,    na.rm = TRUE)
  total_alerts <- sum(x$data$alerts_reported, na.rm = TRUE)
  total_iso    <- sum(x$data$iso_reported,    na.rm = TRUE)

  data_lines <- c(
    sprintf(" Days:      %d", n_days),
    sprintf(" Cases:     %d", total_cases),
    sprintf(" Deaths:    %d", total_deaths),
    sprintf(" ETU:       %d", total_etu),
    sprintf(" Alerts:    %d", total_alerts),
    sprintf(" Isolation: %d", total_iso)
  )

  ## --- 2. Stan fit summary ---
  fit <- x$fit
  if (!inherits(fit, "stanfit")) {
    fit_lines <- "No stanfit object available."
  } else {
    sim       <- fit@sim
    n_chain   <- sim$chains
    iter      <- sim$iter
    warmup    <- sim$warmup
    thin      <- sim$thin
    post_iter <- iter - warmup
    draws_pc  <- post_iter / thin
    total     <- draws_pc * n_chain

    fit_lines <- c(
      sprintf("Chains:    %d", n_chain),
      sprintf("Iter:      %d", iter),
      sprintf("Warmup:    %d", warmup),
      sprintf("Thin:      %d", thin),
      sprintf("Samples:   %d", total)
    )
  }

  ## --- Align side by side ---
  n <- max(length(data_lines), length(fit_lines))
  data_lines <- c(data_lines, rep("", n - length(data_lines)))
  fit_lines  <- c(fit_lines,  rep("", n - length(fit_lines)))

  cat(sprintf("%-20s | %s\n", " Data", "Fit"))
  cat(" ", strrep("-", 38), "\n", sep = "")
  for (i in seq_len(n)) {
    cat(sprintf("%-20s | %s\n", data_lines[i], fit_lines[i]))
  }
  cat("\n")

  invisible(x)

}
