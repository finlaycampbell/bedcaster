#' Calculate delay distributions from linelist data
#'
#' Estimates delay distributions between key epidemiological events using
#' bootstrap resampling. This function calculates the time intervals between
#' different stages of disease progression and healthcare utilization.
#'
#' @param ll A tibble containing linelist data with date columns for different
#'   epidemiological events.
#' @param max_value Numeric value specifying the maximum delay to consider
#'   (default: 30 days). Delays longer than this are excluded from the analysis.
#'
#' @return A list containing delay distribution parameters for each delay type:
#' \itemize{
#'   \item{onset_to_etu}{Delay from symptom onset to ETU admission}
#'   \item{etu_to_survival}{Delay from ETU admission to discharge (survival)}
#'   \item{etu_to_death}{Delay from ETU admission to death}
#'   \item{onset_to_iso}{Delay from symptom onset to contact isolation}
#'   \item{iso_to_release}{Lab turnaround time (isolation to confirmation)}
#'   \item{onset_to_confirmation}{Delay from onset to case confirmation}
#' }
#'
#' Each delay type contains:
#' \itemize{
#'   \item{meanlog}{Log-mean parameter for log-normal distribution}
#'   \item{sdlog}{Log-standard deviation parameter for log-normal distribution}
#' }
#'
#' @details The function uses bootstrap resampling to estimate delay distributions
#' and fits log-normal distributions to the observed delays. This approach
#' accounts for uncertainty in delay estimates and provides robust parameter
#' estimates for downstream modeling.
#'
#' @examples
#' \dontrun{
#' # Calculate delays from linelist data
#' delays <- get_delays(ll)
#'
#' # Use custom maximum delay value
#' delays <- get_delays(ll, max_value = 45)
#' }
#'
#' @importFrom purrr map set_names pluck
#' @export
get_delays <- function(ll, max_value = 30) {
  ll %$%
    map(
      list(
        onset_to_etu = date_of_confirmation - date_of_onset,
        etu_to_survival = date_of_discharge - date_of_isolation,
        etu_to_death = date_of_death - date_of_isolation,
        onset_to_iso = date_of_isolation - date_of_onset,
        iso_to_release = modify_if(
          as.numeric(date_of_confirmation - date_of_isolation),
          ~ !is.na(.x) & .x < 0, ~0
        ),
        onset_to_confirmation = date_of_confirmation - date_of_onset + 1
      ),
      ~ bootstrapped_dist_fit(as.numeric(.x), max_value = max_value)
    ) |>
    map(
      \(delay) unlist(map(
        purrr::set_names(c("meanlog", "sdlog")),
        ~ unlist(pluck(delay, 1, .x, 1))
      ))
    )
}
