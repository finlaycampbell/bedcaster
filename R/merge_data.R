#' Merge different data sources into a single dataframe
#'
#' Combines linelist, alerts, ETU/isolation occupancy, delays, and incidence data
#' into a single time series dataframe suitable for modeling. This function
#' handles missing data and creates a complete time series with proper indexing.
#'
#' @param ll A tibble containing linelist data (used to calculate incidence if not provided).
#' @param alerts A tibble containing daily alerts data.
#' @param etu_iso A tibble containing ETU and isolation occupancy data.
#' @param delays A list containing delay distribution parameters.
#' @param incidence Optional tibble containing pre-calculated incidence data.
#'   If NULL, incidence will be calculated from the linelist data.
#'
#' @return A tibble containing merged data with columns:
#' \itemize{
#'   \item{date}{Date index}
#'   \item{day}{Numeric day index (1, 2, 3, ...)}
#'   \item{cases}{Daily number of cases}
#'   \item{etu}{Daily ETU occupancy}
#'   \item{iso}{Daily isolation occupancy}
#'   \item{alerts}{Daily number of alerts}
#'   \item{prop}{Proportion of cases observed (for case inflation modeling)}
#' }
#'
#' @details The function performs several key operations:
#' \itemize{
#'   \item Calculates daily incidence from linelist if not provided
#'   \item Aggregates alerts data by date
#'   \item Fills gaps in the time series with appropriate values
#'   \item Creates a case inflation factor based on onset-to-confirmation delays
#'   \item Converts to tibble format for analysis
#' }
#'
#' The case inflation factor (prop) is calculated using a log-normal distribution
#' based on the delay from onset to confirmation, which accounts for reporting delays.
#'
#' @examples
#' \dontrun{
#' # Merge all data sources
#' data <- merge_data(ll, alerts, etu_iso, delays)
#'
#' # Merge with pre-calculated incidence
#' data <- merge_data(ll, alerts, etu_iso, delays, incidence)
#' }
#'
#' @importFrom dplyr full_join arrange group_by summarise mutate
#' @importFrom stats setNames
#' @importFrom tidyr replace_na complete
#' @importFrom tibble as_tibble
#' @importFrom distcrete distcrete
#' @importFrom magrittr %>% %$%
#' @export
merge_data <- function(ll, alerts, etu_iso, delays, incidence = NULL) {
  if (is.null(incidence)) {
    incidence <- table(ll$date) %>%
      as.data.frame() %>%
      setNames(c("date", "cases")) %>%
      mutate(date = as.Date(date))
  }

  incidence %>%
    full_join(etu_iso, by = "date") %>%
    full_join(
      summarise(group_by(alerts, date), alerts = sum(alerts)),
      by = "date"
    ) %>%
    arrange(date) %>%
    complete(
      date = seq.Date(min(date), max(date), by = 1),
      fill = list(cases = 0, etu = NA, iso = NA, alerts = NA)
    ) %>%
    as_tibble() %>%
    mutate(
      day = as.numeric(date - min(date) + 1),
      ## fit distribution from onset to confirmation for case inflation
      prop = distcrete(
        name = "lnorm", interval = 1,
        meanlog = delays$onset_to_confirmation["meanlog.mean"],
        sdlog = delays$onset_to_confirmation["sdlog.mean"]
      ) %$%
        d(seq_along(day) - 1) %>% cumsum() %>% rev(),
      .before = cases
    )
}
