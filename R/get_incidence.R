#' Import incidence data from scraped sources
#'
#' Imports daily case incidence data from CSV files that have been scraped
#' from external sources. This function handles the conversion of date formats
#' and case counts.
#'
#' @param file Character string specifying the filename of the incidence data.
#'   Default is "incidence_scraped.csv".
#'
#' @return A tibble containing the processed incidence data with columns:
#' \itemize{
#'   \item{date}{Date of the cases}
#'   \item{cases}{Number of cases on that date (rounded to integers)}
#' }
#'
#' @details The function performs the following operations:
#' \itemize{
#'   \item Converts date column from DD/MM/YYYY format to proper Date format
#'   \item Rounds case counts to integers
#'   \item Renames columns to standard format
#' }
#'
#' @examples
#' \dontrun{
#' # Import default incidence file
#' incidence <- get_incidence()
#'
#' # Import specific file
#' incidence <- get_incidence("custom_incidence.csv")
#' }
#'
#' @importFrom rio import
#' @importFrom here here
#' @importFrom dplyr transmute
#' @importFrom tibble as_tibble
#' @export
get_incidence <- function(file = "incidence_scraped.csv") {
  import(here("data", file)) %>%
    transmute(date = as.Date(x, format = "%d/%m/%Y"), cases = round(y))
}
