#' Import alerts data
#'
#' Imports and processes daily alerts data from Excel files. This function
#' handles the transformation of wide-format data (with districts as columns)
#' into long-format data suitable for analysis.
#'
#' @param file Character string specifying the filename of the alerts data.
#'   Default is "alerts_data.xlsx".
#'
#' @return A tibble containing the processed alerts data with columns:
#' \itemize{
#'   \item{date}{Date of the alerts}
#'   \item{district}{District name}
#'   \item{alerts}{Number of alerts for that district on that date}
#' }
#'
#' @details The function performs the following operations:
#' \itemize{
#'   \item Pivots data from wide format (districts as columns) to long format
#'   \item Converts date column to proper Date format
#'   \item Removes rows with missing values
#'   \item Excludes "TOTAL" rows to avoid double-counting
#' }
#'
#' @examples
#' \dontrun{
#' # Import default alerts file
#' alerts <- get_alerts()
#'
#' # Import specific file
#' alerts <- get_alerts("custom_alerts.xlsx")
#' }
#'
#' @importFrom rio import
#' @importFrom here here
#' @importFrom dplyr pivot_longer transmute drop_na filter
#' @importFrom tibble as_tibble
#' @export
get_alerts <- function(file = "alerts_data.xlsx") {
  import(here("data", file)) %>%
    as_tibble() %>%
    pivot_longer(-Date, names_to = "district", values_to = "alerts") %>%
    transmute(date = as.Date(Date), district, alerts) %>%
    drop_na() %>%
    filter(district != "TOTAL")
}
