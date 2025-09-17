#' Import ETU and isolation occupancy data
#'
#' Imports and processes daily ETU (Ebola Treatment Unit) and isolation
#' occupancy data from Excel files. This data represents the number of beds
#' occupied in treatment and isolation facilities.
#'
#' @param file Character string specifying the filename of the occupancy data.
#'   Default is "sitrep_scrape.xlsx".
#'
#' @return A tibble containing the processed occupancy data with columns:
#' \itemize{
#'   \item{date}{Date of the occupancy data}
#'   \item{etu}{Number of ETU beds occupied}
#'   \item{iso}{Number of isolation beds occupied}
#' }
#'
#' @details The function performs the following operations:
#' \itemize{
#'   \item Converts date column to proper Date format
#'   \item Removes the 'alerts' column (if present) to avoid duplication
#'   \item Removes rows with missing values
#' }
#'
#' @examples
#' \dontrun{
#' # Import default occupancy file
#' etu_iso <- get_etu_iso()
#'
#' # Import specific file
#' etu_iso <- get_etu_iso("custom_occupancy.xlsx")
#' }
#'
#' @importFrom rio import
#' @importFrom here here
#' @importFrom dplyr mutate select
#' @importFrom tidyr drop_na
#' @importFrom tibble as_tibble
#' @export
get_etu_iso <- function(file = "sitrep_scrape.xlsx") {
  import(here("data", file)) %>%
    mutate(date = as.Date(date)) %>%
    select(-alerts) %>%
    drop_na()
}
