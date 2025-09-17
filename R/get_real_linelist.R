#' Import real/raw linelist data
#'
#' Imports and processes real linelist data from Excel files, cleaning column names
#' and standardizing date formats. This function handles the conversion of various
#' date formats commonly found in epidemiological data.
#'
#' @param file Character string specifying the filename of the linelist data.
#'   Default is "EVD linelist_10Nov_0000hrs.xlsx".
#'
#' @return A tibble containing the processed linelist data with standardized
#'   column names and date formats. Key columns include:
#'   \itemize{
#'     \item{subcounty}{Subcounty location}
#'     \item{district}{District location}
#'     \item{case_class}{Classification of the case}
#'     \item{date}{Standardized date (coalesced from onset and sample collection dates)}
#'     \item{date_of_contact_exposure}{Date of contact exposure}
#'   }
#'
#' @details The function performs several data cleaning operations:
#' \itemize{
#'   \item Converts column names to snake_case using janitor::clean_names()
#'   \item Handles both POSIXct and Excel numeric date formats
#'   \item Creates a standardized 'date' column by coalescing onset and sample collection dates
#'   \item Removes N/A values and converts them to proper NA values
#' }
#'
#' @examples
#' \dontrun{
#' # Import default linelist file
#' ll <- get_real_linelist()
#'
#' # Import specific file
#' ll <- get_real_linelist("custom_linelist.xlsx")
#' }
#'
#' @importFrom rio import
#' @importFrom here here
#' @importFrom janitor clean_names
#' @importFrom dplyr select mutate across coalesce
#' @importFrom tibble as_tibble
#' @export
get_real_linelist <- function(file = "EVD linelist_10Nov_0000hrs.xlsx") {
  import(here("data", file)) %>%
    as_tibble() %>%
    clean_names() %>%
    select(subcounty, district, case_class, contains("date")) %>%
    mutate(
      across(
        contains("date"),
        ~ if (inherits(.x, "POSIXct")) {
          as.Date(.x)
        } else {
          as.Date(as.numeric(na_if(.x, "N/A")), origin = "1900-01-01") - 2
        }
      ),
      date = coalesce(date_of_onset, date_of_sample_collection),
      .before = date_of_contact_exposure
    )
}
