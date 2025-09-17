#' Save ggplot objects to files
#'
#' A convenience function for saving ggplot objects to various file formats
#' with consistent settings and automatic directory creation.
#'
#' @param p A ggplot object to save. If NULL, the function returns NULL.
#' @param file Character string specifying the output filename.
#' @param width Numeric width of the plot in the specified units (default: 20.16).
#' @param height Numeric height of the plot in the specified units (default: 13.35).
#' @param units Character string specifying the units for width and height
#'   (default: 'cm').
#' @param dpi Numeric resolution in dots per inch (default: 600).
#' @param folder Character string specifying the output folder (default: "figures").
#' @param ... Additional arguments passed to ggsave().
#'
#' @return The input ggplot object (invisibly) for use in pipes.
#'
#' @details This function:
#' \itemize{
#'   \item Creates the output directory if it doesn't exist
#'   \item Uses high-resolution settings by default (600 DPI)
#'   \item Provides sensible default dimensions
#'   \item Returns the plot object for use in pipes
#' }
#'
#' @examples
#' \dontrun{
#' # Save with default settings
#' p %>% save_plot("my_plot.png")
#'
#' # Save with custom dimensions
#' p %>% save_plot("my_plot.png", width = 30, height = 20)
#'
#' # Save to custom folder
#' p %>% save_plot("my_plot.png", folder = "outputs")
#' }
#'
#' @importFrom here here
#' @importFrom ggplot2 ggsave
#' @export
save_plot <- function(p, file,
                      width = 20.16, height = 13.35,
                      units = "cm",
                      dpi = 600,
                      folder = "figures",
                      ...) {
  if (is.null(p)) {
    return(NULL)
  }

  if (!file.exists(here(folder))) dir.create(here(folder))

  ggsave(
    filename = here(folder, file),
    plot = p,
    width = width,
    height = height,
    units = units,
    dpi = dpi,
    ...
  )
}
