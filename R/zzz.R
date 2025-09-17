# Package setup and configuration
# This file is loaded when the package is attached

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Loading bedcast package for bed occupancy forecasting")
}

.onLoad <- function(libname, pkgname) {
  # Set up any package-specific options here if needed
  invisible()
}
