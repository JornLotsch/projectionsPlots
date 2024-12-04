# utils.R

#' Load Required Packages
#'
#' Load necessary R packages.
#'
#' @param packages A character vector of package names.
#' @return None
#' @examples
#' load_required_packages(c("ggplot2", "dplyr"))
load_required_packages <- function(packages) {
  missing_packages <- packages[!packages %in% installed.packages()[, "Package"]]
  if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages, collapse = ", "))
  invisible(lapply(packages, library, character.only = TRUE))
}