# https://stackoverflow.com/a/38928678

load_libs <- function(bundle = "flipFISH") {
  if (length(bundle) == 1L && bundle == "flipFISH") {
    packages <- c(
      "Rcpp",
      "RcppEigen",
      "ggplot2",
      "roxygen2"
    )
  }

  installed_check <- match(packages, utils::installed.packages()[, 1])

  to_be_installed <- packages[is.na(installed_check)]

  if (length(to_be_installed) > 0L) {
    utils::install.packages(to_be_installed,
      repos = "http://cran.wustl.edu"
    )
  } else {
    print("All required packages already installed")
  }

  for (package in packages) {
    suppressPackageStartupMessages(
      library(package, character.only = TRUE, quietly = TRUE)
    )
  }

}

load_libs("flipFISH")