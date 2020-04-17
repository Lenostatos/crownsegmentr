# Set up renv for package management
install.packages("remotes")
remotes::install_github("rstudio/renv")
renv::init()

# Install packages needed for development
install.packages(c("roxygen2", "testthat", "knitr"))
remotes::install_github("r-lib/devtools")
remotes::install_github("r-lib/usethis")
install.packages("styler")

# Install packages that this package depends on
install.packages("BH")
