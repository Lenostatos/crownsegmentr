# This file is part of crownsegmentr, an R package for identifying tree crowns
# within 3D point clouds.
#
# Copyright (C) 2025 Leon Steinmeier, Nikolai Knapp, UFZ Leipzig
# Contact: timon.miesner@thuenen.de
#
# crownsegmentr is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# crownsegmentr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with crownsegmentr in a file called "COPYING". If not,
# see <http://www.gnu.org/licenses/>.


# Maintenance -------------------------------------------------------------

# Make usethis functions available on the command line without package prefix.
# usethis::use_usethis()

# Package development workflow summary:
# + devtools::document()
# + devtools::load_all()
# + Run some examples interactively.
# + devtools::test()
# + devtools::check()
# Source: https://r-pkgs.org/code.html#constant-health-checks

# Use X after doing Y:
# + Use devtools::load_all() or Ctrl+Shift+L after changing R code to reinstall
#   package.
# + Use "Install and restart" in the Build tab or Ctrl+Shift+B after changing
#   C/C++ code in order to do perform an installation with source compilation.
# + Use devtools::check() to regularly perform R CMD check
# + Use devtools::document() or Ctrl+Shift+D after changing documentation.
# + Use devtools::document() after adding functions in order to add those
#   functions to the NAMESAPCE file.
# + Use devtools::build_readme() after updating README.Rmd.
# + Use usethis::use_test() when editing an R source file in order to create a
#   test file (or switch to it if it already exists).
# + Use usethis::use_r() when editing a test file in order to switch back to the
#   respective R source file.
# + Use devtools::build(binary = TRUE) after updating the package in order to
#   create a compiled bundle of the package that can be installed without
#   devtools and RTools, but only on the OS that it was compiled on.
# + Use devtools::build() after updating the package in order to create a bundle
#   of the package that can be installed on any OS as long as the required
#   compiler toolchain is installed (RTools on Windows and Mac OS and something
#   equivalent on Linux).
# + Use usethis::use_version() for conveniently updating version numbers.
# + Use to style all code according to the tidyverse
#   style guide.
# + Use devtools::dev_sitrep() from time to time to check for dependency
#   updates.
# + Use usethis::use_tidy_description() to format DESCRIPTION.


# Initial development -----------------------------------------------------

# Set up renv for package management
install.packages("renv")
renv::init()

# Install packages needed for development
install.packages(c("devtools", "roxygen2", "knitr", "styler"))

usethis::use_testthat(3)

# Create bare-bones package structure
usethis::create_tidy_package(".")

# Ignore this file
usethis::use_build_ignore("package_development_code.R")

# Use git
usethis::use_git(message = "")

# Add a license
usethis::use_gpl_license(version = 3)

# Add author metadata
usethis::use_author("Leon", "Steinmeier",
  email = "Leon.Steinmeier@posteo.net",
  comment = c(ORCID = "0000-0001-9040-636X", "Created the package as part of his master thesis."),
  role = c("aut", "cre")
)

usethis::use_author("Nikolai", "Knapp",
  email = "nikolai.knapp@thuenen.de",
  comment = "Initialized and motivated the development of the package and provided an R draft of the AMS3D algorithm.",
  role = "aut"
)

usethis::use_author("Timon", "Miesner",
  email = "timon.miesner@thuenen.de",
  comment = c(ORCID = "0000-0001-5091-7456", "Added optimization tools and further options."),
  role = "aut"
)

# Add citation information
usethis::use_citation()
# Read the ?bibentry help topic on how to fill in this information.

# Add a README
usethis::use_readme_md()

# Use unit tests
usethis::use_testthat()

# Require R version with new enough C++ version support
usethis::use_package("R", type = "Depends", min_version = "4.0.0")
# |> operator available from R 4.1.0 onward

# Use C++ with Rcpp
usethis::use_rcpp()

usethis::use_r(name = "crownsegmentr-package.R")
# The following stackoverflow post provides valuable info on the Makevars file:
# https://stackoverflow.com/questions/43597632/understanding-the-contents-of-the-makevars-file-in-r-macros-variables-r-ma

usethis::use_package("Rcpp", min_version = "1.0.0")
usethis::use_package("Rcpp", min_version = "1.0.0", type = "LinkingTo")

# Use other packages
usethis::use_package("BH", type = "LinkingTo", min_version = "1.75.0-0")

usethis::use_package("progress", type = "LinkingTo")

usethis::use_package("assertthat")
usethis::use_import_from("assertthat", "assert_that")

usethis::use_package("methods")

usethis::use_package("dbscan")

usethis::use_package("data.table")
usethis::use_import_from("data.table", ":=")

usethis::use_package("terra")

usethis::use_package("lidR", min_version = "4.0.0")

usethis::use_package("future", type = "Suggests")

usethis::use_tidy_description()

# Create R files
usethis::use_r(name = "try_to_extract_coordinate_data.R")

usethis::use_r(name = "segment_tree_crowns.R")
