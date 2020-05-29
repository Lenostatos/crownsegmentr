# This file is part of crownsegmentr, an R package for identifying tree crowns
# within 3D point clouds.
#
# Copyright (C) 2020 Leon Steinmeier
# Contact: Lenostatos@gmx.de
#
# crownsegmentr is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# crownsegmentr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with crownsegmentr.  If not, see <http://www.gnu.org/licenses/>.

# Set up renv for package management
install.packages("remotes")
remotes::install_github("rstudio/renv")
renv::init()

# Install packages needed for development
install.packages(c("roxygen2", "testthat", "knitr"))
remotes::install_github("r-lib/devtools")
remotes::install_github("r-lib/usethis")

install.packages("styler")

# Create bare-bones package structure
usethis::create_package(".")

# Ignore this file
usethis::use_build_ignore("package_development_code.R")

# Add a license
usethis::use_gpl3_license(name = "Leon Steinmeier")

# Install packages that this package uses
usethis::use_package("assertthat")

install.packages("dbscan")
usethis::use_package("dbscan")

install.packages("data.table")
usethis::use_package("data.table")

usethis::use_package("BH", type = "LinkingTo")

usethis::use_rcpp()
usethis::use_r(name = "crownsegmentr-package.R")

# The following stackoverflow post provides valuable info on the Makevars file:
# https://stackoverflow.com/questions/43597632/understanding-the-contents-of-the-makevars-file-in-r-macros-variables-r-ma

usethis::use_r(name = "try_to_extract_coordinate_data.R")

usethis::use_r(name = "segment_tree_crowns.R")

# Make usethis functions available on the command line without package prefix.
# usethis::use_usethis()

# Use X after doing Y:
# + Use devtools::load_all() or Ctrl+Shift+L after changing R code to reinstall
#   package
# + Use "Install and restart" in the Build tab or Ctrl+Shift+B after changing
#   C/C++ code in order to do perform an installation with source compilation.
# + Use devtools::document() or Ctrl+Shift+D after changing documentation.
# + Use usethis::use_test() when editing an R source file in order to create a
#   test file (or switch to it if it already exists).
# + Use usethis::use_r() when editing a test file in order to switch back to the
#   respective R source file.
# + Use usethis::use_tidy_style() to style all code according to the tidyverse
#   style guide.
