# crownsegmentr 1.1.0
* Fixed a mathematical error in watershed_diameter_raster and li_diameter_raster: diameter is now calculated as 2*sqrt(area/pi), and no longer as sqrt(area).
This causes resulting raster values to be slightly higher in most cases.
* In the light of new research results, the default value for the parameter smoothing_radius in the above functons was changed from 5 to 10.
* These two changes now mean that watershed_diameter_raster and li_diameter_raster produce different results compared to prior versions

# crownsegmentr 1.0.2
* Fixed malfunctioning in diameter_raster functions that stemmed from terra behavior differing on different machines of parallelization regimes

# crownsegmentr 1.0.1
* Submission with all examples CRAN compatible

# crownsegmentr 1.0.0

* Initial CRAN submission.
