# crownsegmentr 1.1.0
* Default values for several parameters in segment_tree_crowns were changed, due to results from a large sensitivity analysis experiment. This applies to these parameters:
segment_crowns_only_above, 
centroid_convergence_distance, 
max_iterations_per_point,
dbscan_neighborhood_radius,
min_num_points_per_crown
* Fixed a mathematical error in watershed_diameter_raster and li_diameter_raster: diameter is now calculated as 2*sqrt(area/pi), and no longer as sqrt(area).
This causes resulting raster values to be slightly higher in most cases.
* In the light of new research results, the default value for the parameter smoothing_radius in the above functons was changed from 5 to 10.
* These changes mean that watershed_diameter_raster and li_diameter_raster now produce different results compared to prior versions, and segment_tree_crowns produces different results if one of the aforementioned parameters is left as default.

# crownsegmentr 1.0.2
* Fixed malfunctioning in diameter_raster functions that stemmed from terra behavior differing on different machines of parallelization regimes

# crownsegmentr 1.0.1
* Submission with all examples CRAN compatible

# crownsegmentr 1.0.0

* Initial CRAN submission.
