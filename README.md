
<!-- README.md is generated from README.Rmd. Please edit that file -->

# crownsegmentr

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/crownsegmentr)](https://CRAN.R-project.org/package=crownsegmentr)

<!-- badges: end -->

crownsegmentr is an R package with a C++ implementation of the AMS3D
algorithm [(Ferraz et. al, 2016)](#ferraz2016) for tree crown
segmentation in airborne lidar data. For a general description of how
the tree segmentation works, see the documentation of the
[`segment_tree_crowns`](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/R/segment_tree_crowns.R)
generic. Pseudo code of the AMS3D algorithm is listed
[below](#pseudo-code-of-the-ams3d-algorithm).

## Installation

You can install the development version of crownsegmentr from
[GitHub](https://github.com/Lenostatos/crownsegmentr) with:

``` r
# install.packages("pak")
pak::pak("Lenostatos/crownsegmentr")
```

## Example

This is a basic example which shows you how to segment a normalized
point cloud and visualize the result:

``` r
# Load a point cloud of some trees included in the lidR package
point_cloud <- lidR::readLAS(system.file(
  "extdata/MixedConifer.laz",
  package = "lidR"
))

# Segment a normalized point cloud
segmented_point_cloud <- crownsegmentr::segment_tree_crowns(
  point_cloud,
  crown_diameter_to_tree_height = 0.25,
  crown_height_to_tree_height = 0.5
)

# Generate random crown colors
crown_colors <- lidR::pastel.colors(
  n = length(unique(segmented_point_cloud@data[["crown_id"]]))
)

# Plot the segmented crown bodies
lidR::plot(
  segmented_point_cloud,
  color = "crown_id",
  pal = crown_colors,
  nbreaks = length(crown_colors),
  size = 3,
  axis = TRUE
)
```

## Code Structure

The package provides four functions: The central function
[`segment_tree_crowns`](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/R/segment_tree_crowns.R),
the preprocessing functions
[`watershed_diameter_raster`](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/R/diameter_raster.R)
and
[`li_diameter_raster`](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/R/diameter_raster.R),
and the postprocessing function
[`remove_small_trees`](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/R/remove_small_trees.R).
While the latter three are implemented only in R, the code base of
`segment_tree_crowns` is split into an R “front-end” and a C++
“back-end”.

#### Design principles

The function `segment_tree_crowns` is a [S4
generic](https://adv-r.hadley.nz/s4.html#s4-generics) , i.e. it can
handle point cloud data stored in different data types and behave
differently according to that type. More specifically, the generic
function chooses one out of several so called “methods” based on the
input data type. There are methods for

- `data.frame`s/`data.table`s,
- [`lidR::LAS`](https://github.com/Jean-Romain/lidR/blob/HEAD/R/Class-LAS.R)
  objects, and
- [`lidR::LAScatalog`](https://github.com/Jean-Romain/lidR/blob/HEAD/R/Class-LAScatalog.R)s.

While the other functions (`li_diameter_raster`,
`watershed_diameter_raster` and `remove_small_trees`) are also
implemented as S4 generics, they can only take point clouds of type
`lidR::LAS`.

### segment_tree_crowns front-end

The data type specific methods deal with specifics of their data type.
The actual segmentation is done by the internal function
[`segment_tree_crowns_core`](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/R/segment_tree_crowns_core.R),
which is used by the `data.frame`/`data.table` and `lidR::LAS` methods.
The `lidR::LAScatalog` method internally calls the `lidR::LAS` method.

#### segment_tree_crowns_core

This function performs the segmentation by first calling the C++
back-end to calculate centroids and by then clustering the terminal centroids
with the DBSCAN algorithm (as implemented in the
[`dbscan::dbscan`](https://cran.r-project.org/package=dbscan) function).
It takes point cloud data in the tabular form of a `data.frame` or
`data.table` and returns a list with at most three elements. The first
element always contains a vector of crown IDs with one ID for each point
(i.e. row) in the input data. The second and third elements are optional
and contain coordinates of terminal and prior centroids together with crown 
IDs and (row) indices of the points they belong to.

#### data.frame/data.table Method

This method just calls the core function and binds the returned crown
IDs to the input table. It returns this extended table and, if
requested, also the terminal centroids or all centroids calculated by the core
function.

#### lidR::LAS Method

Similar to the `data.frame`/`data.table` method in that it extends the
input object with a crown ID attribute and, if requested, returns the
centroids as separate `lidR::LAS` objects. The centroid objects are assigned 
the metadata of the input object.

#### lidR::LAScatalog Method

For context: The [`lidR` R
package](https://cran.r-project.org/package=lidR) offers a framework for
processing point clouds of large areas, possibly stored in multiple
files and referenced by so called
[`LAScatalog`](https://cran.r-project.org/package=lidR/vignettes/lidR-LAScatalog-class.html)s.
`LAScatalog`s organize point clouds in adjacent chunks which are
processed individually. The chunks each get a buffer area around them so
that edge effects can be accounted for. In the case of individual tree
segmentation, edge effects would be that tree crowns are cut off at the
edge of a chunk when not using a buffer. Parallel processing of the
chunks is also supported.

The `segment_tree_crowns` method for `lidR::LAScatalog`s internally
defines a function which segments a chunk of a `LAScatalog`. This
function is then applied to all chunks of the `LAScatalog` provided by
the user. The “chunk function” first passes the chunk to the `lidR::LAS`
method. Afterwards, it excludes both tree crowns and unsegmented points
in the buffer area. Since crown IDs overlap across chunks (the IDs in
each chunk start at 1), the chunk function also calculates unique
replacements for the crown IDs based on the apices’ absolute
coordinates.

#### Other Internal Functionality

- [validation
  functions](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/R/validation_functions.R)
  for method arguments
- helper functions
  [`extract_coordinate_values`](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/R/extract_coordinate_values.R)
  and `collect_scale_n_offset_of_LAScatalog_files`
- [test
  suite](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/tests/testthat)
  for R functions/methods

### C++ Back-End

The back-end is a small C++ library which implements the AMS3D
algorithm. The core functionality can be found in:

- `namespace ams3d`: Functionality for calculating a point's terminal centroid 
  with the AMS3D algorithm
  ([header](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/inst/include/ams3d.h)),
  and
- `namespace spatial`: a facade to the [Boost
  Geometry](https://www.boost.org/doc/libs/1_75_0/libs/geometry/doc/html/geometry/introduction.html)
  library which provides e.g. the spatial index used for finding points
  inside cylinders
  ([header](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/inst/include/spatial.h)).

There is also some R interface code outside of any namespace called

- `ams3d_R_interface`
  ([header](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/inst/include/ams3d_R_interface.h)).

*Note*: The namespaces `ams3d` and `spatial` contain internal functions,
classes, etc. which should not be used in other namespaces. These
internal components are indicated by an underscore at the beginning of
their name.

#### ams3d_R_interface

This code contains the C++ functions that are callable from R. They are
not contained in any namespace, since this is a requirement of the
[`Rcpp`](https://cran.r-project.org/package=Rcpp) package which does the
actual exposition to R. The functions are

- `calculate_centroids_normalized`
  ([source](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/src/ams3d_R_interface_normalized.cpp))
  for processing normalized point clouds,
- `calculate_centroids_terraneous`
  ([source](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/src/ams3d_R_interface_terraneous.cpp))
  for processing non-normalized point clouds, and
- `calculate_centroids_flexible`
  ([source](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/src/ams3d_R_interface_flexible.cpp))
  for processing both normalized and non-normalized while possibly also
  using rasters for the crown diameter and crown height to tree height
  parameters.

They are all doing basically the same thing, which is looping over
points they get from R and call the functionality exposed by `ams3d` and
`spatial` to calculate terminal centroids for these points.

In addition, there are also some helper functions in the
`namespace ams3d_R_interface_util`
([source](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/src/ams3d_R_interface_util.cpp)).

The `ams3d_R_interface` is the only part of the C++ code which calls
R-specific functions (a.o. it manages a progress bar provided by the R
package [`progress`](https://cran.r-project.org/package=progress)). By
separating the core functionality from the R-specific C++ code it is
possible to use the core functionality with other C++ code when not
using R.

#### namespace ams3d

This namespace only exposes two functions
([header](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/inst/include/ams3d.h)):

- `calculate_terminal_centroid`
- `calculate_all_centroids`

They do exactly the same thing, i.e. calculate the terminal centroid of a 
point, except that the `*_plus_centroids` variant also returns all prior 
centroids which were calculated during the process. Both functions are 
overloaded three times:

1.  The most simple overload assumes a normalized point cloud with
    ground height at zero and takes single numbers for the
    `crown_diameter_to_tree_height` and `crown_height_to_tree_height`
    arguments
    ([source](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/src/ams3d_normalized.cpp)).
2.  The second overload can deal with not normalized point clouds by
    reading ground heights from a raster
    ([source](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/src/ams3d_terraneous.cpp)).
3.  The third overload reads both ground height and cylinder dimension
    parameters from rasters
    ([source](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/src/ams3d_flexible.cpp)).

There is also an internal `_Kernel` class
([source](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/src/ams3d_kernel.cpp))
that models the cylinder used to find points in the neighborhood of a
point or centroid. It also contains the logic to calculate the weighted
centroid of all points inside the cylinder. A `_Kernel` object is
instantiated with a point or centroid and the
crown-diameter-to-tree-height and crown-height-to-tree-height
parameters. It features only one public method:
`calculate_centroid_in( point_cloud )`. Cylinder dimensions are
calculated in the constructor and the centroid calculation logic is
implemented in a few private methods.

There is one more very small internal namespace in `ams3d` called
`_math_functions` that contains the gaussian and epanechnikov functions
used for weighing the points inside a cylinder during centroid
calculation.

#### namespace spatial

This namespace mainly exposes functionality of and based on the C++
library
[`Boost Geometry`](https://www.boost.org/doc/libs/1_75_0/libs/geometry/doc/html/index.html)
for dealing with point data. Most of this functionality consists of data
types
([header](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/inst/include/spatial_types.h))
but there are also some functions and one functor
([header](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/inst/include/spatial_util.h)
and
[source](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/src/spatial_util.cpp)).
Additionally, there is some raster functionality
([header](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/inst/include/spatial_raster.h))
and a few custom iterators
([header](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/inst/include/spatial_index_creation.h))
for inserting points into the spatial index provided by the Boost
Geometry library (an R\*-tree).

The data types are:

- `coordinate_t` for coordinate values,
- `distance_t` for distances,
- `point_2d_t` and `point_3d_t` for 2D and 3D points,
- `index_for_3d_points_t` for R\*-tree index structures, and
- `box_t` for 3D boxes.

A few simple functions of `boost::geometry` are directly forwarded:

- `distance( geometry1, geometry2 )` returns the distance between two
  geometric objects,
- `get_x( point )`, `get_y( point )`, and `get_z( point )` return the
  respective coordinate value of a point, and
- `get_xy_point( point )` returns a point’s xy-coordinates as a 2D
  point.

Exposed functions with own logic are

- `get_points_intersecting_vertical_cylinder( <cylinder dimensions and an index structure> )`
  which searches an R\*-tree index and
- `weighted_mean_of( points, weights )` which calculates a weighted
  average position of a collection of points.

In addition there is an internal functor class called
`_within_xy_distance_functor` whose objects are needed for queries to
the R\*-tree index. Functors are function objects, i.e. objects with a
`()`-operator, making it possible to use such an object like a function.

For handling raster data in the C++ back-end, a small set of three
raster classes was set up
([header](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/inst/include/spatial_raster.h)):

- `I_Raster` an interface class which defines methods that every raster
  object should have.
- `Raster` a normal raster class which implements the methods demanded
  by `I_Raster`.
- `Single_value_pseudo_raster` also implements the methods demanded by
  `I_Raster` but internally stores just one value. It can return this
  value faster than a `Raster` object could if it was used for storing
  just one value.

The raster classes were set up like this to make it possible to pass
either a single value or an actual raster to the same function argument.
Without this it would be necessary to code every combination of function
parameters where the ground height and the
`crown_diameter_to_tree_height` and `crown_height_to_tree_height`
parameters can be either a single value or a raster of values. With the
`I_Raster` class it is possible to pass either
`Single_value_pseudo_raster` or `Raster` objects to the same function
parameter of type `I_Raster`.

The spatial index provided by the Boost Geometry library can be
constructed with a packing algorithm that results in faster queries. In
order to use this feature, the index constructor needs to be passed all
of the points it should use in the form of a begin and end iterator. In
order to skip points with non-finite coordinate values and points below
certain heights, custom iterators for `std::vector< point_3d_t >`
([header](https://github.com/Lenostatos/crownsegmentr/blob/HEAD/inst/include/spatial_index_creation.h))
are used by a few public functions to set up spatial indices:

- `create_index_of_finite( points, min_height )` skips points with
  non-finite coordinate values and points below `min_height`
- `create_index_of_above_ground( points, min_height_above_ground, ground_height_grid )`
  same as above but does not assume a normalized point cloud and instead
  reads ground heights from a raster
- `create_index_of_above_ground( points, min_height_above_ground_grid, ground_height_grid )`
  same as above but reads minimum above-ground heights also from a
  raster

## Coding Practices

### C++

There may be a few syntax constructs in the C++ code which appear
unfamiliar to R users. This section gives the rationale for some of
these constructs. Most of it is based on information found at the
website [learncpp.com](https://www.learncpp.com/).

#### Object Instantiation with `{}`

There are a few different ways to create and assign a value to objects
in C++:

- `int foo = 0; // copy initialization`
- `int foo(0);  // direct initialization`
- `int foo{0};  // list initialization`
- `std::vector<int> foos{1, 2, 3}; // list initialization`

According to
[learncpp.com](https://www.learncpp.com/cpp-tutorial/variable-assignment-and-initialization/),
list initialization is the preferred option. However, as you can see in
the examples, it works a little bit differently for array-like objects
like e.g. `std::vector`s or `Rcpp::List`s. So if you want to create a
vector with a certain size instead of with some elements, you need to
use direct initialization instead,
e.g. `Rcpp::NumericVector foos( <size of the vector> );` for `Rcpp`
vectors.

#### Constructors with Member Initializer Lists

Instances (i.e. objects) of a class are initialized by constructors.
Constructors are basically functions without a return value which
internally assign values to the properties of an object. According to
[learncpp.com](https://www.learncpp.com/cpp-tutorial/constructor-member-initializer-lists/),
the most direct way to do these assignments is to use member initializer
lists. The syntax of these initializer lists looks like this:

    some_class (
        <arguments to constructor>
    ):
        property_1{ <initial value> },
        property_2{ <initial value> } // ...and so on
    {
        <the actual body of the constructor with possibly more code>
    }

## Pseudo Code of the AMS3D Algorithm

    # outer loop over the points in the point cloud
    for each point in the point cloud:
        find_mode_of( point )


    # inner loop to find the mode of a single point
    find_mode_of( point ):

        # this assignment is needed in the loop below
        current_centroid = point
        
        # "move" towards the nearest mode by calculating a new centroid at the 
        # location of the previous centroid and repeating this until the
        # centroids converge
        do:
            former_centroid = current_centroid
            current_centroid = calculate_centroid_of (
                points_in_neighborhood_of( former_centroid )
            )
        while( distance_of( former_centroid, current_centroid ) > very_small 
                AND number_of_iterations < too_many )
      
        # the terminal centroid is returned as approximation to the mode
        return( current_centroid )
        

    # the neighborhood of a point is defined by a crown-sized cylinder
    points_in_neighborhood_of( point ):
        return( points_in( vertical_cylinder_at( point ) ) )
        
        
    # the cylinder's size is calculated using its above-ground height and the
    # four main parameters to the algorithm
    vertical_cylinder_at( point ):
        cylinder = new cylinder (
            height   = above_ground_height_of( point ) * cl_2_th + cl_c
            diameter = above_ground_height_of( point ) * cd_2_th + cd_c
        )
        
        return( upper_three_quarters_of( cylinder ) )

        # The numbers "cl_2_th", "cd_2_th", "cl_c" and "cd_c" are parameters to the
        # algorithm and stand for "crown length to tree height",
        # "crown diameter to tree height", "crown length constant" and
        # "crown diameter constant".
        
        
    points_in( cylinder ):
        <use a spatial index for finding the points in cylinder>


    calculate_centroid_of( points ):
        return( weighted_average_position_of (
            points, 
            weights_of( points )
        ) )


    weights_of( points ):
        return( for each point in points:
            exp( -5 * relative_horizontal_distance_of_cylinder_center_to( point )^2 )
             * (  1 - relative_vertical_distance_of_cylinder_center_to( point )^2 )
        )
        # horizontal: gaussian profile
        # vertical  : epanechnikov profile

## References

<a name="ferraz2016"></a> Ferraz, A.; Saatchi, S.; Mallet, C. & Meyer,
V., (2016) “Lidar detection of individual tree size in tropical
forests”, In: *Remote Sensing of Environment*, 183, 318-333, DOI:
[10.1016/j.rse.2016.05.028](https://doi.org/10.1016/j.rse.2016.05.028)
<a name="ferraz2012"></a> Ferraz, A.; Bretar, F.; Jaquemond, S.;
Gonçalves, G.; Pereira, L.;, Tomé, M. & Soares, P., “3-D mapping of a
multi-layered Mediterranean forest using ALS data”, In: *Remote Sensing
of Environment*, 121, 210-223, DOI:
[10.1016/j.rse.2012.01.020](https://doi.org/10.1016/j.rse.2012.01.020)
