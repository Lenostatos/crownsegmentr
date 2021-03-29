# crownsegmentr
An R package with a C++ implementation of the AMS3D algorithm [(Ferraz et. al, 2016)](#ferraz2016) for tree crown segmentation in airborne lidar data. For a general description of how the tree segmentation works, see the documentation of the [`segment_tree_crowns` generic](R/segment_tree_crowns.R). Pseudo code of the AMS3D algorithm is listed [below](pseudo-code-of-the-ams3d-algorithm).

## Code Structure
The code base is split into an R "front-end" and a C++ "back-end".

### R Front-End
The front-end exposes just one R function called [`segment_tree_crowns`](R/segment_tree_crowns.R). This function is an [S4 generic](https://adv-r.hadley.nz/s4.html#s4-generics), i.e. it can be passed point cloud data stored in different data types and behaves differently according to that type. More specifically, the generic function chooses one out of several so called "methods" based on the input data type. There currently are methods for

- `data.frame`s/`data.table`s,
- `lidR::LAS` objects, and
- `lidR::LAScatalog`s.

All of these methods just deal with specifics of their data type. The actual segmentation is done by the internal core function [`segment_tree_crowns_core`](R/segment_tree_crowns_core.R), which is used by the `data.frame`/`data.table` and `lidR::LAS` methods. The `lidR::LAScatalog` method internally calls the `lidR::LAS` method.

#### `segment_tree_crowns_core`
This function performs the segmentation by first calling the C++ back-end to calculate modes and by then clustering these modes with the DBSCAN algorithm (as implemented in the [`dbscan::dbscan` function](https://cran.r-project.org/package=dbscan)). It takes point cloud data in the tabular form of a `data.frame` or `data.table` and returns a list with at most three elements. The first element always contains a vector of crown IDs with one ID for each point (i.e. row) in the input data. The second and third elements are optional and contain mode and centroid coordinates together with crown IDs and (row) indices of the points they belong to.

#### `data.frame`/`data.table` Method
This method just calls the core function and binds the returned crown IDs to the input table. It returns this extended table and, if requested, also the modes and/or centroids returned by the core function.

#### `lidR::LAS` Method
Similar to the `data.frame`/`data.table` method in that it extends the input object with a crown ID attribute and, if requested, returns the modes and centroids as separate `lidR::LAS` objects. The mode and centroid objects are assigned the metadata of the input point cloud.

#### `lidR::LAScatalog` Method
For context: The [`lidR` R package](https://cran.r-project.org/package=lidR) offers a framework for processing point clouds of large areas, possibly stored in multiple files and referenced by so called [`LAScatalog`s](https://cran.r-project.org/package=lidR/vignettes/lidR-LAScatalog-class.html). `LAScatalog`s organize point clouds in adjacent chunks which are processed individually. The chunks each get a buffer area around them so that edge effects can be accounted for. In this case, edge effects would be that tree crowns are cut off at the edge of a chunk when not using a buffer. Parallel processing of the chunks is also supported.

The `segment_tree_crowns` method for `lidR::LAScatalog`s internally defines a function which segments a chunk of a `LAScatalog`. This function is then applied to all chunks of the `LAScatalog` provided by the user. The "chunk function" first passes the chunk to the `lidR::LAS` method. Afterwards, it excludes both tree crowns and unsegmented points in the buffer area. Since crown IDs overlap across chunks (the IDs in each chunk start at 1), the chunk function also calculates unique replacements for the crown IDs based on the apices' absolute coordinates.

#### Other Internal Functionality
- [validation functions](R/validation_functions.R) for method arguments
- helper functions [`extract_coordinate_values`](R/extract_coordinate_values.R) and `collect_scale_n_offset_of_LAScatalog_files`
- [test suite](tests/testthat) for R functions/methods


### C++ Back-End
The back-end is a small C++ library which implements the AMS3D algorithm. The core functionality can be found in:

- `namespace ams3d`: An implementation for calculating a single mode with the AMS3D algorithm ([header](inst/include/ams3d.h) and [source](src/ams3d.cpp)), and
- `namespace spatial`: a facade to the [Boost Geometry](https://www.boost.org/doc/libs/1_75_0/libs/geometry/doc/html/geometry/introduction.html) library which provides e.g. the spatial index used for finding points inside cylinders ([header](inst/include/spatial.h) and [source](src/spatial.cpp)).

There is also some interface code outside of any namespace called "ams3d_R_interface" ([header](inst/include/ams3d_R_interface.h) and [source](src/ams3d_R_interface.cpp)). This code loops over the points which it gets from R and calls the C++ functions exposed by `ams3d` and `spatial` to calculate modes for these points. This interface is not contained in any namespace, since this is a requirement of the [`Rcpp`](https://cran.r-project.org/package=Rcpp) package which is used to actually connect the interface to R. It is also the only part of the C++ code which calls R-specific functions (a.o. it manages a progress bar provided by the [R package `progress`](https://cran.r-project.org/package=progress)). By separating the core functionality from the R-specific C++ code it is possible to use the core functionality with other C++ code when not using R.

*Note*: The namespaces `àms3d` and `spatial` contain internal functions, classes, etc. which should not be used in other namespaces. These internal components are indicated by an underscore at the beginning of their name.

#### `namespace ams3d`
This namespace only exposes two functions:

- `calculate_a_single_mode`
- `calculate_a_single_mode_plus_centroids`

They do exactly the same thing, i.e. calculate the mode of a point, except that the `*_plus_centroids` variant also returns the centroids which were calculated during the process. 

There is also an internal `_Kernel` class that models the cylinder used to find points in the neighborhood of a point or centroid. It also contains the logic to calculate the weighted centroid of all points inside the cylinder. A `Kernel` object is instantiated with a point or centroid and the crown-diameter-to-tree-height and crown-height-to-tree-height parameters. It features only one public method: `calculate_centroid_in( point_cloud )`. Cylinder dimensions are calculated in the constructor and the centroid calculation logic is implemented in a few private methods.

There is one more very small internal namespace in `ams3d` called `_math_functions` that contains the gaussian and epanechnikov functions used for weighing the points inside a cylinder during centroid calculation.

#### `namespace spatial`
This namespace exposes functionality of and based on `boost::geometry` for dealing with point data. Most of this functionality consists of data types but there are also some functions and one functor.

The data types are:

- `coordinate_t` for coordinate values,
- `distance_t` for distances,
- `point_2d_t` and `point_3d_t` for 2D and 3D points,
- `index_for_3d_points_t` for R*-tree index structures, and
- `box_t` for 3D boxes.

There are also a few simple functions which are directly forwarded from `boost::geometry`:

- `distance( geometry1, geometry2 )` returns the distance between two geometric objects and
- `get_x( point )`, `get_y( point )`, and `get_z( point )` return the respective coordinate values of a point.

Exposed functions with own logic are

- `create_index_of( points )` which creates an R*-tree index,
- `get_points_intersecting_vertical_cylinder( <cylinder dimensions and an index structure> )` which searches an R*-tree index, and
- `weighted_mean_of( points, weights )` which calculates a weighted average position of a collection of points.

There is one more internal functor class called `_within_xy_distance_functor` whose objects are needed for queries to the R*-tree index. Functors are objects with a `()`-operator, i.e. they are some kind of function objects because they can be used like functions with that operator.


## Coding Practices

### C++
There may be a few syntax constructs in the C++ code which appear unfamiliar to R users. This section gives the rationale for some of these constructs. Most of it is based on information found at the [learncpp.com](https://www.learncpp.com/) website.

#### Object Instantiation with `{}`
There are a few different ways to create and assign a value to objects in C++:

- `int foo = 0; // copy initialization`
- `int foo(0);  // direct initialization`
- `int foo{0};  // list initialization`
- `std::vector<int> foos{1, 2, 3}; // list initialization`

According to [learncpp.com](https://www.learncpp.com/cpp-tutorial/variable-assignment-and-initialization/), list initialization is the preferred option. However, as you can see in the examples, it works a little bit differently for array-like objects like e.g. `std::vector`s or `Rcpp::List`s. So if you want to create a vector with a certain size instead of with some elements, you need to use direct initialization instead, e.g. `Rcpp::NumericVector foos( <size of the vector> );` for `Rcpp` vectors.

#### Constructors with Member Initializer Lists
Instances (i.e. objects) of a class are initialized by constructors. Constructors are basically functions without a return value which internally assign values to the properties of an object. According to [learncpp.com](https://www.learncpp.com/cpp-tutorial/constructor-member-initializer-lists/) the most direct way to do these assignments is to use member initializer lists. The syntax of these initializers looks like this:

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
      
        # the last centroid is returned as the original point's mode
        return( current_centroid )
        
    
    # the neighborhood of a point is defined by a crown-sized cylinder
    points_in_neighborhood_of( point ):
        return( points_in( vertical_cylinder_at( point ) ) )
        
        
    # the cylinder's size is calculated using its above-ground height and the
    # two main parameters to the algorithm
    vertical_cylinder_at( point ):
        cylinder = new cylinder (
            height   = above_ground_height_of( point ) * ch_2_th
            diameter = above_ground_height_of( point ) * cd_2_th
        )
        
        return( upper_three_quarters_of( cylinder ) )

        # ch_2_th and cd_2_th are available as parameters to the
        # algorithm and stand for "crown height to tree height" and
        # "crown diameter to tree height"
        
        
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
<a name="ferraz2016"></a>
    Ferraz, A.; Saatchi, S.; Mallet, C. & Meyer, V.,
    "Lidar detection of individual tree size in tropical forests"",
    In: *Remote Sensing of Environment*, Elsevier BV, 2016, 183, 318-333,
    DOI: [10.1016/j.rse.2016.05.028](https://doi.org/10.1016/j.rse.2016.05.028)
