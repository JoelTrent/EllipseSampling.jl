```@index
Pages = ["user_interface.md"]
```

# User Interface

## Ellipse Construction

For ease of use, all ellipse sampling \[generation\] functions can either be called directly with ellipse parameters such as a, b and α, or we can construct an [`EllipseSampling.Ellipse`](@ref) struct first with those parameters, and use this struct when calling the generation functions.

We construct an [`EllipseSampling.Ellipse`](@ref) struct using [`construct_ellipse`](@ref).

```@docs
construct_ellipse
```

## Sampling Methods

Three sampling methods, [`generate_N_equally_spaced_points`](@ref), [`generate_N_clustered_points`](@ref) and [`generate_perimeter_point`](@ref), are provided for use.

1. [`generate_N_equally_spaced_points`](@ref) will generate points uniformly with respect to arc length along the perimeter of an ellipse. 
2. [`generate_N_clustered_points`](@ref) will generate points uniformly with respect to arc length along the perimeter of a distorted version of the input ellipse and then map these onto the input ellipse. The distorted version of the input ellipse varies the major axis radius of the input ellipse between being equal to the minor axis radius and the original major axis radius, allowing a stronger clustering of points around the vertexes of the ellipse on the major axis (i.e. the region of greatest curvature). This means it is a generalised version of [`generate_N_equally_spaced_points`](@ref). 
3. [`generate_perimeter_point`](@ref) will generate a singular point on the ellipse's perimeter given a specified distance along the normalised perimeter ∈ \[0,1\]. This is useful for defining a custom sampling method not covered by the previous two sampling methods. 

Example use of these functions and visualisation of their outputs is shown in [Quick Start](@ref)

```@docs
generate_N_equally_spaced_points
generate_N_clustered_points
generate_perimeter_point
```

## Ellipse Angle Calculation and Parameteric Equations

```@docs 
t_from_arclength
t_from_arclength_general
x_parametric_equation 
y_parametric_equation
```