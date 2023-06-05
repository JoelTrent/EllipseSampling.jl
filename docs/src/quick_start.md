```@setup quick_start
using EllipseSampling
using Plots; gr()
Plots.reset_defaults()
default(palette=:seaborn_colorblind6, msw=0, markeralpha=0.7, aspect_ratio=:equal, label=nothing)
```
# Quick Start

## Equally Spaced Points

To generate 9 equally spaced points on an ellipse, with x radius of 1.0 and y radius of 0.5, no rotation and a centre of ``x,y = [2,1]`` we use:

```example quick_start
using EllipseSampling, Plots

e=construct_ellipse(1.0, 0.5, 0.0, 2.0, 1.0)
points=generate_N_equally_spaced_points(9, e; start_point_shift=0.0) 
scatter(points[1,:], points[2,:])
```

Or alternatively:
```example quick_start
using EllipseSampling, Plots

points=generate_N_equally_spaced_points(9, 1.0, 0.5, 0.0, 2.0, 1.0; start_point_shift=0.0) 
scatter(points[1,:], points[2,:])
```

### Rotated Ellipses

If our ellipse has an anticlockwise rotation of ``\\frac{1}{3}\\pi`` radians or 60 degrees then we modify that argument to [`construct_ellipse`](@ref).

```example quick_start
using EllipseSampling, Plots

e=construct_ellipse(1.0, 0.5, pi/3.0, 2.0, 1.0)
points=generate_N_equally_spaced_points(9, e; start_point_shift=0.0) 
scatter(points[1,:], points[2,:])
```

## Clustered Points

To more easily see the clustering effect we will increase the number of points generated and decrease the y radius. Note, the closer in magnitude the major and minor axis radii are, the weaker the clustering effect. 
To generate 30 clustered points on an ellipse, with x radius of 1.0 and y radius of 0.1, no rotation and a centre of ``x,y = [2,1]`` we use:

```example quick_start
using EllipseSampling, Plots

e=construct_ellipse(1.0, 0.1, 0.0, 2.0, 1.0)
points=generate_N_clustered_points(30, e; start_point_shift=0.0, sqrt_distortion=0.0) 
scatter(points[1,:], points[2,:])
```

The clustering effect becomes weaker when we increase the parameter `sqrt_distortion` towards 1.0:

```example quick_start
using EllipseSampling, Plots

plot()
e=construct_ellipse(1.0, 0.1, 0.0, 2.0, 1.0)
for sqrt_distortion in 0.0:0.2:1.0
    points=generate_N_clustered_points(10, e; start_point_shift=0.0, sqrt_distortion=sqrt_distortion) 
    scatter!(points[1,:], points[2,:], label=string("sqrt_distortion=",sqrt_distortion))
end
```

The clustering effect is completely gone if our ellipse is a circle:

```example quick_start
using EllipseSampling, Plots

plt=plot()
e=construct_ellipse(1.0, 1.0, 0.0, 2.0, 1.0)
for sqrt_distortion in 0.0:0.5:1.0
    points=generate_N_clustered_points(10, e; start_point_shift=0.0, sqrt_distortion=sqrt_distortion) 
    scatter!(plt, points[1,:], points[2,:], label=string("sqrt_distortion=",sqrt_distortion),
            markersize=7-sqrt_distortion*3, markeralpha=0.5)
end
plt
```

## Custom Sampling Method

If we want to use a custom sampling method different to equal spacing and clustered methods then we can use [`generate_perimeter_point`](@ref) with any arbitrary distribution defined on \[0,1\]. For example, if we want to take 100 uniform random samples (uniform with respect to arc length) of our ellipse perimeter we would use:

```example quick_start
using EllipseSampling, Plots
e=construct_ellipse(1.0, 0.5, 0.0, 2.0, 1.0)
N = 100
samples = rand(N)

# wrap e in Ref so that the function correctly broadcasts across samples
points = generate_perimeter_point.(samples, Ref(e)) 
points = reduce(hcat, points)

scatter(points[1,:], points[2,:])
```

