```@meta
CurrentModule = EllipseSampling
```

# EllipseSampling.jl

[EllipseSampling.jl](https://github.com/JoelTrent/EllipseSampling.jl) is a lightweight package for generating points on the boundary defined by an 2D ellipse built around Julia implementations of two functions by [John D. Cook](https://www.johndcook.com/blog/2022/11/02/ellipse-rng/). It handles cases where the ellipse has been rotated and/or translated, and where either of the x and y axes is the major or minor axis. It provides methods ([`generate_N_equally_spaced_points`](@ref)) and ([`generate_N_clustered_points`](@ref)) to generate N equally-spaced points and clustered points respectively on the ellipse's boundary as well as a method ([`generate_perimeter_point`](@ref)) to sample singular points uniformly on the boundary. Resultantly, any distribution defined on [0, 1] can be used to sample points. Calculation of the arc length and circumference of an ellipse uses [Elliptic.jl](https://github.com/nolta/Elliptic.jl).

This package's creation was motivated by the need to sample points on the boundary defined by an elliptical approximation of the log-likelihood function around the maximum likelihood estimate of a mechanistic model at a particular confidence level. Resultantly, it also provides a method ([`EllipseSampling.calculate_ellipse_parameters`](@ref)) to convert the matrix representation of this approximation into the parameters of the equivalent ellipse. Points can then be sampled from this ellipse.

To get started with the package see [Quick Start](@ref).

A deeper dive into the user interface and internal library can be found in [User Interface](@ref) and [Internal Library](@ref), respectively.

## Getting Started: Installation And First Steps

To install the package, use the following command inside the Julia REPL:

```julia
using Pkg
Pkg.add("EllipseSampling")
```

To load the package, use the command:

```julia
using EllipseSampling
```