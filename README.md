# EllipseSampling.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JoelTrent.github.io/EllipseSampling.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JoelTrent.github.io/EllipseSampling.jl/dev/)
[![Build Status](https://github.com/JoelTrent/EllipseSampling.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JoelTrent/EllipseSampling.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/JoelTrent/EllipseSampling.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JoelTrent/EllipseSampling.jl/)

# Description

A lightweight package for generating points on the boundary defined by an 2D ellipse built around Julia implementations of two functions by [John D. Cook](https://www.johndcook.com/blog/2022/11/02/ellipse-rng/). It handles cases where the ellipse has been rotated and/or translated, and where either of the x and y axes is the major or minor axis. It provides methods to generate N equally-spaced and clustered points on the ellipse's boundary as well as a method to sample singular points uniformly on the boundary. Resultantly, any distribution defined on [0, 1] can be used to sample points. Calculation of the arc length and circumference of an ellipse uses [Elliptic.jl](https://github.com/nolta/Elliptic.jl).

This package's creation was motivated by the need to sample points on the boundary defined by an elliptical approximation of the log-likelihood function around the maximum likelihood estimate of a mechanistic model at a particular confidence level. Resultantly, it also provides a method to convert the matrix representation of this approximation into the parameters of the equivalent ellipse. Points can then be sampled from this ellipse.

## Getting Started: Installation And First Steps

To install the package, use the following command inside the Julia REPL:

```julia
using Pkg
Pkg.add(url="https://github.com/JoelTrent/EllipseSampling.jl")
```

To load the package, use the command:

```julia
using EllipseSampling
```
