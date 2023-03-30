# EllipseSampling

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JoelTrent.github.io/EllipseSampling.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JoelTrent.github.io/EllipseSampling.jl/dev/)
[![Build Status](https://github.com/JoelTrent/EllipseSampling.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JoelTrent/EllipseSampling.jl/actions/workflows/CI.yml?query=branch%3Amain)

A set of helpful functions for calculating points on the boundary defined by an ellipse built around Julia implementations of two functions by [John D. Cook](https://www.johndcook.com/blog/2022/11/02/ellipse-rng/). It handles cases where the ellipse has been rotated and/or translated, and where either of the x and y axes is the major or minor axis. Built around the desire to generate N equally-spaced points on the boundary defined by an ellipse. 