# Julia version of functions from https://www.johndcook.com/blog/2022/11/02/ellipse-rng/
"""
    E_inverse(em::T, z::T, m::T) where T<:Float64

Julia version of the python function `t_from_length` by [John D. Cook](https://www.johndcook.com/blog/2022/11/02/ellipse-rng/).

# Arguments
- `em`: complete elliptic integral of the second kind evaluated for the eccentricity squared of the ellipse: `Elliptic.E(m)`. See: [Elliptic.jl](https://github.com/nolta/Elliptic.jl).
- `z`: Difference between `em` and the quotient of arc length and the major axis radius, `em - arc_len/a`.
- `m`: the eccentricity of the ellipse squared. See [`EllipseSampling.eccentricity_squared`](@ref).
"""
function E_inverse(em::T, z::T, m::T) where T<:Float64
    t = (z/em)*(pi/2)
    f(y) = Elliptic.E(y, m) - z
    r = Roots.find_zero(f, t, Roots.Order0())
    return r
end

"""
    t_from_arclength(arc_len::Float64, e::Ellipse)

Calculates the angle t, between 0 and 2π radians, of the location on an unrotated ellipse, given an arc length, `arc_len`, anticlockwise from the positive major axis along the perimeter of the ellipse. The ellipse's x axis is the major axis. It is recommended to call [`t_from_arclength(arc_len::T, a::T, b::T) where T<:Float64`] rather than this function as it handles the case where the major axis of the ellipse is the y axis.

Julia version of the python function `t_from_length` by [John D. Cook](https://www.johndcook.com/blog/2022/11/02/ellipse-rng/).

# Arguments
- `arc_len`: arc length, between 0.0 and the circumference of the ellipse, anticlockwise from the positive major axis along the perimeter of the ellipse.
- `e`: a valid [`EllipseSampling.Ellipse`](@ref) struct which defines an ellipse.
"""
function t_from_arclength(arc_len::Float64, e::Ellipse)
    em = Elliptic.E(e.m)
    t = 0.5*pi - E_inverse(em, em - arc_len/e.a, e.m)
    return t
end

"""
    t_from_arclength(arc_len::T, a::T, b::T) where T<:Float64

An alternate way to call [`t_from_arclength(arc_len::Float64, e::Ellipse)`](@ref).

# Arguments
- `arc_len`: arc length, between 0.0 and the circumference of the ellipse, anticlockwise from the positive major axis along the perimeter of the ellipse.
- `a`: radius of the ellipse's major axis.
- `b`: radius of the ellipse's minor axis.
"""
function t_from_arclength(arc_len::T, a::T, b::T) where T<:Float64
    m = eccentricity_squared(a,b)
    em = Elliptic.E(m)
    t = 0.5*pi - E_inverse(em, em - arc_len/a, m)
    return t
end
#######################################################################################

"""
    t_from_arclength_general(arc_len::Float64, e::Ellipse)

Generalised version of [`t_from_arclength(arc_len::Float64, e::Ellipse)`](@ref) which handles cases where either of the x and y axes are the major axis.

# Arguments
- `arc_len`: arc length, between 0.0 and the circumference of the ellipse, anticlockwise from the positive major axis along the perimeter of the ellipse.
- `e`: a valid [`EllipseSampling.Ellipse`](@ref) struct which defines an ellipse.
"""
function t_from_arclength_general(arc_len::Float64, e::Ellipse)
    if e.x_radius < e.y_radius
        return t_from_arclength(arc_len, e) + 0.5*pi
    else
        return t_from_arclength(arc_len, e) 
    end
end

"""
    t_from_arclength_general(arc_len::T, a::T, b::T, x_radius::T, y_radius::T) where T<:Float64

Generalised version of [t_from_arclength(arc_len::T, a::T, b::T) where T<:Float64](@ref) which handles cases where either of the x and y axes are the major axis.

# Arguments
- `arc_len`: arc length, between 0.0 and the circumference of the ellipse, anticlockwise from the positive major axis along the perimeter of the ellipse.
- `a`: radius of the ellipse's major axis.
- `b`: radius of the ellipse's minor axis.
- `x_radius`: radius of the ellipse in the x axis (i.e. when the rotation, `α`, is zero).
- `y_radius`: radius of the ellipse in the y axis (i.e. when the rotation, `α`, is zero).
"""
function t_from_arclength_general(arc_len::T, a::T, b::T, x_radius::T, y_radius::T) where T<:Float64
    if x_radius < y_radius
        return t_from_arclength(arc_len, a, b) + 0.5*pi
    else
        return t_from_arclength(arc_len, a, b) 
    end
end

"""
    generate_perimeter_point(norm_distance_on_perimeter::Float64, e::Ellipse)

Generates a single point on an ellipse defined by the parameters contained within `e`, at distance `norm_distance_on_perimeter` ``\\times`` `e.circumference` around the circumference. The point is returned as a vector of length two.

# Arguments
- `norm_distance_on_perimeter`: a number ∈ [0,1] which represents the normalised distance on the perimeter of an ellipse. A value of `0.5` corresponds to a point halfway along the ellipse's perimeter, while a value of `0.7` corresponds to a point 70% along the ellipse's perimeter.
- `e`: a valid [`EllipseSampling.Ellipse`](@ref) struct which defines an ellipse.

# Details

Points are sampled anti-clockwise around the ellipse starting from the most positive vertex on the major axis prior to any rotation. After a rotation is applied, this vertex may no longer be the most positive vertex on the major axis. Vertex here means the two endpoints that lie on the major axis and intersect the ellipse. 

Consider a ellipse that has no translation or rotational component. If the x radius is `2.0` and the y radius is `1.0`, then the x axis is the major axis. Resultantly, a value of `norm_distance_on_perimeter=0.0` will correspond to point `[x,y] = [2.0, 0.0]`. A value of `norm_distance_on_perimeter=0.25`, which corresponds to a ``\\frac{π}{2}`` anti-clockwise rotation around the ellipse will correspond to point `[x,y] = [0.0, 1.0]`. 

```julia
using EllipseSampling
e = construct_ellipse(2.0, 1.0)
generate_perimeter_point(0.0, e)
generate_perimeter_point(0.25, e)
```

Similarly, if the x radius is `1.0` and the y radius is `2.0`, then the y axis is the major axis. Resultantly, a value of `norm_distance_on_perimeter=0.0` will correspond to point `[x,y] = [0.0, 2.0]`. A value of `norm_distance_on_perimeter=0.25`, which corresponds to a ``\\frac{π}{2}`` anti-clockwise rotation around the ellipse will correspond to point `[x,y] = [-1.0, 0.0]`. 

```julia
e = construct_ellipse(1.0, 2.0)
generate_perimeter_point(0.0, e)
generate_perimeter_point(0.25, e)
```

In the event that the ellipse is a circle (the major and minor axis have the same radius), points will be sampled anti-clockwise around the ellipse starting from the most positive point on the x radius (prior to any rotation). 

To demonstrate how sampling works when the ellipse is rotated, for simplicity, consider an ellipse that is a circle with radius 1.0, no translational component and an anti-clockwise rotation of ``\\frac{π}{2}``. In this case the x axis has been rotated anti-clockwise by 90 degrees and so a value of `norm_distance_on_perimeter=0.0` will correspond to point `[x,y] = [0.0, 1.0]`. A value of `norm_distance_on_perimeter=0.25`, which corresponds to a ``\\frac{π}{2}`` anti-clockwise rotation around the rotated ellipse will correspond to point `[x,y] = [-1.0, 0.0]`. 

```julia
e = construct_ellipse(1.0, 1.0, pi/2)
generate_perimeter_point(0.0, e)
generate_perimeter_point(0.25, e)
```

## Random sampling

This function can be easily used to generate uniform random samples from an ellipse by first sampling N points from a uniform distribution defined on [0,1] and then calling this function for each point 1:N.

For example using:
```julia
using EllipseSampling
e = construct_ellipse(2.0,1.0)
N = 100
samples = rand(N)
points = generate_perimeter_point.(samples, Ref(e))
```

Note, here we wrap the ellipse struct `e` in `Ref` so that Julia does not try to broadcast over `e` as well. To get the points generated in this fashion in the same column-wise format as points generated by [`generate_N_equally_spaced_points`](@ref) we can use:

```julia
points = reduce(hcat, points)
```

Other distributions defined on [0,1] can be used to generate points on the ellipse's perimeter in a similar fashion.
"""
function generate_perimeter_point(norm_distance_on_perimeter::Float64, e::Ellipse)

    (0.0 ≤ norm_distance_on_perimeter && norm_distance_on_perimeter ≤ 1.0) || throw(DomainError("norm_distance_on_perimeter must be between 0.0 and 1.0."))
    
    point = zeros(2)
    angle = t_from_arclength_general(e.circumference*norm_distance_on_perimeter, e)

    point[:] .= x_parametric_equation(angle, e), y_parametric_equation(angle, e)

    return point
end

"""
    generate_perimeter_point(norm_distance_on_perimeter::T, x_radius::T, y_radius::T, α::T=0.0, Cx::T=0.0, Cy::T=0.0) where T<:Float64

An alternate way to call [`generate_perimeter_point(norm_distance_on_perimeter::Float64, e::Ellipse)`](@ref), by supplying the parameters of the ellipse to generate a single point on.

# Arguments
- `norm_distance_on_perimeter`: a number ∈ [0,1] which represents the normalised distance on the perimeter of an ellipse. A value of 0.5 corresponds to a point halfway along the ellipse's perimeter, while a value of 0.7 corresponds to a point 70% along the ellipse's perimeter.
- `x_radius`: radius of the ellipse in the x axis (i.e. when the rotation, `α`, is zero).
- `y_radius`: radius of the ellipse in the y axis (i.e. when the rotation, `α`, is zero).
- `α`: an angle in radians (0 to 2π) that the ellipse has been rotated by. A positive value represents an anti-clockwise rotation. Default is `0.0`.
- `Cx`: the x coordinate of the centre of the ellipse (the translation of the ellipse in the x axis). Default is `0.0`.
- `Cy`: the y coordinate of the centre of the ellipse (the translation of the ellipse in the y axis). Default is `0.0`.
"""
function generate_perimeter_point(norm_distance_on_perimeter::T, x_radius::T, y_radius::T, α::T=0.0, Cx::T=0.0, Cy::T=0.0) where T<:Float64

    e = construct_ellipse(x_radius, y_radius, α, Cx, Cy)
    return generate_perimeter_point(norm_distance_on_perimeter, e)
end

"""
    generate_perimeter_point(norm_distance_on_perimeter::Float64, Γ::Matrix{Float64}, θmle::Vector{Float64}, ind1::Int, ind2::Int; confidence_level::Float64=0.01)

An alternate way to call [`generate_perimeter_point(norm_distance_on_perimeter::Float64, e::Ellipse)`](@ref), by supplying a square matrix Γ, the inverse of the Hessian of a log-likelihood function at its maximum likelihood estimate, indexes of the two variables of interest and the confidence level that represent a 2D ellipse approximation of the log-likelihood function.

# Arguments
- `norm_distance_on_perimeter`: a number ∈ [0,1] which represents the normalised distance on the perimeter of an ellipse. A value of 0.5 corresponds to a point halfway along the ellipse's perimeter, while a value of 0.7 corresponds to a point 70% along the ellipse's perimeter.
- `Γ`: a square matrix (2D) which is the inverse of the Hessian of a log-likelihood function at its maximum likelihood estimate.
- `θmle`: the maximum likelihood estimate for the parameters.
- `ind1`: index of the first parameter of interest (corresponds to the row and column index of `Γ`)
- `ind2`: index of the second parameter of interest (corresponds to the row and column index of `Γ`).

# Keyword Arguments
- `confidence_level`: the confidence level ∈ [0.0,1.0] at which the ellipse approximation is constructed. Default is `0.01`.
"""
function generate_perimeter_point(norm_distance_on_perimeter::Float64, Γ::Matrix{Float64}, θmle::Vector{Float64}, ind1::Int, ind2::Int; confidence_level::Float64=0.01)

    _, _, x_radius, y_radius, α = calculate_ellipse_parameters(Γ, ind1, ind2, confidence_level)
    return generate_perimeter_point(norm_distance_on_perimeter, x_radius, y_radius, α, θmle[ind1], θmle[ind2])
end

"""
    generate_N_equally_spaced_points(num_points::Int, e::Ellipse; start_point_shift::Float64=rand())

Generates `num_points` equally spaced points on an ellipse defined by the parameters contained within `e`. The points are returned as an array with two rows and `num_points` columns, with each point stored in a column.

# Arguments
- `num_points`: a positive integer number of points to generate that are equally spaced on the ellipse. 
- `e`: a valid [`EllipseSampling.Ellipse`](@ref) struct which defines an ellipse.

# Keyword Arguments
- `start_point_shift`: a number ∈ [0,1]. Default is `rand()` (defined on [0,1]), meaning that, by default, every time this function is called a different set of points will be generated.

# Details

Points are sampled anti-clockwise around the ellipse starting from the most positive vertex on the major axis prior to any rotation. After a rotation is applied, this vertex may no longer be the most positive vertex on the major axis. Vertex here means the two endpoints that lie on the major axis and intersect the ellipse. If `start_point_shift=0.0` then the first point on the equally spaced points will be placed on the most positive vertex on the major axis prior to any rotation.

Equal spacing of points on the ellipse is with respect to arc length. The position of the first point placed on the ellipse can be shifted by `start_point_shift`, defined on [0.0,1.0], allowing the set of possible points generated for a given `num_points` to cover the full perimeter. This shift is normalised so that when `start_point_shift=1.0`, the position of the first point placed on the ellipse is equal to the position of the second point placed on the ellipse when `start_point_shift=0.0`.
"""
function generate_N_equally_spaced_points(num_points::Int, e::Ellipse; 
    start_point_shift::Float64=rand())

    (0.0 ≤ start_point_shift && start_point_shift ≤ 1.0) || throw(DomainError("start_point_shift must be between 0.0 and 1.0."))
    0 < num_points || throw(DomainError("the number of points, num_points, to generate on the ellipse must be positive"))

    points = zeros(2,num_points)

    shift = start_point_shift/num_points

    lengths = collect(LinRange((shift)*e.circumference, (1+shift)*e.circumference, num_points+1))[1:end-1]
    angles = t_from_arclength_general.(lengths, Ref(e))
    
    for i in 1:num_points
        points[:,i] .= x_parametric_equation(angles[i], e), 
            y_parametric_equation(angles[i], e)
    end

    return points
end

"""
    generate_N_equally_spaced_points(num_points::Int, x_radius::T, y_radius::T, α::T=0.0, Cx::T=0.0, Cy::T=0.0; start_point_shift::Float64=rand()) where T<:Float64

An alternate way to call [`generate_N_equally_spaced_points(num_points::Int, e::Ellipse; start_point_shift::Float64=rand())`](@ref), by supplying the parameters of the ellipse to generate points on.

# Arguments
- `num_points`: a positive integer number of points to generate that are equally spaced on the ellipse. 
- `x_radius`: radius of the ellipse in the x axis (i.e. when the rotation, `α`, is zero).
- `y_radius`: radius of the ellipse in the y axis (i.e. when the rotation, `α`, is zero).
- `α`: an angle in radians (0 to 2π) that the ellipse has been rotated by. A positive value represents an anti-clockwise rotation. Default is `0.0`.
- `Cx`: the x coordinate of the centre of the ellipse (the translation of the ellipse in the x axis). Default is `0.0`.
- `Cy`: the y coordinate of the centre of the ellipse (the translation of the ellipse in the y axis). Default is `0.0`.

# Keyword Arguments
- `start_point_shift`: a number ∈ [0,1]. Default is `rand()` (defined on [0,1]), meaning that, by default, every time this function is called a different set of points will be generated.
"""
function generate_N_equally_spaced_points(num_points::Int, x_radius::T, y_radius::T, α::T=0.0, Cx::T=0.0, Cy::T=0.0; 
    start_point_shift::Float64=rand()) where T<:Float64

    e = construct_ellipse(x_radius, y_radius, α, Cx, Cy)
    return generate_N_equally_spaced_points(num_points, e, start_point_shift=start_point_shift)
end

"""
    generate_N_equally_spaced_points(num_points::Int, Γ::Matrix{Float64}, θmle::Vector{Float64}, ind1::Int, ind2::Int; confidence_level::Float64=0.01, start_point_shift::Float64=rand())

An alternate way to call [`generate_N_equally_spaced_points(num_points::Int, e::Ellipse; start_point_shift::Float64=rand())`](@ref), by supplying a square matrix Γ, the inverse of the Hessian of a log-likelihood function at its maximum likelihood estimate, indexes of the two variables of interest and the confidence level that represent a 2D ellipse approximation of the log-likelihood function.

# Arguments
- `num_points`: a positive integer number of points to generate that are equally spaced on the ellipse. 
- `Γ`: a square matrix (2D) which is the inverse of the Hessian of a log-likelihood function at its maximum likelihood estimate.
- `θmle`: the maximum likelihood estimate for the parameters.
- `ind1`: index of the first parameter of interest (corresponds to the row and column index of `Γ`)
- `ind2`: index of the second parameter of interest (corresponds to the row and column index of `Γ`).

# Keyword Arguments
- `confidence_level`: the confidence level ∈ [0.0,1.0] at which the ellipse approximation is constructed. Default is `0.01`.
- `start_point_shift`: a number ∈ [0,1]. Default is `rand()` (defined on [0,1]), meaning that, by default, every time this function is called a different set of points will be generated.
"""
function generate_N_equally_spaced_points(num_points::Int, Γ::Matrix{Float64}, θmle::Vector{Float64}, ind1::Int, ind2::Int; 
    confidence_level::Float64=0.01, start_point_shift::Float64=rand())

    _, _, x_radius, y_radius, α = calculate_ellipse_parameters(Γ, ind1, ind2, confidence_level)
    return generate_N_equally_spaced_points(num_points, x_radius, y_radius, α, θmle[ind1], θmle[ind2],
                                    start_point_shift=start_point_shift)
end

"""
    generate_N_clustered_points(num_points::Int, e::Ellipse; start_point_shift::Float64=rand(), sqrt_distortion::Float64=0.0)

Generates `num_points` spaced points on an ellipse defined by the parameters contained within `e`. The points are returned as an array with two rows and `num_points` columns, with each point stored in a column.

# Arguments
- `num_points`: a positive integer number of points to generate that are spaced on the ellipse. 
- `e`: a valid [`EllipseSampling.Ellipse`](@ref) struct which defines an ellipse.

# Keyword Arguments
- `start_point_shift`: a number ∈ [0,1]. Default is `rand()` (defined on [0,1]), meaning that, by default, every time this function is called a different set of points will be generated.
- `sqrt_distortion`: a number ∈ [0,1]. Default is `0.0`, meaning that, by default, this function will evenly space points on the the ellipse `e` with respect to the parameter `t`.

# Details

Points are sampled in the same fashion as in [`generate_N_equally_spaced_points(num_points::Int, e::Ellipse; start_point_shift::Float64=rand())`](@ref), except that the spacing of them is dependent on the parameter `sqrt_distortion`. If `sqrt_distortion == 1.0` then these points are equally spaced with respect to the arc length. If `sqrt_distortion == 0.0` then these points are equally spaced with respect to the parameter `t` in the parameteric equations for the x and y locations on an ellipse ([`x_parametric_equation(t::T, e::Ellipse) where T<:Float64`](@ref) and [`y_parametric_equation(t::T, e::Ellipse) where T<:Float64`](@ref)). Varying this parameter between `0.0` and `1.0` allows all of the spacings between these two extremes to be obtained.

The impact of using this function is that for all values of `sqrt_distortion < 1.0`, points generated will be more clustered around the vertexes of the ellipse on the major axis (i.e. the region of greatest curvature). The strength of this clustering increases as `sqrt_distortion` ``\\to 0.0``. The strength of the clustering is also dependent on the relative magnitudes of the major and minor axis radii. If the major and minor axis have the same radii, then the ellipse is a circle and no clustering will be observed, irrespective of the value of `sqrt_distortion`. As the magnitude of the major axis radius increases relative to the minor axis radius, the strength of the clustering will increase. 
    
This effect is valuable when using an ellipse as a starting point to find a new level set in 2D (that may be an ellipse) that contains the starting ellipse. If we seek find the new level set by pushing out from the starting ellipse tangentially at each of the generated points then the points that diverge the fastest are located at the region of greatest curvature. If generated points are equally spaced with respect to arc length, then the new level set is likely to be well defined in regions that are approximately parallel to the major axis of the starting ellipse (and with length on a parallel line of similar length to the major axis, particularly for starting ellipses with a significantly larger major axis), but poorly defined in all other regions. The clustering effect is then valuable as it helps to better define the new level set.   

The function works by defining a new ellipse, `e_new` with minor axis radius equal to the supplied ellipse's minor axis radius (`e_new.b == e.b`) and major axis radius as a function of the supplied ellipse's, `e`, major and minor axis radii and the parameter `sqrt_distortion`. Namely, `e_new.a == e.b + sqrt_distortion^2 * (e.a - e.b)`. `e_new` is contained within `e` and can be varied between a circle and `e` using the parameter `sqrt_distortion` ∈ [0,1]. After defining `e_new` we determine the angle parameters `t_vector` that equally spaces `num_points` on `e_new` with respect to arc length and then find the points on `e` parameterised by the elements of `t_vector`.
"""
function generate_N_clustered_points(num_points::Int, 
                                    e::Ellipse; 
                                    start_point_shift::Float64=rand(),
                                    sqrt_distortion::Float64=0.0)

    0 < num_points || throw(DomainError("the number of points, num_points, to generate on the ellipse must be positive"))
    (0.0 <= sqrt_distortion && sqrt_distortion <= 1.0) || throw(DomainError("sqrt_distortion must be in the closed interval [0.0, 1.0]"))

    points = zeros(2, num_points)

    if sqrt_distortion == 1.0 # original ellipse
        return generate_N_equally_spaced_points(num_points, e,
                                                start_point_shift=start_point_shift)

    elseif sqrt_distortion == 0.0 # circle
        shift = start_point_shift/num_points
        angles = collect(LinRange(shift, (1+shift)*2.0*pi, num_points+1)[1:(end-1)])
    else
        d_a = (e.b + sqrt_distortion^2*(e.a-e.b))
        d_e = construct_ellipse(d_a, e.b)

        shift = start_point_shift/num_points

        lengths = collect(LinRange((shift)*d_e.circumference, (1+shift)*d_e.circumference, num_points+1))[1:(end-1)]
        angles = t_from_arclength_general.(lengths, Ref(d_e))
    end

    if e.x_radius < e.y_radius
        angles .= angles .+ 0.5*pi
    end
    for i in 1:num_points
        points[:,i] .= x_parametric_equation(angles[i], e), 
                y_parametric_equation(angles[i], e)
    end

    return points
end

"""
    generate_N_clustered_points(num_points::Int, x_radius::T, y_radius::T, α::T=0.0, Cx::T=0.0, Cy::T=0.0; start_point_shift::Float64=rand(), sqrt_distortion::Float64=0.0) where T<:Float64

An alternate way to call [`generate_N_clustered_points(num_points::Int, e::Ellipse; start_point_shift::Float64=rand(), sqrt_distortion::Float64=0.)`](@ref), by supplying the parameters of the ellipse to generate points on.

# Arguments
- `num_points`: a positive integer number of points to generate that are equally spaced on the ellipse. 
- `x_radius`: radius of the ellipse in the x axis (i.e. when the rotation, `α`, is zero).
- `y_radius`: radius of the ellipse in the y axis (i.e. when the rotation, `α`, is zero).
- `α`: an angle in radians (0 to 2π) that the ellipse has been rotated by. A positive value represents an anti-clockwise rotation. Default is `0.0`.
- `Cx`: the x coordinate of the centre of the ellipse (the translation of the ellipse in the x axis). Default is `0.0`.
- `Cy`: the y coordinate of the centre of the ellipse (the translation of the ellipse in the y axis). Default is `0.0`.

# Keyword Arguments
- `start_point_shift`: a number ∈ [0,1]. Default is `rand()` (defined on [0,1]), meaning that, by default, every time this function is called a different set of points will be generated.
- `sqrt_distortion`: a number ∈ [0,1]. Default is `0.0`, meaning that, by default, this function will evenly space points on the the ellipse `e` with respect to the parameter `t`.
"""
function generate_N_clustered_points(num_points::Int, x_radius::T, y_radius::T, α::T=0.0, Cx::T=0.0, Cy::T=0.0; start_point_shift::Float64=rand(), sqrt_distortion::Float64=0.0) where T<:Float64

    e = construct_ellipse(x_radius, y_radius, α, Cx, Cy)
    return generate_N_clustered_points(num_points, e, start_point_shift=start_point_shift,
                                        sqrt_distortion=sqrt_distortion)
end

"""
    generate_N_clustered_points(num_points::Int, Γ::Matrix{Float64}, θmle::Vector{Float64}, ind1::Int, ind2::Int; confidence_level::Float64=0.01, start_point_shift::Float64=rand(), sqrt_distortion::Float64=0.0)

An alternate way to call [`generate_N_clustered_points(num_points::Int, e::Ellipse; start_point_shift::Float64=rand(), sqrt_distortion::Float64=0.)`](@ref), by supplying a square matrix Γ, the inverse of the Hessian of a log-likelihood function at its maximum likelihood estimate, indexes of the two variables of interest and the confidence level that represent a 2D ellipse approximation of the log-likelihood function.

# Arguments
- `num_points`: a positive integer number of points to generate that are equally spaced on the ellipse. 
- `Γ`: a square matrix (2D) which is the inverse of the Hessian of a log-likelihood function at its maximum likelihood estimate.
- `θmle`: the maximum likelihood estimate for the parameters.
- `ind1`: index of the first parameter of interest (corresponds to the row and column index of `Γ`)
- `ind2`: index of the second parameter of interest (corresponds to the row and column index of `Γ`).

# Keyword Arguments
- `confidence_level`: the confidence level ∈ [0.0,1.0] at which the ellipse approximation is constructed. Default is `0.01`.
- `start_point_shift`: a number ∈ [0,1]. Default is `rand()` (defined on [0,1]), meaning that, by default, every time this function is called a different set of points will be generated.
- `sqrt_distortion`: a number ∈ [0,1]. Default is `0.0`, meaning that, by default, this function will evenly space points on the the ellipse `e` with respect to the parameter `t`.
"""
function generate_N_clustered_points(num_points::Int, Γ::Matrix{Float64}, θmle::Vector{Float64}, ind1::Int, ind2::Int; 
    confidence_level::Float64=0.01, start_point_shift::Float64=rand(), sqrt_distortion::Float64=0.0)

    _, _, x_radius, y_radius, α = calculate_ellipse_parameters(Γ, ind1, ind2, confidence_level)
    return generate_N_clustered_points(num_points, x_radius, y_radius, α, θmle[ind1], θmle[ind2],
                                    start_point_shift=start_point_shift, sqrt_distortion=sqrt_distortion)
end