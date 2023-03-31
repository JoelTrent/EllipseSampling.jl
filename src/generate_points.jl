# Julia version of functions from https://www.johndcook.com/blog/2022/11/02/ellipse-rng/
function E_inverse(em::T, z::T, m::T) where T<:Float64
    t = (z/em)*(pi/2)
    f(y) = Elliptic.E(y, m) - z
    r = Roots.find_zero(f, t, Roots.Order0())
    return r
end

function t_from_arclength(arc_len::T, a::T, b::T) where T<:Float64
    m = eccentricity_squared(a,b)
    em = Elliptic.E(m)
    t = 0.5*pi - E_inverse(em, em - arc_len/a, m)
    return t
end

function t_from_arclength(arc_len::Float64, e::Ellipse)
    em = Elliptic.E(e.m)
    t = 0.5*pi - E_inverse(em, em - arc_len/e.a, e.m)
    return t
end
#######################################################################################

function t_from_arclength_general(arc_len::T, a::T, b::T, x_radius::T, y_radius::T) where T<:Float64
    if x_radius < y_radius
        return t_from_arclength(arc_len, a, b) + 0.5*pi
    else
        return t_from_arclength(arc_len, a, b) 
    end
end

function t_from_arclength_general(arc_len::Float64, e::Ellipse)
    if x_radius < y_radius
        return t_from_arclength(arc_len, e) + 0.5*pi
    else
        return t_from_arclength(arc_len, e) 
    end
end

"""
    generate_point_on_perimeter(norm_distance_on_perimeter::Float64, e::Ellipse)

Generates a single point on an ellipse defined by the parameters contained within `e`, at distance `norm_distance_on_perimeter` ``\\times`` `e.circumference` around the circumference. 

# Arguments
- `norm_distance_on_perimeter`: A number ∈ [0,1] which represents the normalised distance on the perimeter of an ellipse. A value of `0.5` corresponds to a point halfway along the ellipse's perimeter, while a value of `0.7` corresponds to a point 70% along the ellipse's perimeter.
- `e`: A valid [`Ellipse`](@ref) struct which defines an ellipse.

# Details

This function can be easily used to generate uniform random samples from an ellipse by first sampling N points from a uniform distribution defined on [0,1] and then calling this function for each point 1:N.

For example using:
```julia
e = construct_ellipse(2,1)
N = 100
norm_samples = rand(N)
points = generate_point_on_perimeter.(samples, e)
```

Other distributions defined on [0,1] can be used to generate points on the ellipse's perimeter in a similar fashion.
"""
function generate_point_on_perimeter(norm_distance_on_perimeter::Float64, e::Ellipse)

    @assert (0.0 ≤ norm_distance_on_perimeter && norm_distance_on_perimeter ≤ 1.0) "The value of `norm_distance_on_perimeter` is not between 0.0 and 1.0."
    
    point = zeros(2)
    angle = t_from_arclength_general(e.circumference*norm_distance_on_perimeter, e)

    point[:] .= x_parametric_equation(angle, e), y_parametric_equation(angle, e)

    return point
end

"""
    generate_point_on_perimeter(norm_distance_on_perimeter::T, x_radius::T, y_radius::T, α::T=0.0, Cx::T=0.0, Cy::T=0.0) where T<:Float64

An alternative way to call [`generate_point_on_perimeter(norm_distance_on_perimeter::Float64, e::Ellipse)`](@ref), by supplying the parameters of the ellipse to generate a single point on.

# Arguments
- `norm_distance_on_perimeter`: A number ∈ [0,1] which represents the normalised distance on the perimeter of an ellipse. A value of 0.5 corresponds to a point halfway along the ellipse's perimeter, while a value of 0.7 corresponds to a point 70% along the ellipse's perimeter.
- `x_radius`: radius of the ellipse in the x axis (i.e. when the rotation, `α`, is zero).
- `y_radius`: radius of the ellipse in the y axis (i.e. when the rotation, `α`, is zero).
- `α`: an angle in radians (0 to 2π) that the ellipse has been rotated by. A positive value represents an anti-clockwise rotation. Default is `0.0`.
- `Cx`: the x coordinate of the centre of the ellipse (the translation of the ellipse in the x axis). Default is `0.0`.
- `Cy`: the y coordinate of the centre of the ellipse (the translation of the ellipse in the y axis). Default is `0.0`.
"""
function generate_point_on_perimeter(norm_distance_on_perimeter::T, x_radius::T, y_radius::T, α::T=0.0, Cx::T=0.0, Cy::T=0.0) where T<:Float64

    e = construct_ellipse(x_radius, y_radius, α, Cx, Cy)
    return generate_point_on_perimeter(norm_distance_on_perimeter, e, )
end

"""
    generate_point_on_perimeter(norm_distance_on_perimeter::T, Γ::Matrix{Float64}, θmle::Vector{Float64}, ind1::Int, ind2::Int; 
        confidence_level::Float64=0.01)

An alternative way to call [`generate_point_on_perimeter(norm_distance_on_perimeter::Float64, e::Ellipse)`](@ref), by supplying a square matrix Γ, the inverse of the Hessian of a log-likelihood function at its maximum likelihood estimate, indexes of the two variables of interest and the confidence level that represent a 2D ellipse approximation of the log-likelihood function.

# Arguments
- `norm_distance_on_perimeter`: A number ∈ [0,1] which represents the normalised distance on the perimeter of an ellipse. A value of 0.5 corresponds to a point halfway along the ellipse's perimeter, while a value of 0.7 corresponds to a point 70% along the ellipse's perimeter.
- `Γ`: A square matrix (2D) which is the inverse of the Hessian of a log-likelihood function at its maximum likelihood estimate.
- `θmle`: The maximum likelihood estimate for the parameters.
- `ind1`: Index of the first parameter of interest (corresponds to the row and column index of `Γ`)
- `ind2`: Index of the second parameter of interest (corresponds to the row and column index of `Γ`).

# Keyword Arguments
- `confidence_level`: The confidence level ∈[0.0,1.0] at which the ellipse approximation is constructed. Default is `0.01`.
"""
function generate_point_on_perimeter(norm_distance_on_perimeter::T, Γ::Matrix{Float64}, θmle::Vector{Float64}, ind1::Int, ind2::Int; 
    confidence_level::Float64=0.01)

    _, _, x_radius, y_radius, α = calculate_ellipse_parameters(Γ, ind1, ind2, confidence_level)
    return generate_point_on_perimeter(norm_distance_on_perimeter, x_radius, y_radius, α, θmle[ind1], θmle[ind2])
end

"""
    generateN_equally_spaced_points(num_points::Int, e::Ellipse; 
        start_point_shift::Float64=rand())

Generates `num_points` equally spaced points on an ellipse defined by the parameters contained within `e`. 

# Arguments
- `num_points`: A positive integer number of points to generate that are equally spaced on the ellipse. 
- `e`: A valid [`Ellipse`](@ref) struct which defines an ellipse.

# Keyword Arguments
- `start_point_shift`: A number ∈ [0,1]. Default is `rand()` (defined on [0,1]), meaning that, by default, every time this function is called a different set of points will be generated.

# Details

Equal spacing of points on the ellipse is with respect to arc length. The position of the first point placed on the ellipse can be shifted by `start_point_shift`, defined on [0.0,1.0], allowing the set of possible points generated for a given `num_points` to cover the full perimeter. This shift is normalised so that when `start_point_shift=1.0`, the position of the first point placed on the ellipse is equal to the position of the second point placed on the ellipse when `start_point_shift=0.0`.
"""
function generateN_equally_spaced_points(num_points::Int, e::Ellipse; 
    start_point_shift::Float64=rand())

    @assert (0.0 ≤ start_point_shift && start_point_shift ≤ 1.0) "The value of `start_point_shift` is not between 0.0 and 1.0."
    @assert 0 < num_points "The number of points, `num_points`, to generate on the ellipse must be positive."

    points = zeros(2,num_points)

    shift = start_point_shift/num_points

    lengths = collect(LinRange((shift)*e.circumference, (1+shift)*e.circumference, num_points+1))[1:end-1]
    angles = t_from_arclength_general.(lengths, e)
    
    for i in 1:num_points
        points[:,i] .= x_parametric_equation(angles[i], e), 
            y_parametric_equation(angles[i], e)
    end

    return points
end

"""
    generateN_equally_spaced_points(num_points::Int, x_radius::T, y_radius::T, α::T=0.0, Cx::T=0.0, Cy::T=0.0; 
        start_point_shift::Float64=rand()) where T<:Float64

An alternative way to call [`generateN_equally_spaced_points(num_points::Int, e::Ellipse; start_point_shift::Float64=rand())`](@ref), by supplying the parameters of the ellipse to generate points on.

# Arguments
- `num_points`: A positive integer number of points to generate that are equally spaced on the ellipse. 
- `x_radius`: radius of the ellipse in the x axis (i.e. when the rotation, `α`, is zero).
- `y_radius`: radius of the ellipse in the y axis (i.e. when the rotation, `α`, is zero).
- `α`: an angle in radians (0 to 2π) that the ellipse has been rotated by. A positive value represents an anti-clockwise rotation. Default is `0.0`.
- `Cx`: the x coordinate of the centre of the ellipse (the translation of the ellipse in the x axis). Default is `0.0`.
- `Cy`: the y coordinate of the centre of the ellipse (the translation of the ellipse in the y axis). Default is `0.0`.

# Keyword Arguments
- `start_point_shift`: A number ∈ [0,1]. Default is `rand()` (defined on [0,1]), meaning that, by default, every time this function is called a different set of points will be generated.
"""
function generateN_equally_spaced_points(num_points::Int, x_radius::T, y_radius::T, α::T=0.0, Cx::T=0.0, Cy::T=0.0; 
    start_point_shift::Float64=rand()) where T<:Float64

    e = construct_ellipse(x_radius, y_radius, α, Cx, Cy)
    return generateN_equally_spaced_points(num_points, e, start_point_shift=start_point_shift)
end

"""
    generateN_equally_spaced_points(num_points::Int, Γ::Matrix{Float64}, θmle::Vector{Float64}, ind1::Int, ind2::Int; 
        confidence_level::Float64=0.01, start_point_shift::Float64=rand())

An alternative way to call [`generateN_equally_spaced_points(num_points::Int, e::Ellipse; start_point_shift::Float64=rand())`](@ref), by supplying a square matrix Γ, the inverse of the Hessian of a log-likelihood function at its maximum likelihood estimate, indexes of the two variables of interest and the confidence level that represent a 2D ellipse approximation of the log-likelihood function.

# Arguments
- `num_points`: A positive integer number of points to generate that are equally spaced on the ellipse. 
- `Γ`: A square matrix (2D) which is the inverse of the Hessian of a log-likelihood function at its maximum likelihood estimate.
- `θmle`: The maximum likelihood estimate for the parameters.
- `ind1`: Index of the first parameter of interest (corresponds to the row and column index of `Γ`)
- `ind2`: Index of the second parameter of interest (corresponds to the row and column index of `Γ`).

# Keyword Arguments
- `confidence_level`: The confidence level ∈[0.0,1.0] at which the ellipse approximation is constructed. Default is `0.01`.
- `start_point_shift`: A number ∈ [0,1]. Default is `rand()` (defined on [0,1]), meaning that, by default, every time this function is called a different set of points will be generated.
"""
function generateN_equally_spaced_points(num_points::Int, Γ::Matrix{Float64}, θmle::Vector{Float64}, ind1::Int, ind2::Int; 
    confidence_level::Float64=0.01, start_point_shift::Float64=rand())

    _, _, x_radius, y_radius, α = calculate_ellipse_parameters(Γ, ind1, ind2, confidence_level)
    return generateN_equally_spaced_points(num_points, x_radius, y_radius, α, θmle[ind1], θmle[ind2],
                                    start_point_shift=start_point_shift)
end