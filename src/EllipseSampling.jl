module EllipseSampling

import Roots, Elliptic, Distributions

export generateN_equally_spaced_points, generate_point_on_perimeter
    x_parametric_equation, y_parametric_equation, 
    t_from_arclength, t_from_arclength_general

# Julia version of functions from https://www.johndcook.com/blog/2022/11/02/ellipse-rng/
function E_inverse(em::T, z::T, m::T) where T<:Float64
    t = (z/em)*(pi/2)
    f(y) = E(y, m) - z
    r = Roots.find_zero(f, t, Roots.Order0())
    return r
end

function t_from_arclength(arc_len::T, a::T, b::T) where T<:Float64
    m = 1 - (b/a)^2
    em = Elliptic.E(m)
    t = 0.5*pi - E_inverse(em, em - arc_len/a, m)
    return t
end
#######################################################################################

function t_from_arclength_general(arc_len::T, a::T, b::T, x_radius::T, y_radius::T) where T<:Float64
    if x_radius < y_radius
        return t_from_length(arc_len, a, b) + 0.5*pi
    else
        return t_from_length(arc_len, a, b) 
    end
end

"""
    x_parametric_equation(t::T, x_radius::T, y_radius::T, α::T=0.0, Cx::T=0.0) where T<:Float64

Implements the parametric equation for variable x of a translated and rotated ellipse, ``x(t)``, where t is an angle between 0 and 2π radians.
    
# Arguments
- `t`: an angle between 0 and 2π radians that defines the x location on the ellipse.
- `x_radius`: radius of the ellipse in the x axis (i.e. when the rotation, `α`, is zero).
- `y_radius`: radius of the ellipse in the y axis (i.e. when the rotation, `α`, is zero).
- `α`: an angle between 0 and 2π radians that the ellipse has been rotated by. Default is `0.0`.
- `Cx`: the x coordinate of the centre of the ellipse (the translation of the ellipse in the x axis). Default is `0.0`.
"""
function x_parametric_equation(t::T, x_radius::T, y_radius::T, α::T=0.0, Cx::T=0.0) where T<:Float64
    return x_radius*(cos(t)*cos(α)) - y_radius*(sin(t)*sin(α)) + Cx
end

"""
    y_parametric_equation(t::T, x_radius::T, y_radius::T, α::T=0.0, Cx::T=0.0) where T<:Float64

Implements the parametric equation for variable y of a translated and rotated ellipse, ``y(t)``, where `t` is an angle between 0 and 2π radians.
    
# Arguments
- `t`: an angle between 0 and 2π radians that defines the y location on the ellipse.
- `x_radius`: radius of the ellipse in the x axis (i.e. when the rotation, `α`, is zero).
- `y_radius`: radius of the ellipse in the y axis (i.e. when the rotation, `α`, is zero).
- `α`: an angle between 0 and 2π radians that the ellipse has been rotated by. Default is `0.0`.
- `Cy`: the y coordinate of the centre of the ellipse (the translation of the ellipse in the y axis). Default is `0.0`.
"""
function y_parametric_equation(t::T, x_radius::T, y_radius::T, α::T=0.0, Cy::T=0.0) where T<:Float64
    return x_radius*(cos(t)*sin(α)) + y_radius*(sin(t)*cos(α)) + Cy
end

"""
    calculate_ellipse_parameters(Γ::Matrix{Float64}, ind1::Int, ind2::Int,
        confidence_level::Float64)

Given a square matrix Γ, the inverse of the Hessian of a log-likelihood function at its maximum likelihood estimate, indexes of the two variables of interest and the confidence level to construct a 2D ellipse approximation of the log-likelihood function, return the parameters of that ellipse; `a` and `b`, the radius of the major and minor axis respectively, `x_radius` and `y_radius`, the radius of the ellipse in the x and y axis respectively (i.e. the radius when the rotation `α` is zero) and `α`, an angle between 0 and 2π radians that the ellipse has been rotated by. `a` is equal to the maximum of `x_radius` and `y_radius`, while b is equal to the minimum of `x_radius` and `y_radius`.

References and equation(s) to come.

# Arguments
- `Γ`: A square matrix (2D) which is the inverse of the Hessian of a log-likelihood function at its maximum likelihood estimate.
- `ind1`: Index of the first parameter of interest (corresponds to the row and column index of `Γ`)
- `ind2`: Index of the second parameter of interest (corresponds to the row and column index of `Γ`).
- `confidence_level`: The confidence level ∈[0.0,1.0] at which the ellipse approximation is constructed.
"""
function calculate_ellipse_parameters(Γ::Matrix{Float64}, ind1::Int, ind2::Int,
    confidence_level::Float64)

    @assert (0.0 ≤ confidence_level && confidence_level ≤ 1.0) "The value of `confidence_level` is not between 0.0 and 1.0."
    @assert size(Γ)[1] == size(Γ)[2] "`Γ` must be a square matrix."
    @assert 0 < ind1 && ind1 ≤ size(Γ)[1] "`ind1` must be a valid row index in `Γ`." 
    @assert 0 < ind2 && ind2 ≤ size(Γ)[1] "`ind2` must be a valid row index in `Γ`." 

    Hw = inv(Γ[[ind1, ind2], [ind1, ind2]]) .* 0.5 ./ (quantile(Distributions.Chisq(2), confidence_level)*0.5) # normalise Hw so that the RHS of the ellipse equation == 1

    α = atan(2*Hw[1,2]/(Hw[1,1]-Hw[2,2]))/2
    y_radius = sqrt( (cos(α)^4 - sin(α)^4) / (Hw[2,2]*(cos(α)^2) - Hw[1,1]*(sin(α)^2))  )
    x_radius = sqrt( (cos(α)^2) / (Hw[1,1] - (sin(α)^2)/y_radius^2))

    a = max(x_radius, y_radius)
    b = min(x_radius, y_radius)

    return a, b, x_radius, y_radius, α 
end

function perimeter_length(a::Float64, b::Float64); return Elliptic.E(1 - (b/a)^2) * 4 * a end

function assert_parameters_are_valid(a::T, b::T, x_radius::T, y_radius::T)
    
    @assert (a ≥ b) "The radius of the major axis, `a`, must be greater than or equal to the radius of the minor axis, `b`. I.e. `a=max(x_radius, y_radius), b=min(x_radius, y_radius)`."
    @assert (x_radius > 0) "The `x_radius` must be strictly positive `(>0)`."
    @assert (y_radius > 0) "The `y_radius` must be strictly positive `(>0)`."
    return nothing
end

# norm_distance_on_perimeter ∈ [0,1]
function generate_point_on_perimeter(norm_distance_on_perimeter::T, a::T, b::T, x_radius::T, y_radius::T, α::T=0.0, Cx::T=0.0, Cy::T=0.0) where T<:Float64

    test_parameters_are_valid(a, b, x_radius, y_radius)
    @assert (0.0 ≤ start_point_shift && start_point_shift ≤ 1.0) "The value of `norm_distance_on_perimeter` is not between 0.0 and 1.0."
    
    point = zeros(2)
    perimeter_len = perimeter_length(a,b)
    angle = t_from_arclength_general(perimeter_len*norm_distance_on_perimeter, a, b, x_radius, y_radius)

    point[:] .= x_parametric_equation(angle, x_radius, y_radius, α, Cx), 
                    y_parametric_equation(angle, x_radius, y_radius, α, Cy)

    return point
end

# start_point_shift ∈ [0,1] (random by default)
function generateN_equally_spaced_points(num_points::Int, a::T, b::T, x_radius::T, y_radius::T, α::T=0.0, Cx::T=0.0, Cy::T=0.0; 
    start_point_shift::Float64=rand()) where T<:Float64

    assert_parameters_are_valid(a, b, x_radius, y_radius)
    @assert (0.0 ≤ start_point_shift && start_point_shift ≤ 1.0) "The value of `start_point_shift` is not between 0.0 and 1.0."

    points = zeros(2,num_points)
    perimeter_len = perimeter_length(a,b)

    shift = start_point_shift/num_points

    lengths = collect(LinRange((shift)*perimeter_len, (1+shift)*perimeter_len, num_points+1))[1:end-1]
    angles = t_from_arclength_general.(lengths, a, b, x_radius, y_radius)
    
    for i in 1:num_points
        points[:,i] .= x_parametric_equation(angles[i], x_radius, y_radius, α, Cx), 
            y_parametric_equation(angles[i], x_radius, y_radius, α, Cy)
    end

    return points
end

# start_point_shift ∈ [0,1] (random by default)
function generateN_equally_spaced_points(num_points::Int, Γ::Matrix{Float64}, θmle::Vector{Float64}, ind1::Int, ind2::Int; confidence_level::Float64=0.01,             
    start_point_shift::Float64=rand())

    a, b, x_radius, y_radius, α = calculate_ellipse_parameters(Γ, ind1, ind2, confidence_level)
    
    return generateN_equally_spaced_points(num_points, a, b, x_radius, y_radius, α, θmle[ind1], θmle[ind2],
                                    start_point_shift=start_point_shift)
end

end
