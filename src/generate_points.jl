# Julia version of functions from https://www.johndcook.com/blog/2022/11/02/ellipse-rng/
function E_inverse(em::T, z::T, m::T) where T<:Float64
    t = (z/em)*(pi/2)
    f(y) = E(y, m) - z
    r = Roots.find_zero(f, t, Roots.Order0())
    return r
end

function t_from_arclength(arc_len::T, a::T, b::T) where T<:Float64
    m = eccentricity_squared(a,b)
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

# norm_distance_on_perimeter ∈ [0,1]
function generate_point_on_perimeter(norm_distance_on_perimeter::T, a::T, b::T, x_radius::T, y_radius::T, α::T=0.0, Cx::T=0.0, Cy::T=0.0) where T<:Float64

    test_parameters_are_valid(a, b, x_radius, y_radius)
    @assert (0.0 ≤ start_point_shift && start_point_shift ≤ 1.0) "The value of `norm_distance_on_perimeter` is not between 0.0 and 1.0."
    
    point = zeros(2)
    perimeter_len = circumference(a,b)
    angle = t_from_arclength_general(perimeter_len*norm_distance_on_perimeter, a, b, x_radius, y_radius)

    point[:] .= x_parametric_equation(angle, x_radius, y_radius, α, Cx), 
                    y_parametric_equation(angle, x_radius, y_radius, α, Cy)

    return point
end

function generate_point_on_perimeter(norm_distance_on_perimeter::T, Γ::Matrix{Float64}, θmle::Vector{Float64}, ind1::Int, ind2::Int; confidence_level::Float64=0.01) where T<:Float64

    a, b, x_radius, y_radius, α = calculate_ellipse_parameters(Γ, ind1, ind2, confidence_level)

    return generate_point_on_perimeter(norm_distance_on_perimeter, a, b, x_radius, y_radius, α, θmle[ind1], θmle[ind2])
end

# start_point_shift ∈ [0,1] (random by default)
function generateN_equally_spaced_points(num_points::Int, a::T, b::T, x_radius::T, y_radius::T, α::T=0.0, Cx::T=0.0, Cy::T=0.0; 
    start_point_shift::Float64=rand()) where T<:Float64

    assert_parameters_are_valid(a, b, x_radius, y_radius)
    @assert (0.0 ≤ start_point_shift && start_point_shift ≤ 1.0) "The value of `start_point_shift` is not between 0.0 and 1.0."
    @assert 0 < num_points "The number of points, `num_points`, to generate on the ellipse must be positive."

    points = zeros(2,num_points)
    perimeter_len = circumference(a,b)

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