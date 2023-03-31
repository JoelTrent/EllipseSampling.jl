"""
    calculate_ellipse_parameters(Γ::Matrix{Float64}, ind1::Int, ind2::Int,
        confidence_level::Float64)

Given a square matrix Γ, the inverse of the Hessian of a log-likelihood function at its maximum likelihood estimate, indexes of the two variables of interest and the confidence level to construct a 2D ellipse approximation of the log-likelihood function, return the parameters of that ellipse.

# Arguments
- `Γ`: A square matrix (2D) which is the inverse of the Hessian of a log-likelihood function at its maximum likelihood estimate.
- `ind1`: Index of the first parameter of interest (corresponds to the row and column index of `Γ`)
- `ind2`: Index of the second parameter of interest (corresponds to the row and column index of `Γ`).
- `confidence_level`: The confidence level ∈[0.0,1.0] at which the ellipse approximation is constructed.

# Details
The parameters of interest are `a` and `b`, the radius of the major and minor axis respectively, `x_radius` and `y_radius`, the radius of the ellipse in the x and y axis respectively (i.e. the radius when the rotation `α` is zero) and `α`, an angle between 0 and 2π radians that the ellipse has been rotated by. `a` is equal to the maximum of `x_radius` and `y_radius`, while b is equal to the minimum of `x_radius` and `y_radius`.

References and equation(s) to come.
"""
function calculate_ellipse_parameters(Γ::Matrix{Float64}, ind1::Int, ind2::Int,
    confidence_level::Float64)

    @assert (0.0 ≤ confidence_level && confidence_level ≤ 1.0) "The value of `confidence_level` is not between 0.0 and 1.0."
    @assert size(Γ)[1] == size(Γ)[2] "`Γ` must be a square matrix."
    @assert 0 < ind1 && ind1 ≤ size(Γ)[1] "`ind1` must be a valid row index in `Γ`." 
    @assert 0 < ind2 && ind2 ≤ size(Γ)[1] "`ind2` must be a valid row index in `Γ`." 

    Hw = inv(Γ[[ind1, ind2], [ind1, ind2]]) .* 0.5 ./ (Distributions.quantile(Distributions.Chisq(2), confidence_level)*0.5) # normalise Hw so that the RHS of the ellipse equation == 1

    α = atan(2*Hw[1,2]/(Hw[1,1]-Hw[2,2]))/2
    y_radius = sqrt( (cos(α)^4 - sin(α)^4) / (Hw[2,2]*(cos(α)^2) - Hw[1,1]*(sin(α)^2))  )
    x_radius = sqrt( (cos(α)^2) / (Hw[1,1] - (sin(α)^2)/y_radius^2))

    a = max(x_radius, y_radius)
    b = min(x_radius, y_radius)

    return a, b, x_radius, y_radius, α 
end