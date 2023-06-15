"""
    calculate_ellipse_parameters(Γ::Matrix{Float64}, ind1::Int, ind2::Int,
        confidence_level::Float64)

Given a square matrix Γ, the inverse of the Hessian of a log-likelihood function at its maximum likelihood estimate, indexes of the two variables of interest and the confidence level to construct a 2D ellipse approximation of the log-likelihood function, return the parameters of that ellipse.

# Arguments
- `Γ`: a square matrix (2D) which is the inverse of the Hessian of a log-likelihood function at its maximum likelihood estimate.
- `ind1`: index of the first parameter of interest (corresponds to the row and column index of `Γ`)
- `ind2`: index of the second parameter of interest (corresponds to the row and column index of `Γ`).
- `confidence_level`: the confidence level ∈ [0.0,1.0] at which the ellipse approximation is constructed.

# Details
The parameters of interest are `a` and `b`, the radius of the major and minor axis respectively, `x_radius` and `y_radius`, the radius of the ellipse in the x and y axis respectively (i.e. the radius when the rotation `α` is zero) and `α`, an angle between 0 and 2π radians that the ellipse has been rotated by. `a` is equal to the maximum of `x_radius` and `y_radius`, while b is equal to the minimum of `x_radius` and `y_radius`.

References and equation(s) to come.
"""
function calculate_ellipse_parameters(Γ::Matrix{Float64}, ind1::Int, ind2::Int,
    confidence_level::Float64)

    (0.0 ≤ confidence_level && confidence_level ≤ 1.0) || throw(DomainError(confidence_level, "confidence_level must be between 0.0 and 1.0."))
    size(Γ)[1] == size(Γ)[2] || throw(DimensionMismatch("Γ must be a square matrix."))
    (0 < ind1 && ind1 ≤ size(Γ)[1]) || throw(BoundsError("ind1 must be a valid row index in Γ.")) 
    (0 < ind2 && ind2 ≤ size(Γ)[1]) || throw(BoundsError("ind2 must be a valid row index in Γ."))

    # normalise Hw so that the RHS of the ellipse equation == 1
    Hw = inv(Γ[[ind1, ind2], [ind1, ind2]]) .* 0.5 ./ (Distributions.quantile(Distributions.Chisq(2), confidence_level)*0.5)

    α = atan(2*Hw[1,2]/(Hw[1,1]-Hw[2,2]))/2
    
    # if α close to +/-0.25pi, +/-1.25pi, then switch to BigFloat precision
    if isapprox(abs(rem(α/pi, 1)), 0.25, atol=1e-2)
        
        # convert values to BigFloat for enhanced precision - required for correct results when α → 0.25pi or 1.25pi.
        Hw = inv(BigFloat.(Γ[[ind1, ind2], [ind1, ind2]], RoundUp, precision=64)) .* 0.5 ./ (Distributions.quantile(Distributions.Chisq(2), confidence_level)*0.5) 
        
        α = atan(2*Hw[1,2]/(Hw[1,1]-Hw[2,2]))/2
        y_radius = sqrt( (cos(α)^4 - sin(α)^4) / (Hw[2,2]*(cos(α)^2) - Hw[1,1]*(sin(α)^2))  )
        x_radius = sqrt( (cos(α)^2) / (Hw[1,1] - (sin(α)^2)/y_radius^2))

        α, x_radius, y_radius = convert(Float64, α), convert(Float64, x_radius), convert(Float64, y_radius)
    else
        y_radius = sqrt( (cos(α)^4 - sin(α)^4) / (Hw[2,2]*(cos(α)^2) - Hw[1,1]*(sin(α)^2))  )
        x_radius = sqrt( (cos(α)^2) / (Hw[1,1] - (sin(α)^2)/y_radius^2))
    end

    a = max(x_radius, y_radius)
    b = min(x_radius, y_radius)

    # a, b and α could also be calculated by finding the sqrt of the reciprocal of the eigenvalues of the normalised inv(Γ)
    # using LinearAlgebra
    # sqrt.(1.0 ./ eigvals(inv(Γ[[ind1, ind2], [ind1, ind2]]) .* 0.5 ./ (Distributions.quantile(Distributions.Chisq(2), confidence_level)*0.5)))

    # this method for finding α assumes that a is the x axis
    # eigs = eigvecs(inv(Γ[[2,3], [2,3]]) .* 0.5 ./ (quantile(Chisq(2), 0.95)*0.5))
    # atan(eigs[2,1], eigs[1,1]) # if a is x axis
    # atan(eigs[2,1], eigs[1,1]) + pi/2 # if a is y axis

    return a, b, x_radius, y_radius, α 
end