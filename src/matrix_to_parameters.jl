"""
    calculate_ellipse_parameters(Γ::Matrix{Float64}, ind1::Int, ind2::Int,
        confidence_level::Float64)

Given a square matrix Γ, the inverse of the Hessian of a log-likelihood function at its maximum likelihood estimate, indexes of the two variables of interest and the confidence level to construct a 2D ellipse approximation of the log-likelihood function, return the parameters of that ellipse.

# Arguments
- `Γ`: a square matrix which is the inverse of the Hessian of a log-likelihood function at its maximum likelihood estimate.
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
    Hw = inv(Γ[[ind1, ind2], [ind1, ind2]]) .* 0.5 ./ (Distributions.quantile(Distributions.Chisq(2), confidence_level) * 0.5)
    eigs = eigen(Hw)
    a_eig, b_eig = sqrt.(1.0 ./ eigs.values)

    # α_eig constructed such that it is the angle in radians from the x axis to the major axis (i.e. x_radius = a_eig) https://cookierobotics.com/007/
    if Hw[1,2]==0 
        α_eig = Hw[1,1] ≥ Hw[2,2] ? 0.0 : 0.5π
    else
        α_eig = atan(eigs.values[1]-Hw[1,1], Hw[1,2])

        if α_eig < 0.0
            α_eig += π # value in interval [0, 2pi]
        end
    end

    x_radius = a_eig; y_radius=b_eig
    return a_eig, b_eig, x_radius, y_radius, α_eig
end