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
The parameters of interest are `a` and `b`, the radius of the major and minor axis respectively, `x_radius` and `y_radius`, the radius of the ellipse in the x and y axis respectively (i.e. the radius when the rotation `α` is zero) and `α`, an angle between 0 and π radians that the major axis of the ellipse has been rotated by from the positive x axis. `a` is equal to the maximum of `x_radius` and `y_radius`, while `b` is equal to the minimum of `x_radius` and `y_radius`.

## Observed Fisher Information Matrix and Approximation of the Log-Likelihood Function

We can approximate a log-likelihood function from multiple parameters, ``\\theta``, by considering the observed Fisher information matrix (FIM). The observed FIM is
a quadratic approximation of the curvature of the log-likelihood function at the maximum likelihood estimate, ``\\hat{\\theta}``. The observed FIM is the matrix of second derivatives (the Hessian) of the log-likelihood function evaluated at the MLE with elements [pawitanall2001](@cite): 
```math
H_{jk}(\\hat{\\theta}) \\equiv -\\frac{\\partial^2}{\\partial \\theta_j \\partial \\theta_k} \\ell (\\hat{\\theta} \\, ; \\, y_{1:I}^{\\textrm{o}} ).
```

This then allows us to define the following approximation of the normalised log-likelihood function
using a second-order Taylor expansion at the MLE [pawitanall2001](@cite):
```math
    \\hat{\\ell} (\\theta \\, ; \\, y_{1:I}^{\\textrm{o}} ) \\approx -\\frac{1}{2} (\\theta-\\hat{\\theta})' H(\\hat{\\theta}) (\\theta-\\hat{\\theta}).
```

Similarly, for two parameters, ``\\psi``, from a larger number of parameters, we first invert the observed FIM, ``\\Gamma(\\hat{\\theta}) = H^{-1}(\\hat{\\theta})``, and then select just the rows and columns relating to the parameters of interest, before again inverting the matrix:
```math
    \\hat{\\ell}_p (\\psi \\, ; \\, y_{1:I}^{\\textrm{o}} )  \\approx -\\frac{1}{2} (\\psi-\\hat{\\psi})' ([e_j, e_k]' \\, \\Gamma(\\hat{\\theta}) \\, [e_j, e_k])^{-1} (\\psi-\\hat{\\psi}), \\hspace{0.2cm} \\theta_j \\cup \\theta_k = \\psi, 
```
where ``e_j`` and ``e_k`` are the ``j``th and ``k``th canonical vectors of ``\\mathbb{R}^{|\\theta|}``.

## Obtaining Ellipse parameters

By normalising our log-likelihood approximation equation for two parameters by our target confidence threshold of interest ``\\ell_c`` (at `confidence_level`) so that one side of the equation is equal to 1 we obtain the equation of an ellipse [friendlyelliptical2013](@cite):
```math
1 = -\\frac{1}{2\\ell_c} (\\psi-\\hat{\\psi})' ([e_j, e_k]' \\, \\Gamma(\\hat{\\theta}) \\, [e_j, e_k])^{-1} (\\psi-\\hat{\\psi}) = (\\psi-\\hat{\\psi})' \\mathcal{C} (\\psi-\\hat{\\psi}),
```

The major and minor axis radii, `a` and `b` respectively, can then be evaluated by considering the inverse of the square roots of the eigenvalues of ``\\mathcal{C}`` (ordered from largest to smallest) [friendlyelliptical2013](@cite). To determine the rotation, `α` of the major axis of the ellipse from the positive ``x`` axis we calculate the inverse tangent of the division of the ``y`` and ``x`` components of the eigenvector corresponding to the largest eigenvalue [friendlyelliptical2013](@cite). 
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
    a_eig, b_eig = 1.0 ./ sqrt.(eigs.values)

    # α_eig constructed such that it is the angle in radians from the x axis to the major axis, 
    # i.e. angle between the x axis and the largest eigenvector
    α_eig = atan(eigs.vectors[2,1], eigs.vectors[1,1]) # value in interval [-pi, pi]
    if α_eig < 0.0; α_eig += π end # value in interval [0, pi]

    x_radius = a_eig; y_radius=b_eig
    return a_eig, b_eig, x_radius, y_radius, α_eig
end