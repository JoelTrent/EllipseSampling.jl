module EllipseSampling

import Roots, Elliptic, Distributions

export construct_ellipse,
    generateN_equally_spaced_points, generate_point_on_perimeter
    x_parametric_equation, y_parametric_equation, 
    t_from_arclength, t_from_arclength_general

include("mathematic_relationships.jl")
include("ellipse_struct.jl")
include("parametric_equations.jl")
include("generate_points.jl")
include("matrix_to_parameters.jl")

import SnoopPrecompile

SnoopPrecompile.@precompile_all_calls begin

    Γ = [7.7862e-6 -0.0506896 -0.0141446; -0.00506896 20.2146 6.61578; -0.0141446 6.61578 30.222]
    generateN_equally_spaced_points(1, Γ, [0.0, 0.0], 2, 3, confidence_level=0.01, start_point_shift=0.0)
    generate_point_on_perimeter(0.0, Γ, [0.0, 0.0], 2, 3, confidence_level=0.01)

    ellipse = construct_ellipse(2,1)
    # generate_point_on_perimeter(0.0, ellipse)
    # generateN_equally_spaced_points(1, ellipse)


end
end
