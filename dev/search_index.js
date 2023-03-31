var documenterSearchIndex = {"docs":
[{"location":"internal_library/","page":"Internal Library","title":"Internal Library","text":"Pages = [\"internal_libary.md\"]","category":"page"},{"location":"internal_library/#Internal-Library","page":"Internal Library","title":"Internal Library","text":"","category":"section"},{"location":"internal_library/","page":"Internal Library","title":"Internal Library","text":"Documentation for structs and functions not covered within User Interface documentation.","category":"page"},{"location":"internal_library/#Structs","page":"Internal Library","title":"Structs","text":"","category":"section"},{"location":"internal_library/","page":"Internal Library","title":"Internal Library","text":"EllipseSampling.Ellipse","category":"page"},{"location":"internal_library/#EllipseSampling.Ellipse","page":"Internal Library","title":"EllipseSampling.Ellipse","text":"Ellipse(x_radius::Float64,\n    y_radius::Float64,\n    α::Float64,\n    Cx::Float64,\n    Cy::Float64,\n    a::Float64,\n    b::Float64,\n    m::Float64,\n    circumference::Float64)\n\nContains the information required to define an ellipse which may have been rotated and translated. See construct_ellipse\n\nArguments\n\nx_radius: radius of the ellipse in the x axis (i.e. when the rotation, α, is zero).\ny_radius: radius of the ellipse in the y axis (i.e. when the rotation, α, is zero).\nα: an angle in radians (0 to 2π) that the ellipse has been rotated by. A positive value represents an anti-clockwise rotation.\nCx: the x coordinate of the centre of the ellipse (the translation of the ellipse in the x axis).\nCy: the y coordinate of the centre of the ellipse (the translation of the ellipse in the y axis). \na: the major radius of the ellipse.\nb: the minor radius of the ellipse.\nm: the eccentricity of the ellipse squared. See eccentricity_squared\ncircumference: the circumference of the ellipse.\n\n\n\n\n\n","category":"type"},{"location":"internal_library/#Functions","page":"Internal Library","title":"Functions","text":"","category":"section"},{"location":"internal_library/","page":"Internal Library","title":"Internal Library","text":"EllipseSampling.eccentricity_squared\nEllipseSampling.circumference\nEllipseSampling.assert_parameters_are_valid\nEllipseSampling.E_inverse\nEllipseSampling.calculate_ellipse_parameters","category":"page"},{"location":"internal_library/#EllipseSampling.eccentricity_squared","page":"Internal Library","title":"EllipseSampling.eccentricity_squared","text":"eccentricity_squared(a::T, b::T) where T<:Float64\n\nComputes m, the eccentricity of the ellipse squared, m=e^2, where m=1-big(fracbabig)^2.\n\nArguments\n\na: the major radius of the ellipse.\nb: the minor radius of the ellipse.\n\nDetails\n\nThis relationship between m and e is seen by considering the equation for e: e=fracca, where c^2=big(a^2-b^2big) and ab.\n\nReplacing c in the equation for e: \n\ne=fracsqrta^2-b^2a\n\nSubstituting the equation for e into m: \n\nm = bigg(fracsqrta^2-b^2abigg)^2 = fraca^2-b^2a^2 = 1 - fracb^2a^2 = 1 - bigg(fracbabigg)^2\n\n\n\n\n\n","category":"function"},{"location":"internal_library/#EllipseSampling.circumference","page":"Internal Library","title":"EllipseSampling.circumference","text":"circumference(a::Float64, b::Float64)\n\nCalculates the circumference of an ellipse using the elliptic integral of the second kind. Uses Elliptic.jl.\n\nArguments\n\na: the major radius of the ellipse.\nb: the minor radius of the ellipse.\n\n\n\n\n\n","category":"function"},{"location":"internal_library/#EllipseSampling.assert_parameters_are_valid","page":"Internal Library","title":"EllipseSampling.assert_parameters_are_valid","text":"assert_parameters_are_valid(a::T, b::T, x_radius::T, y_radius::T) where T<:Float64\n\nAsserts that the parameters relate to a valid ellipse. I.e. that a ≥ b and, x_radius and y_radius are positive. Note: a=max(x_radius, y_radius), b=min(x_radius, y_radius).\n\nArguments\n\na: the major radius of the ellipse.\nb: the minor radius of the ellipse.\nx_radius: radius of the ellipse in the x axis (i.e. when the rotation, α, is zero).\ny_radius: radius of the ellipse in the y axis (i.e. when the rotation, α, is zero).\n\n\n\n\n\n","category":"function"},{"location":"internal_library/#EllipseSampling.calculate_ellipse_parameters","page":"Internal Library","title":"EllipseSampling.calculate_ellipse_parameters","text":"calculate_ellipse_parameters(Γ::Matrix{Float64}, ind1::Int, ind2::Int,\n    confidence_level::Float64)\n\nGiven a square matrix Γ, the inverse of the Hessian of a log-likelihood function at its maximum likelihood estimate, indexes of the two variables of interest and the confidence level to construct a 2D ellipse approximation of the log-likelihood function, return the parameters of that ellipse.\n\nArguments\n\nΓ: A square matrix (2D) which is the inverse of the Hessian of a log-likelihood function at its maximum likelihood estimate.\nind1: Index of the first parameter of interest (corresponds to the row and column index of Γ)\nind2: Index of the second parameter of interest (corresponds to the row and column index of Γ).\nconfidence_level: The confidence level ∈[0.0,1.0] at which the ellipse approximation is constructed.\n\nDetails\n\nThe parameters of interest are a and b, the radius of the major and minor axis respectively, x_radius and y_radius, the radius of the ellipse in the x and y axis respectively (i.e. the radius when the rotation α is zero) and α, an angle between 0 and 2π radians that the ellipse has been rotated by. a is equal to the maximum of x_radius and y_radius, while b is equal to the minimum of x_radius and y_radius.\n\nReferences and equation(s) to come.\n\n\n\n\n\n","category":"function"},{"location":"user_interface/","page":"User Interface","title":"User Interface","text":"Pages = [\"user_interface.md\"]","category":"page"},{"location":"user_interface/#User-Interface","page":"User Interface","title":"User Interface","text":"","category":"section"},{"location":"user_interface/","page":"User Interface","title":"User Interface","text":"construct_ellipse\ngenerateN_equally_spaced_points\ngenerate_point_on_perimeter\nx_parametric_equation \ny_parametric_equation\nt_from_arclength\nt_from_arclength_general","category":"page"},{"location":"user_interface/#EllipseSampling.construct_ellipse","page":"User Interface","title":"EllipseSampling.construct_ellipse","text":"construct_ellipse(x_radius::T, y_radius::T, α::T=0.0, Cx::T=0.0, Cy::T=0.0) where T<:Float64\n\nConstructs a Ellipse struct which contains the information required to define an ellipse which may have been rotated and translated.\n\nArguments\n\nx_radius: radius of the ellipse in the x axis (i.e. when the rotation, α, is zero).\ny_radius: radius of the ellipse in the y axis (i.e. when the rotation, α, is zero).\nα: an angle in radians (0 to 2π) that the ellipse has been rotated by. A positive value represents an anti-clockwise rotation. Default is 0.0.\nCx: the x coordinate of the centre of the ellipse (the translation of the ellipse in the x axis). Default is 0.0.\nCy: the y coordinate of the centre of the ellipse (the translation of the ellipse in the y axis). Default is 0.0.\n\nDetails\n\nThe general equation for a rotated and translated ellipse is given by:\n\n1 = A(x-C_x)^2 + B(x-C_x)(y-C_y) + C(y-C_y)^2\n\nWhere:\n\nA = bigg(fraccos^2(α)a^2 + fracsin^2(α)b^2bigg)\n\nB = 2cos(α)sin(α)bigg(frac1a^2 - frac1b^2bigg)\n\nC = bigg(fracsin^2(α)a^2 + fraccos^2(α)b^2bigg)\n\n\n\n\n\n","category":"function"},{"location":"user_interface/#EllipseSampling.generateN_equally_spaced_points","page":"User Interface","title":"EllipseSampling.generateN_equally_spaced_points","text":"generateN_equally_spaced_points(num_points::Int, e::Ellipse; \n    start_point_shift::Float64=rand())\n\nGenerates num_points equally spaced points on an ellipse defined by the parameters contained within e. \n\nArguments\n\nnum_points: A positive integer number of points to generate that are equally spaced on the ellipse. \ne: A valid Ellipse struct which defines an ellipse.\n\nKeyword Arguments\n\nstart_point_shift: A number ∈ [0,1]. Default is rand() (defined on [0,1]), meaning that, by default, every time this function is called a different set of points will be generated.\n\nDetails\n\nEqual spacing of points on the ellipse is with respect to arc length. The position of the first point placed on the ellipse can be shifted by start_point_shift, defined on [0.0,1.0], allowing the set of possible points generated for a given num_points to cover the full perimeter. This shift is normalised so that when start_point_shift=1.0, the position of the first point placed on the ellipse is equal to the position of the second point placed on the ellipse when start_point_shift=0.0.\n\n\n\n\n\ngenerateN_equally_spaced_points(num_points::Int, x_radius::T, y_radius::T, α::T=0.0, Cx::T=0.0, Cy::T=0.0; \n    start_point_shift::Float64=rand()) where T<:Float64\n\nAn alternative way to call generateN_equally_spaced_points(num_points::Int, e::Ellipse; start_point_shift::Float64=rand()), by supplying the parameters of the ellipse to generate points on.\n\nArguments\n\nnum_points: A positive integer number of points to generate that are equally spaced on the ellipse. \nx_radius: radius of the ellipse in the x axis (i.e. when the rotation, α, is zero).\ny_radius: radius of the ellipse in the y axis (i.e. when the rotation, α, is zero).\nα: an angle in radians (0 to 2π) that the ellipse has been rotated by. A positive value represents an anti-clockwise rotation. Default is 0.0.\nCx: the x coordinate of the centre of the ellipse (the translation of the ellipse in the x axis). Default is 0.0.\nCy: the y coordinate of the centre of the ellipse (the translation of the ellipse in the y axis). Default is 0.0.\n\nKeyword Arguments\n\nstart_point_shift: A number ∈ [0,1]. Default is rand() (defined on [0,1]), meaning that, by default, every time this function is called a different set of points will be generated.\n\n\n\n\n\ngenerateN_equally_spaced_points(num_points::Int, Γ::Matrix{Float64}, θmle::Vector{Float64}, ind1::Int, ind2::Int; \n    confidence_level::Float64=0.01, start_point_shift::Float64=rand())\n\nAn alternative way to call generateN_equally_spaced_points(num_points::Int, e::Ellipse; start_point_shift::Float64=rand()), by supplying a square matrix Γ, the inverse of the Hessian of a log-likelihood function at its maximum likelihood estimate, indexes of the two variables of interest and the confidence level that represent a 2D ellipse approximation of the log-likelihood function.\n\nArguments\n\nnum_points: A positive integer number of points to generate that are equally spaced on the ellipse. \nΓ: A square matrix (2D) which is the inverse of the Hessian of a log-likelihood function at its maximum likelihood estimate.\nθmle: The maximum likelihood estimate for the parameters.\nind1: Index of the first parameter of interest (corresponds to the row and column index of Γ)\nind2: Index of the second parameter of interest (corresponds to the row and column index of Γ).\n\nKeyword Arguments\n\nconfidence_level: The confidence level ∈[0.0,1.0] at which the ellipse approximation is constructed. Default is 0.01.\nstart_point_shift: A number ∈ [0,1]. Default is rand() (defined on [0,1]), meaning that, by default, every time this function is called a different set of points will be generated.\n\n\n\n\n\n","category":"function"},{"location":"user_interface/#EllipseSampling.generate_point_on_perimeter","page":"User Interface","title":"EllipseSampling.generate_point_on_perimeter","text":"generate_point_on_perimeter(norm_distance_on_perimeter::Float64, e::Ellipse)\n\nGenerates a single point on an ellipse defined by the parameters contained within e, at distance norm_distance_on_perimeter times e.circumference around the circumference. \n\nArguments\n\nnorm_distance_on_perimeter: A number ∈ [0,1] which represents the normalised distance on the perimeter of an ellipse. A value of 0.5 corresponds to a point halfway along the ellipse's perimeter, while a value of 0.7 corresponds to a point 70% along the ellipse's perimeter.\ne: A valid Ellipse struct which defines an ellipse.\n\nDetails\n\nThis function can be easily used to generate uniform random samples from an ellipse by first sampling N points from a uniform distribution defined on [0,1] and then calling this function for each point 1:N.\n\nFor example using:\n\ne = construct_ellipse(2,1)\nN = 100\nnorm_samples = rand(N)\npoints = generate_point_on_perimeter.(samples, Ref(e))\n\nNote, here we wrap the ellipse struct e in Ref so that Julia does not try to broadcast over e as well.\n\nOther distributions defined on [0,1] can be used to generate points on the ellipse's perimeter in a similar fashion.\n\n\n\n\n\ngenerate_point_on_perimeter(norm_distance_on_perimeter::T, x_radius::T, y_radius::T, α::T=0.0, Cx::T=0.0, Cy::T=0.0) where T<:Float64\n\nAn alternative way to call generate_point_on_perimeter(norm_distance_on_perimeter::Float64, e::Ellipse), by supplying the parameters of the ellipse to generate a single point on.\n\nArguments\n\nnorm_distance_on_perimeter: A number ∈ [0,1] which represents the normalised distance on the perimeter of an ellipse. A value of 0.5 corresponds to a point halfway along the ellipse's perimeter, while a value of 0.7 corresponds to a point 70% along the ellipse's perimeter.\nx_radius: radius of the ellipse in the x axis (i.e. when the rotation, α, is zero).\ny_radius: radius of the ellipse in the y axis (i.e. when the rotation, α, is zero).\nα: an angle in radians (0 to 2π) that the ellipse has been rotated by. A positive value represents an anti-clockwise rotation. Default is 0.0.\nCx: the x coordinate of the centre of the ellipse (the translation of the ellipse in the x axis). Default is 0.0.\nCy: the y coordinate of the centre of the ellipse (the translation of the ellipse in the y axis). Default is 0.0.\n\n\n\n\n\ngenerate_point_on_perimeter(norm_distance_on_perimeter::Float64, Γ::Matrix{Float64}, θmle::Vector{Float64}, ind1::Int, ind2::Int; \n    confidence_level::Float64=0.01)\n\nAn alternative way to call generate_point_on_perimeter(norm_distance_on_perimeter::Float64, e::Ellipse), by supplying a square matrix Γ, the inverse of the Hessian of a log-likelihood function at its maximum likelihood estimate, indexes of the two variables of interest and the confidence level that represent a 2D ellipse approximation of the log-likelihood function.\n\nArguments\n\nnorm_distance_on_perimeter: A number ∈ [0,1] which represents the normalised distance on the perimeter of an ellipse. A value of 0.5 corresponds to a point halfway along the ellipse's perimeter, while a value of 0.7 corresponds to a point 70% along the ellipse's perimeter.\nΓ: A square matrix (2D) which is the inverse of the Hessian of a log-likelihood function at its maximum likelihood estimate.\nθmle: The maximum likelihood estimate for the parameters.\nind1: Index of the first parameter of interest (corresponds to the row and column index of Γ)\nind2: Index of the second parameter of interest (corresponds to the row and column index of Γ).\n\nKeyword Arguments\n\nconfidence_level: The confidence level ∈[0.0,1.0] at which the ellipse approximation is constructed. Default is 0.01.\n\n\n\n\n\n","category":"function"},{"location":"user_interface/#EllipseSampling.x_parametric_equation","page":"User Interface","title":"EllipseSampling.x_parametric_equation","text":"x_parametric_equation(t::T, x_radius::T, y_radius::T, α::T=0.0, Cx::T=0.0) where T<:Float64\n\nImplements the parametric equation for variable x of a translated and rotated ellipse, x(t), where t is an angle between 0 and 2π radians.\n\nArguments\n\nt: an angle between 0 and 2π radians that defines the x location on the ellipse.\nx_radius: radius of the ellipse in the x axis (i.e. when the rotation, α, is zero).\ny_radius: radius of the ellipse in the y axis (i.e. when the rotation, α, is zero).\nα: an angle in radians (0 to 2π) that the ellipse has been rotated by. A positive value represents an anti-clockwise rotation. Default is 0.0.\nCx: the x coordinate of the centre of the ellipse (the translation of the ellipse in the x axis). Default is 0.0.\n\n\n\n\n\nx_parametric_equation(t::T, e::Ellipse) where T<:Float64\n\nImplements the parametric equation for variable x of a translated and rotated ellipse, x(t), where t is an angle between 0 and 2π radians.\n\nArguments\n\nt: an angle between 0 and 2π radians that defines the x location on the ellipse.\ne: A valid Ellipse struct which defines an ellipse.\n\n\n\n\n\n","category":"function"},{"location":"user_interface/#EllipseSampling.y_parametric_equation","page":"User Interface","title":"EllipseSampling.y_parametric_equation","text":"y_parametric_equation(t::T, x_radius::T, y_radius::T, α::T=0.0, Cx::T=0.0) where T<:Float64\n\nImplements the parametric equation for variable y of a translated and rotated ellipse, y(t), where t is an angle between 0 and 2π radians.\n\nArguments\n\nt: an angle between 0 and 2π radians that defines the y location on the ellipse.\nx_radius: radius of the ellipse in the x axis (i.e. when the rotation, α, is zero).\ny_radius: radius of the ellipse in the y axis (i.e. when the rotation, α, is zero).\nα: an angle in radians (0 to 2π) that the ellipse has been rotated by. A positive value represents an anti-clockwise rotation. Default is 0.0.\nCy: the y coordinate of the centre of the ellipse (the translation of the ellipse in the y axis). Default is 0.0.\n\n\n\n\n\ny_parametric_equation(t::T, e::Ellipse) where T<:Float64\n\nImplements the parametric equation for variable y of a translated and rotated ellipse, y(t), where t is an angle between 0 and 2π radians.\n\nArguments\n\nt: an angle between 0 and 2π radians that defines the y location on the ellipse.\ne: A valid Ellipse struct which defines an ellipse.\n\n\n\n\n\n","category":"function"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = EllipseSampling","category":"page"},{"location":"#EllipseSampling.jl","page":"Home","title":"EllipseSampling.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"EllipseSampling.jl is a lightweight package for generating points on the boundary defined by an 2D ellipse built around Julia implementations of two functions by John D. Cook. It handles cases where the ellipse has been rotated and/or translated, and where either of the x and y axes is the major or minor axis. It provides a method (generateN_equally_spaced_points) to generate N equally-spaced points on the ellipse's boundary as well as a method (generate_point_on_perimeter) to sample uniformly on the boundary. Resultantly, any distribution defined on [0, 1] can be used to sample points. Calculation of the arc length and circumference of an ellipse uses Elliptic.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package's creation was motivated by the need to sample points on the boundary defined by an elliptical approximation of the log-likelihood function around the maximum likelihood estimate of a mechanistic model at a particular confidence level. Resultantly, it also provides a method (EllipseSampling.calculate_ellipse_parameters) to convert the matrix representation of this approximation into the parameters of the equivalent ellipse. Points can then be sampled from this ellipse.","category":"page"},{"location":"","page":"Home","title":"Home","text":"For a tutorial on how to use the package see Tutorial.","category":"page"},{"location":"","page":"Home","title":"Home","text":"A deeper dive into the user interface and internal library can be found in User Interface and Internal Library, respectively.","category":"page"},{"location":"#Getting-Started:-Installation-And-First-Steps","page":"Home","title":"Getting Started: Installation And First Steps","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install the package, use the following command inside the Julia REPL:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(url=\"https://github.com/JoelTrent/EllipseSampling.jl\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"To load the package, use the command:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using EllipseSampling","category":"page"},{"location":"tutorial/#Tutorial","page":"Tutorial","title":"Tutorial","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"TODO:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Information on how to call generate_point_on_perimeter.\nInformation on how to call generateN_equally_spaced_points.\nPlots of different points generated on an ellipse, as part of demonstrating how the above functions work.\nNote on, for ease of use, these generation functions can either be called directly with ellipse parameters such as a, b and α, or we can construct an ellipse struct first with those parameters, and use the ellipse struct when calling the generation functions.","category":"page"}]
}
