using EllipseSampling
using Test
using LinearAlgebra
import Distributions

function isapprox_ellipsesampling(vec1::Vector{<:Float64}, vec2::Vector{<:Float64})
    return isapprox(vec1, vec2, atol=1e-14)
end

function isapprox_ellipsesampling(x1::Float64, x2::Float64)
    return isapprox(x1, x2, atol=1e-14)
end

@testset "EllipseSampling.jl" begin

    @testset "NEquallySpacedPointsTest" begin
        e = construct_ellipse(2.0, 1.0)
        correct_output = [2.0 -2.0; 0.0 0.0]
        num_points=2
        output = generate_N_equally_spaced_points(num_points, e, start_point_shift=0.0)
        for i in 1:num_points
            @test isapprox_ellipsesampling(correct_output[:, i], output[:, i]) 
        end

        e = construct_ellipse(1.0, 1.0)
        correct_output = [1.0 -0.5 -0.5; 0.0 sqrt(3.0)/2.0 -sqrt(3.0)/2.0]
        num_points=3
        output = generate_N_equally_spaced_points(num_points, e, start_point_shift=0.0)
        for i in 1:num_points
            @test isapprox_ellipsesampling(correct_output[:, i], output[:, i]) 
        end

        correct_output = [1.0 0.0 -1.0 0.0; 0.0 1.0 0.0 -1.0]
        num_points=4
        output = generate_N_equally_spaced_points(num_points, e, start_point_shift=0.0)
        for i in 1:num_points
            @test isapprox_ellipsesampling(correct_output[:, i], output[:, i]) 
        end
    end
    
    @testset "NClusteredPointsTest" begin
        e = construct_ellipse(2.0, 1.0)
        correct_output = [2.0 -2.0; 0.0 0.0]
        num_points=2
        output = generate_N_clustered_points(num_points, e, start_point_shift=0.0, sqrt_distortion=0.0)
        for i in 1:num_points
            @test isapprox_ellipsesampling(correct_output[:, i], output[:, i]) 
        end

        output = generate_N_clustered_points(num_points, e, start_point_shift=0.0, sqrt_distortion=0.5)
        for i in 1:num_points
            @test isapprox_ellipsesampling(correct_output[:, i], output[:, i]) 
        end

        output = generate_N_clustered_points(num_points, e, start_point_shift=0.0, sqrt_distortion=1.0)
        for i in 1:num_points
            @test isapprox_ellipsesampling(correct_output[:, i], output[:, i]) 
        end

        e = construct_ellipse(1.0, 1.0)
        correct_output = [1.0 -0.5 -0.5; 0.0 sqrt(3.0)/2.0 -sqrt(3.0)/2.0]
        num_points=3
        output = generate_N_clustered_points(num_points, e, start_point_shift=0.0, sqrt_distortion=0.1)
        for i in 1:num_points
            @test isapprox_ellipsesampling(correct_output[:, i], output[:, i]) 
        end

        correct_output = [1.0 0.0 -1.0 0.0; 0.0 1.0 0.0 -1.0]
        num_points=4
        output = generate_N_clustered_points(num_points, e, start_point_shift=0.0, sqrt_distortion=0.0)
        for i in 1:num_points
            @test isapprox_ellipsesampling(correct_output[:, i], output[:, i]) 
        end

        output = generate_N_clustered_points(num_points, e, start_point_shift=0.0, sqrt_distortion=0.1)
        for i in 1:num_points
            @test isapprox_ellipsesampling(correct_output[:, i], output[:, i]) 
        end
    end

    # Generate point on perimeter tests
    @testset "PointOnPerimeterTest" begin
        e = construct_ellipse(1.0, 1.0)
        @test isapprox_ellipsesampling(generate_perimeter_point(0.0, e), [1.0, 0.0])
        @test isapprox_ellipsesampling(generate_perimeter_point(0.25, e), [0.0, 1.0])

        e = construct_ellipse(2.0, 1.0)
        @test isapprox_ellipsesampling(generate_perimeter_point(0.0, e), [2.0, 0.0])
        @test isapprox_ellipsesampling(generate_perimeter_point(0.25, e), [0.0, 1.0])

        e = construct_ellipse(1.0, 2.0)
        @test isapprox_ellipsesampling(generate_perimeter_point(0.0, e), [0.0, 2.0])
        @test isapprox_ellipsesampling(generate_perimeter_point(0.25, e), [-1.0, 0.0])

        e = construct_ellipse(1.0, 1.0, pi/2)
        @test isapprox_ellipsesampling(generate_perimeter_point(0.0, e), [0.0, 1.0])
        @test isapprox_ellipsesampling(generate_perimeter_point(0.25, e), [-1.0, 0.0])

        @test isapprox_ellipsesampling(generate_perimeter_point(0.0, 1.0, 1.0, pi/2), [0.0, 1.0])
        @test isapprox_ellipsesampling(generate_perimeter_point(0.25, 1.0, 1.0, pi/2), [-1.0, 0.0])
    end

    @testset "ParametricEquationsTest" begin
        x_radius=1.0; y_radius=1.0; α=0.0; Cx=1.0
        @test isapprox_ellipsesampling(x_parametric_equation(0.0, x_radius, y_radius, α, Cx), 2.0)
        @test isapprox_ellipsesampling(x_parametric_equation(pi/2.0, x_radius, y_radius, α, Cx), 1.0)
        @test isapprox_ellipsesampling(x_parametric_equation(pi*1.0, x_radius, y_radius, α, Cx), 0.0)

        Cy=1.0
        @test isapprox_ellipsesampling(y_parametric_equation(0.0, x_radius, y_radius, α, Cy), 1.0)
        @test isapprox_ellipsesampling(y_parametric_equation(pi/2.0, x_radius, y_radius, α, Cy), 2.0)
        @test isapprox_ellipsesampling(y_parametric_equation(pi*1.0, x_radius, y_radius, α, Cy), 1.0)
    end

    @testset "TFromArcLengthTest" begin
        a=1.0; b=1.0; x_radius=1.0; y_radius=1.0
        @test isapprox_ellipsesampling(t_from_arclength_general(pi/2.0, a, b, x_radius, y_radius), pi/2.0)
        @test isapprox_ellipsesampling(t_from_arclength_general(pi*1.0, a, b, x_radius, y_radius), pi*1.0)
    end

    @testset "CalculateEllipseParametersTest" begin
        Γ = [7.7862e-6 -0.0506896 -0.0141446; -0.00506896 20.2146 6.61578; -0.0141446 6.61578 30.222]
        ind1, ind2 = 2, 3
        confidence_level = 0.01
        Hw = inv(Γ[[ind1, ind2], [ind1, ind2]]) .* 0.5 ./ (Distributions.quantile(Distributions.Chisq(2), confidence_level)*0.5)
        eigs = eigen(Hw)
        a_eig, b_eig = sqrt.(1.0 ./ eigs.values)
        eigvectors = eigs.vectors

        a, b, x_radius, y_radius, α  = EllipseSampling.calculate_ellipse_parameters(Γ, ind1, ind2, confidence_level)

        @test isapprox_ellipsesampling(a_eig, a)
        @test isapprox_ellipsesampling(b_eig, b)

        α_eig = atan(eigvectors[2,1], eigvectors[1,1])
        if α_eig < 0.0; α_eig += π end

        @test isapprox_ellipsesampling(α_eig, α)

        # For issue #30 when α → +/- 0.25 or +/- 1.25
        a, b = 2.0, 1.0 
        α = 0.25*π

        Hw11 = (cos(α)^2 / a^2 + sin(α)^2 / b^2)
        Hw22 = (sin(α)^2 / a^2 + cos(α)^2 / b^2)
        Hw12 = cos(α)*sin(α)*(1/a^2 - 1/b^2)
        Hw_norm = [Hw11 Hw12; Hw12 Hw22]

        confidence_level=0.95
        Hw = Hw_norm ./ (0.5 ./ (Distributions.quantile(Distributions.Chisq(2), confidence_level)*0.5))
        Γ = convert.(Float64, inv(BigFloat.(Hw, precision=64)))

        a_calc, b_calc, _, _, α_calc = EllipseSampling.calculate_ellipse_parameters(Γ, 1, 2, confidence_level)

        @test isapprox_ellipsesampling(a, a_calc)
        @test isapprox_ellipsesampling(b, b_calc)
        @test isapprox_ellipsesampling(α, α_calc)
    end

    @testset "CalculateEllipseParametersTest_HigherDof" begin
        Γ = [7.7862e-6 -0.0506896 -0.0141446; -0.00506896 20.2146 6.61578; -0.0141446 6.61578 30.222]
        ind1, ind2 = 2, 3
        confidence_level = 0.01
        dof=4
        Hw = inv(Γ[[ind1, ind2], [ind1, ind2]]) .* 0.5 ./ (Distributions.quantile(Distributions.Chisq(dof), confidence_level)*0.5)
        eigs = eigen(Hw)
        a_eig, b_eig = sqrt.(1.0 ./ eigs.values)
        eigvectors = eigs.vectors

        a, b, x_radius, y_radius, α  = EllipseSampling.calculate_ellipse_parameters(Γ, ind1, ind2, confidence_level, dof)

        @test isapprox_ellipsesampling(a_eig, a)
        @test isapprox_ellipsesampling(b_eig, b)

        α_eig = atan(eigvectors[2,1], eigvectors[1,1])
        if α_eig < 0.0; α_eig += π end

        @test isapprox_ellipsesampling(α_eig, α)

        # For issue #30 when α → +/- 0.25 or +/- 1.25
        a, b = 2.0, 1.0 
        α = 0.25*π

        Hw11 = (cos(α)^2 / a^2 + sin(α)^2 / b^2)
        Hw22 = (sin(α)^2 / a^2 + cos(α)^2 / b^2)
        Hw12 = cos(α)*sin(α)*(1/a^2 - 1/b^2)
        Hw_norm = [Hw11 Hw12; Hw12 Hw22]

        confidence_level=0.95
        Hw = Hw_norm ./ (0.5 ./ (Distributions.quantile(Distributions.Chisq(dof), confidence_level)*0.5))
        Γ = convert.(Float64, inv(BigFloat.(Hw, precision=64)))

        a_calc, b_calc, _, _, α_calc = EllipseSampling.calculate_ellipse_parameters(Γ, 1, 2, confidence_level, dof)

        @test isapprox_ellipsesampling(a, a_calc)
        @test isapprox_ellipsesampling(b, b_calc)
        @test isapprox_ellipsesampling(α, α_calc)
    end

end
