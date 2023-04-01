using EllipseSampling
using Test

function equality_of_2D_coordinates(vec1::Vector{<:Float64}, vec2::Vector{<:Float64})
    return ((abs(vec1[1]-vec2[1]) < 1e-14) + (abs(vec1[2]-vec2[2]) < 1e-14)) == 2
end

function equality_of_1D_coordinates(x1::Float64, x2::Float64)
    return abs(x1-x2) < 1e-14
end

@testset "EllipseSampling.jl" begin

    @testset "NEquallySpacedPointsTest" begin
        e = construct_ellipse(2.0, 1.0)
        correct_output = [2.0 -2.0; 0.0 0.0]
        num_points=2
        output = generateN_equally_spaced_points(num_points, e, start_point_shift=0.0)
        for i in 1:num_points
            @test equality_of_2D_coordinates(correct_output[:, i], output[:, i]) 
        end

        e = construct_ellipse(1.0, 1.0)
        correct_output = [1.0 -0.5 -0.5; 0.0 sqrt(3.0)/2.0 -sqrt(3.0)/2.0]
        num_points=3
        output = generateN_equally_spaced_points(num_points, e, start_point_shift=0.0)
        for i in 1:num_points
            @test equality_of_2D_coordinates(correct_output[:, i], output[:, i]) 
        end

        correct_output = [1.0 0.0 -1.0 0.0; 0.0 1.0 0.0 -1.0]
        num_points=4
        output = generateN_equally_spaced_points(num_points, e, start_point_shift=0.0)
        for i in 1:num_points
            @test equality_of_2D_coordinates(correct_output[:, i], output[:, i]) 
        end
    end
    
    # Generate point on perimeter tests
    @testset "PointOnPerimeterTest" begin
        e = construct_ellipse(1.0, 1.0)
        @test equality_of_2D_coordinates(generate_point_on_perimeter(0.0, e), [1.0, 0.0])
        @test equality_of_2D_coordinates(generate_point_on_perimeter(0.25, e), [0.0, 1.0])

        e = construct_ellipse(2.0, 1.0)
        @test equality_of_2D_coordinates(generate_point_on_perimeter(0.0, e), [2.0, 0.0])
        @test equality_of_2D_coordinates(generate_point_on_perimeter(0.25, e), [0.0, 1.0])

        e = construct_ellipse(1.0, 2.0)
        @test equality_of_2D_coordinates(generate_point_on_perimeter(0.0, e), [0.0, 2.0])
        @test equality_of_2D_coordinates(generate_point_on_perimeter(0.25, e), [-1.0, 0.0])

        e = construct_ellipse(1.0, 1.0, pi/2)
        @test equality_of_2D_coordinates(generate_point_on_perimeter(0.0, e), [0.0, 1.0])
        @test equality_of_2D_coordinates(generate_point_on_perimeter(0.25, e), [-1.0, 0.0])

        @test equality_of_2D_coordinates(generate_point_on_perimeter(0.0, 1.0, 1.0, pi/2), [0.0, 1.0])
        @test equality_of_2D_coordinates(generate_point_on_perimeter(0.25, 1.0, 1.0, pi/2), [-1.0, 0.0])
    end

    @testset "ParametricEquationsTest" begin
        x_radius=1.0; y_radius=1.0; α=0.0; Cx=1.0
        @test equality_of_1D_coordinates(x_parametric_equation(0.0, x_radius, y_radius, α, Cx), 2.0)
        @test equality_of_1D_coordinates(x_parametric_equation(pi/2.0, x_radius, y_radius, α, Cx), 1.0)
        @test equality_of_1D_coordinates(x_parametric_equation(pi*1.0, x_radius, y_radius, α, Cx), 0.0)

        Cy=1.0
        @test equality_of_1D_coordinates(y_parametric_equation(0.0, x_radius, y_radius, α, Cy), 1.0)
        @test equality_of_1D_coordinates(y_parametric_equation(pi/2.0, x_radius, y_radius, α, Cy), 2.0)
        @test equality_of_1D_coordinates(y_parametric_equation(pi*1.0, x_radius, y_radius, α, Cy), 1.0)
    end

    @testset "TFromArcLengthTest" begin
        a=1.0; b=1.0; x_radius=1.0; y_radius=1.0
        @test equality_of_1D_coordinates(t_from_arclength_general(pi/2.0, a, b, x_radius, y_radius), pi/2.0)
        @test equality_of_1D_coordinates(t_from_arclength_general(pi*1.0, a, b, x_radius, y_radius), pi*1.0)
    end

    # @testset "CalculateEllipseParametersTest" begin
    #     calculate_ellipse_parameters() # or use the generateN_equally_spaced_points(Γ) version
    # end

end
