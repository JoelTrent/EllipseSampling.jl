using EllipseSampling
using Test

function equality_of_2D_coordinates(vec1::Vector{<:Float64}, vec2::Vector{<:Float64})
    return ((abs(vec1[1]-vec2[1]) < 1e-14) + (abs(vec1[2]-vec2[2]) < 1e-14)) == 2
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
end
