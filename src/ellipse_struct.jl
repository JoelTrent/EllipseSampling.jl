
"""
    assert_parameters_are_valid(a::T, b::T, x_radius::T, y_radius::T) where T<:Float64

Asserts that the parameters relate to a valid ellipse. I.e. that `a ≥ b` and, `x_radius` and `y_radius` are positive. Note: `a=max(x_radius, y_radius), b=min(x_radius, y_radius)`.

# Arguments
- `a`: the major radius of the ellipse.
- `b`: the minor radius of the ellipse.
- `x_radius`: radius of the ellipse in the x axis (i.e. when the rotation, `α`, is zero).
- `y_radius`: radius of the ellipse in the y axis (i.e. when the rotation, `α`, is zero).
"""
function assert_parameters_are_valid(a::T, b::T, x_radius::T, y_radius::T) where T<:Float64
    
    @assert (a ≥ b) "The radius of the major axis, `a`, must be greater than or equal to the radius of the minor axis, `b`. I.e. `a=max(x_radius, y_radius), b=min(x_radius, y_radius)`."
    @assert (x_radius > 0) "The `x_radius` must be strictly positive `(>0)`."
    @assert (y_radius > 0) "The `y_radius` must be strictly positive `(>0)`."
    return nothing
end

"""
    Ellipse(x_radius::Float64,
        y_radius::Float64,
        α::Float64,
        Cx::Float64,
        Cy::Float64,
        a::Float64,
        b::Float64,
        m::Float64,
        circumference::Float64)

Contains the information required to define an ellipse which may have been rotated and translated. See [`construct_ellipse`](@ref)

# Arguments
- `x_radius`: radius of the ellipse in the x axis (i.e. when the rotation, `α`, is zero).
- `y_radius`: radius of the ellipse in the y axis (i.e. when the rotation, `α`, is zero).
- `α`: an angle in radians (0 to 2π) that the ellipse has been rotated by. A positive value represents an anti-clockwise rotation.
- `Cx`: the x coordinate of the centre of the ellipse (the translation of the ellipse in the x axis).
- `Cy`: the y coordinate of the centre of the ellipse (the translation of the ellipse in the y axis). 
- `a`: the major radius of the ellipse.
- `b`: the minor radius of the ellipse.
- `m`: the eccentricity of the ellipse squared. See [`eccentricity_squared`](@ref)
- `circumference`: the circumference of the ellipse.
"""
struct Ellipse
    x_radius::Float64
    y_radius::Float64
    α::Float64
    Cx::Float64
    Cy::Float64
    a::Float64
    b::Float64
    m::Float64
    circumference::Float64
end

"""
    construct_ellipse(x_radius::T, y_radius::T, α::T=0.0, Cx::T=0.0, Cy::T=0.0) where T<:Float64

Constructs a [`Ellipse`](@ref) `struct` which contains the information required to define an ellipse which may have been rotated and translated.

# Arguments
- `x_radius`: radius of the ellipse in the x axis (i.e. when the rotation, `α`, is zero).
- `y_radius`: radius of the ellipse in the y axis (i.e. when the rotation, `α`, is zero).
- `α`: an angle in radians (0 to 2π) that the ellipse has been rotated by. A positive value represents an anti-clockwise rotation. Default is `0.0`.
- `Cx`: the x coordinate of the centre of the ellipse (the translation of the ellipse in the x axis). Default is `0.0`.
- `Cy`: the y coordinate of the centre of the ellipse (the translation of the ellipse in the y axis). Default is `0.0`.

# Details
The general equation for a rotated and translated ellipse is given by:

```math
1 = Ax^2 + Bxy + Cy^2
```

Where:
```math
A = (\\frac{\\cos^2(α)}{a^2} + \\frac{\\sin^2(α)}{b^2})x^2
```
```math
B = 2\\cos(α)sin(α)(\\frac{1}{a^2} - \\frac{1}{b^2}x^2)xy
```
```math
C = (\\frac{\\sin^2(α)}{a^2} + \\frac{\\cos^2(α)}{b^2})y^2
```
"""
function construct_ellipse(x_radius::T, y_radius::T, α::T=0.0, Cx::T=0.0, Cy::T=0.0) where T<:Float64

    a = max(x_radius, y_radius)
    b = min(x_radius, y_radius)
    assert_parameters_are_valid(a, b, x_radius, y_radius)

    m = eccentricity_squared(a, b)
    C = circumference(a, b)

    return Ellipse(x_radius, y_radius, α, Cx, Cy, a, b, m, C)
end