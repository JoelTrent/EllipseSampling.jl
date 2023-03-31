"""
    eccentricity_squared(a::T, b::T) where T<:Float64

Computes `m`, the eccentricity of the ellipse squared, ``m=e^2``, where ``m=1-\\big(\\frac{b}{a}\\big)^2``.

# Arguments
- `a`: the major radius of the ellipse.
- `b`: the minor radius of the ellipse.

# Details
This relationship between m and e is seen by considering the equation for e: ``e=\\frac{c}{a}``, where ``c^2=\\big(a^2-b^2\\big)`` and ``a>b``.

Replacing c in the equation for e: 
    
```math
e=\\frac{\\sqrt{a^2-b^2}}{a}
```

Substituting the equation for e into m: 
```math
m = \\bigg(\\frac{\\sqrt{a^2-b^2}}{a}\\bigg)^2 = \\frac{a^2-b^2}{a^2} = 1 - \\frac{b^2}{a^2} = 1 - \\bigg(\\frac{b}{a}\\bigg)^2
```
"""
function eccentricity_squared(a::T, b::T) where T<:Float64
    return 1 - (b/a)^2
end

"""
    circumference(a::Float64, b::Float64)

Calculates the circumference of an ellipse using the elliptic integral of the second kind. Uses `Elliptic.jl`.

# Arguments
- `a`: the major radius of the ellipse.
- `b`: the minor radius of the ellipse.
"""
function circumference(a::Float64, b::Float64); return Elliptic.E(1 - (b/a)^2) * 4 * a end
