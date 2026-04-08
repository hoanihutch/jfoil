# utils.jl -- utility functions (Phase 0)

"""
    vprint(param, verb, args...)

Conditional printing based on verbosity level.
"""
function vprint(param, verb, args...)
    if verb <= param.verb
        println(args...)
    end
end

"""
    norm2(x)

2-norm of a 2-vector.
"""
@inline function norm2(x::AbstractVector)
    return sqrt(x[1]^2 + x[2]^2)
end

"""
    dist(a, b)

Distance = sqrt(a^2 + b^2).
"""
@inline function dist(a::Real, b::Real)
    return sqrt(a^2 + b^2)
end
