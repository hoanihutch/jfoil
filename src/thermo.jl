# thermo.jl -- thermodynamics functions (Phase 3)

"""
    init_thermo!(M)

Initialize thermodynamic variables in param structure from operating conditions.
"""
function init_thermo!(M)
    g = M.param.gam
    gmi = g - 1.0
    rhoinf = M.oper.rho
    Vinf = M.oper.Vinf
    M.param.Vinf = Vinf
    M.param.muinf = rhoinf * Vinf * M.geom.chord / M.oper.Re
    Minf = M.oper.Ma
    M.param.Minf = Minf

    if Minf > 0
        M.param.KTb = sqrt(1.0 - Minf^2)
        M.param.KTl = Minf^2 / (1.0 + M.param.KTb)^2
        M.param.H0 = (1.0 + 0.5 * gmi * Minf^2) * Vinf^2 / (gmi * Minf^2)
        Tr = 1.0 - 0.5 * Vinf^2 / M.param.H0
        finf = Tr^1.5 * (1.0 + M.param.Tsrat) / (Tr + M.param.Tsrat)
        M.param.cps = 2.0 / (g * Minf^2) * (((1.0 + 0.5 * gmi * Minf^2) / (1.0 + 0.5 * gmi))^(g / gmi) - 1.0)
    else
        finf = 1.0
    end

    M.param.mu0 = M.param.muinf / finf
    M.param.rho0 = rhoinf * (1.0 + 0.5 * gmi * Minf^2)^(1.0 / gmi)
end


"""
    get_cp(u, param)

Pressure coefficient from speed with Karman-Tsien compressibility correction.
Returns (cp, cp_u).
"""
function get_cp(u::Real, param)
    Vinf = param.Vinf
    cp = 1.0 - (u / Vinf)^2
    cp_u = -2.0 * u / Vinf^2

    if param.Minf > 0
        l, b = param.KTl, param.KTb
        den = b + 0.5 * l * (1.0 + b) * cp
        den_cp = 0.5 * l * (1.0 + b)
        cp_u = cp_u * (1.0 - cp / den * den_cp) / den
        cp = cp / den
    end

    return cp, cp_u
end

# Vectorized version
function get_cp(u::AbstractVector, param)
    n = length(u)
    cp = Vector{Float64}(undef, n)
    cp_u = Vector{Float64}(undef, n)
    for i in 1:n
        cp[i], cp_u[i] = get_cp(u[i], param)
    end
    return cp, cp_u
end


"""
    get_uk(u, param)

Karman-Tsien corrected speed from incompressible speed.
Returns (uk, uk_u).
"""
function get_uk(u::Real, param)
    if param.Minf > 0
        l = param.KTl
        Vinf = param.Vinf
        den = 1.0 - l * (u / Vinf)^2
        den_u = -2.0 * l * u / Vinf^2
        uk = u * (1.0 - l) / den
        uk_u = (1.0 - l) / den - (uk / den) * den_u
    else
        uk = u
        uk_u = 1.0
    end

    return uk, uk_u
end

# Vectorized version
function get_uk(u::AbstractVector, param)
    n = length(u)
    uk = Vector{Float64}(undef, n)
    uk_u = Vector{Float64}(undef, n)
    for i in 1:n
        uk[i], uk_u[i] = get_uk(u[i], param)
    end
    return uk, uk_u
end
