# closures_shape.jl -- shape parameter closure relations (Phase 4)
# get_H, get_Hw, get_Mach2, get_Hk, get_Hs, get_Hss, get_Ret, get_rho, get_de

"""
    get_H(U)

Shape parameter H = delta*/theta from state vector U = [th, ds, sa, ue].
Returns (H, H_U) where H_U is 1x4 linearization.
"""
@inline function get_H(U::AbstractVector)
    H = U[2] / U[1]
    H_U = SVector(-H / U[1], 1.0 / U[1], 0.0, 0.0)
    return H, H_U
end


"""
    get_Hw(U, wgap)

Wake gap shape parameter Hw = wgap/theta.
"""
@inline function get_Hw(U::AbstractVector, wgap::Real)
    Hw = wgap / U[1]
    Hw_U = SVector(-Hw / U[1], 0.0, 0.0, 0.0)
    return Hw, Hw_U
end


"""
    get_Mach2(U, param)

Squared Mach number from state. Returns (M2, M2_U).
"""
@inline function get_Mach2(U::AbstractVector, param)
    if param.Minf > 0
        H0, g = param.H0, param.gam
        uk, uk_u = get_uk(U[4], param)
        c2 = (g - 1) * (H0 - 0.5 * uk^2)
        c2_uk = (g - 1) * (-uk)
        M2 = uk^2 / c2
        M2_uk = 2 * uk / c2 - M2 / c2 * c2_uk
        M2_U = SVector(0.0, 0.0, 0.0, M2_uk * uk_u)
    else
        M2 = 0.0
        M2_U = zero(SVector{4,Float64})
    end
    return M2, M2_U
end


"""
    get_Hk(U, param)

Kinematic shape parameter Hk (H corrected for compressibility).
"""
@inline function get_Hk(U::AbstractVector, param)
    H, H_U = get_H(U)

    if param.Minf > 0
        M2, M2_U = get_Mach2(U, param)
        den = 1.0 + 0.113 * M2
        den_M2 = 0.113
        Hk = (H - 0.29 * M2) / den
        Hk_U = (H_U .- 0.29 .* M2_U) ./ den .- Hk / den * den_M2 .* M2_U
    else
        Hk, Hk_U = H, H_U
    end

    return Hk, Hk_U
end


"""
    get_Hs(U, param)

Kinetic energy shape parameter H* (Hstar).
"""
@inline function get_Hs(U::AbstractVector, param)
    Hk, Hk_U = get_Hk(U, param)

    if param.wake && Hk < 1.00005
        Hk = 1.00005; Hk_U = zero(SVector{4,Float64})
    end
    if !param.wake && Hk < 1.05
        Hk = 1.05; Hk_U = zero(SVector{4,Float64})
    end

    if param.turb
        Hsmin = 1.5; dHsinf = 0.015
        Ret, Ret_U = get_Ret(U, param)

        Ho = 4.0; Ho_U = zero(SVector{4,Float64})
        if Ret > 400
            Ho = 3.0 + 400.0 / Ret
            Ho_U = -400.0 / Ret^2 .* Ret_U
        end
        Reb = Ret; Reb_U = Ret_U
        if Ret < 200
            Reb = 200.0; Reb_U = zero(SVector{4,Float64})
        end

        if Hk < Ho
            Hr = (Ho - Hk) / (Ho - 1.0)
            Hr_U = (Ho_U .- Hk_U) ./ (Ho - 1.0) .- (Ho - Hk) / (Ho - 1.0)^2 .* Ho_U
            aa = (2.0 - Hsmin - 4.0 / Reb) * Hr^2
            aa_U = (4.0 / Reb^2 .* Reb_U) .* Hr^2 .+ (2.0 - Hsmin - 4.0 / Reb) * 2.0 * Hr .* Hr_U
            Hs = Hsmin + 4.0 / Reb + aa * 1.5 / (Hk + 0.5)
            Hs_U = -4.0 / Reb^2 .* Reb_U .+ aa_U .* 1.5 / (Hk + 0.5) .- aa * 1.5 / (Hk + 0.5)^2 .* Hk_U
        else
            lrb = log(Reb); lrb_U = 1.0 / Reb .* Reb_U
            aa = Hk - Ho + 4.0 / lrb
            aa_U = Hk_U .- Ho_U .- 4.0 / lrb^2 .* lrb_U
            bb = 0.007 * lrb / aa^2 + dHsinf / Hk
            bb_U = 0.007 .* (lrb_U ./ aa^2 .- 2.0 * lrb / aa^3 .* aa_U) .- dHsinf / Hk^2 .* Hk_U
            Hs = Hsmin + 4.0 / Reb + (Hk - Ho)^2 * bb
            Hs_U = -4.0 / Reb^2 .* Reb_U .+ 2.0 * (Hk - Ho) .* (Hk_U .- Ho_U) .* bb .+ (Hk - Ho)^2 .* bb_U
        end

        M2, M2_U = get_Mach2(U, param)
        den = 1.0 + 0.014 * M2; den_M2 = 0.014
        Hs = (Hs + 0.028 * M2) / den
        Hs_U = (Hs_U .+ 0.028 .* M2_U) ./ den .- Hs / den * den_M2 .* M2_U
    else
        a = Hk - 4.35
        if Hk < 4.35
            num = 0.0111 * a^2 - 0.0278 * a^3
            Hs = num / (Hk + 1.0) + 1.528 - 0.0002 * (a * Hk)^2
            Hs_Hk = (0.0111 * 2 * a - 0.0278 * 3 * a^2) / (Hk + 1.0) - num / (Hk + 1.0)^2 - 0.0002 * 2 * a * Hk * (Hk + a)
        else
            Hs = 0.015 * a^2 / Hk + 1.528
            Hs_Hk = 0.015 * 2 * a / Hk - 0.015 * a^2 / Hk^2
        end
        Hs_U = Hs_Hk .* Hk_U
    end

    return Hs, Hs_U
end


"""
    get_Hss(U, param)

Density shape parameter H** from state.
"""
@inline function get_Hss(U::AbstractVector, param)
    M2, M2_U = get_Mach2(U, param)
    Hk, Hk_U = get_Hk(U, param)
    num = 0.064 / (Hk - 0.8) + 0.251
    num_U = -0.064 / (Hk - 0.8)^2 .* Hk_U
    Hss = M2 * num
    Hss_U = M2_U .* num .+ M2 .* num_U
    return Hss, Hss_U
end


"""
    get_Ret(U, param)

Momentum-thickness Reynolds number Re_theta.
"""
@inline function get_Ret(U::AbstractVector, param)
    if param.Minf > 0
        M2, M2_U = get_Mach2(U, param)
        uk, uk_u = get_uk(U[4], param)
        H0 = param.H0; gmi = param.gam - 1.0; Ts = param.Tsrat
        Tr = 1.0 - 0.5 * uk^2 / H0
        Tr_uk = -uk / H0
        f = Tr^1.5 * (1.0 + Ts) / (Tr + Ts)
        f_Tr = 1.5 * f / Tr - f / (Tr + Ts)
        mu = param.mu0 * f
        mu_uk = param.mu0 * f_Tr * Tr_uk
        den = 1.0 + 0.5 * gmi * M2; den_M2 = 0.5 * gmi
        rho = param.rho0 / den^(1.0 / gmi)
        rho_U = (-1.0 / gmi) * rho / den * den_M2 .* M2_U
        Ret = rho * uk * U[1] / mu
        Ret_U = rho_U .* (uk * U[1] / mu) .+ (rho * U[1] / mu - Ret / mu * mu_uk) .* SVector(0.0, 0.0, 0.0, uk_u) .+ rho * uk / mu .* SVector(1.0, 0.0, 0.0, 0.0)
    else
        Ret = param.rho0 * U[1] * U[4] / param.mu0
        Ret_U = SVector(U[4], 0.0, 0.0, U[1]) ./ param.mu0
    end
    return Ret, Ret_U
end


"""
    get_rho(U, param)

Density from isentropic relations.
"""
@inline function get_rho(U::AbstractVector, param)
    if param.Minf > 0
        M2, M2_U = get_Mach2(U, param)
        gmi = param.gam - 1.0
        den = 1.0 + 0.5 * gmi * M2; den_M2 = 0.5 * gmi
        rho = param.rho0 / den^(1.0 / gmi)
        rho_U = (-1.0 / gmi) * rho / den * den_M2 .* M2_U
    else
        rho = param.rho0
        rho_U = zero(SVector{4,Float64})
    end
    return rho, rho_U
end


"""
    get_de(U, param)

Simplified BL thickness measure delta.
"""
@inline function get_de(U::AbstractVector, param)
    Hk, Hk_U = get_Hk(U, param)
    aa = 3.15 + 1.72 / (Hk - 1.0)
    aa_U = -1.72 / (Hk - 1.0)^2 .* Hk_U
    de = U[1] * aa + U[2]
    de_U = SVector(aa, 1.0, 0.0, 0.0) .+ U[1] .* aa_U
    dmx = 12.0
    if de > dmx * U[1]
        de = dmx * U[1]
        de_U = SVector(dmx, 0.0, 0.0, 0.0)
    end
    return de, de_U
end
