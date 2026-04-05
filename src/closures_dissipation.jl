# closures_dissipation.jl -- dissipation & shear stress closures (Phase 4)
# get_cDi, get_cDi_lam, get_cDi_lamwake, get_cDi_outer, get_cDi_lamstress,
# get_cDi_turbwall, get_cDixt, get_cdutstag, get_cteq, get_cttr

"""
    get_cDi(U, param)

Total dissipation function 2*cD/H*. Returns (cDi, cDi_U).
"""
function get_cDi(U::AbstractVector, param)
    if param.turb
        cDi = 0.0; cDi_U = zeros(4)

        if !param.wake
            cDi0, cDi0_U = get_cDi_turbwall(U, param)
            cDi += cDi0; cDi_U .+= cDi0_U
            cDil, cDil_U = get_cDi_lam(U, param)
        else
            cDil, cDil_U = get_cDi_lamwake(U, param)
        end

        cDi0, cDi0_U = get_cDi_outer(U, param)
        cDi += cDi0; cDi_U .+= cDi0_U

        cDi0, cDi0_U = get_cDi_lamstress(U, param)
        cDi += cDi0; cDi_U .+= cDi0_U

        if cDil > cDi
            cDi = cDil; cDi_U = cDil_U
        end

        if param.wake
            cDi *= 2.0; cDi_U = cDi_U .* 2.0
        end
    else
        cDi, cDi_U = get_cDi_lam(U, param)
    end

    return cDi, cDi_U
end


"""
    get_cDi_turbwall(U, param)

Turbulent wall contribution to dissipation.
"""
function get_cDi_turbwall(U::AbstractVector, param)
    if param.wake
        return 0.0, zeros(4)
    end

    cf, cf_U = get_cf(U, param)
    Hk, Hk_U = get_Hk(U, param)
    Hs, Hs_U = get_Hs(U, param)
    Us, Us_U = get_Us(U, param)
    Ret, Ret_U = get_Ret(U, param)

    lr = log(Ret); lr_U = Ret_U ./ Ret
    Hmin = 1.0 + 2.1 / lr; Hmin_U = -2.1 / lr^2 .* lr_U
    aa = tanh((Hk - 1.0) / (Hmin - 1.0))
    fac = 0.5 + 0.5 * aa
    fac_U = 0.5 * (1.0 - aa^2) .* (Hk_U ./ (Hmin - 1.0) .- (Hk - 1.0) / (Hmin - 1.0)^2 .* Hmin_U)

    cDi = 0.5 * cf * Us * (2.0 / Hs) * fac
    cDi_U = cf_U .* (Us / Hs * fac) .+ cf .* Us_U ./ Hs .* fac .- cDi / Hs .* Hs_U .+ cf * Us / Hs .* fac_U

    return cDi, cDi_U
end


"""
    get_cDi_lam(U, param)

Laminar dissipation function.
"""
function get_cDi_lam(U::AbstractVector, param)
    Hk, Hk_U = get_Hk(U, param)
    Ret, Ret_U = get_Ret(U, param)

    if Hk < 4.0
        num = 0.00205 * (4.0 - Hk)^5.5 + 0.207
        num_Hk = 0.00205 * 5.5 * (4.0 - Hk)^4.5 * (-1.0)
    else
        Hk1 = Hk - 4.0
        num = -0.0016 * Hk1^2 / (1.0 + 0.02 * Hk1^2) + 0.207
        num_Hk = -0.0016 * (2.0 * Hk1 / (1.0 + 0.02 * Hk1^2) - Hk1^2 / (1.0 + 0.02 * Hk1^2)^2 * 0.02 * 2.0 * Hk1)
    end

    cDi = num / Ret
    cDi_U = num_Hk / Ret .* Hk_U .- num / Ret^2 .* Ret_U
    return cDi, cDi_U
end


"""
    get_cDi_lamwake(U, paramin)

Laminar wake dissipation function.
"""
function get_cDi_lamwake(U::AbstractVector, paramin)
    param = deepcopy(paramin)
    param.turb = false

    Hk, Hk_U = get_Hk(U, param)
    Hs, Hs_U = get_Hs(U, param)
    Ret, Ret_U = get_Ret(U, param)
    HsRet = Hs * Ret
    HsRet_U = Hs_U .* Ret .+ Hs .* Ret_U

    num = 2.0 * 1.1 * (1.0 - 1.0 / Hk)^2 * (1.0 / Hk)
    num_Hk = 2.0 * 1.1 * (2.0 * (1.0 - 1.0 / Hk) * (1.0 / Hk^2) * (1.0 / Hk) + (1.0 - 1.0 / Hk)^2 * (-1.0 / Hk^2))
    cDi = num / HsRet
    cDi_U = num_Hk .* Hk_U ./ HsRet .- num / HsRet^2 .* HsRet_U
    return cDi, cDi_U
end


"""
    get_cDi_outer(U, param)

Turbulent outer layer contribution to dissipation.
"""
function get_cDi_outer(U::AbstractVector, param)
    if !param.turb
        return 0.0, zeros(4)
    end

    Hs, Hs_U = get_Hs(U, param)
    Us, Us_U = get_Us(U, param)

    ct = U[3]^2; ct_U = [0.0, 0.0, 2.0 * U[3], 0.0]

    cDi = ct * (0.995 - Us) * 2.0 / Hs
    cDi_U = ct_U .* ((0.995 - Us) * 2.0 / Hs) .+ ct .* (-Us_U) .* 2.0 / Hs .- ct * (0.995 - Us) * 2.0 / Hs^2 .* Hs_U

    return cDi, cDi_U
end


"""
    get_cDi_lamstress(U, param)

Laminar stress contribution to dissipation.
"""
function get_cDi_lamstress(U::AbstractVector, param)
    Hs, Hs_U = get_Hs(U, param)
    Us, Us_U = get_Us(U, param)
    Ret, Ret_U = get_Ret(U, param)
    HsRet = Hs * Ret
    HsRet_U = Hs_U .* Ret .+ Hs .* Ret_U

    num = 0.15 * (0.995 - Us)^2 * 2.0
    num_Us = 0.15 * 2.0 * (0.995 - Us) * (-1.0) * 2.0
    cDi = num / HsRet
    cDi_U = num_Us .* Us_U ./ HsRet .- num / HsRet^2 .* HsRet_U
    return cDi, cDi_U
end


"""
    get_cDixt(U, x, param)

cDi*x/theta from state. Returns (cDixt, cDixt_U, cDixt_x).
"""
function get_cDixt(U::AbstractVector, x::Real, param)
    cDi, cDi_U = get_cDi(U, param)
    cDixt = cDi * x / U[1]
    cDixt_U = cDi_U .* (x / U[1])
    cDixt_U[1] -= cDixt / U[1]
    cDixt_x = cDi / U[1]
    return cDixt, cDixt_U, cDixt_x
end


"""
    get_cdutstag(Uin, param)

cDi*ue*theta at stagnation (laminar only). Returns (D, D_U).
"""
function get_cdutstag(Uin::AbstractVector, param)
    U = copy(Uin)
    U[4] = 0.0
    Hk, Hk_U = get_Hk(U, param)

    if Hk < 4.0
        num = 0.00205 * (4.0 - Hk)^5.5 + 0.207
        num_Hk = 0.00205 * 5.5 * (4.0 - Hk)^4.5 * (-1.0)
    else
        Hk1 = Hk - 4.0
        num = -0.0016 * Hk1^2 / (1.0 + 0.02 * Hk1^2) + 0.207
        num_Hk = -0.0016 * (2.0 * Hk1 / (1.0 + 0.02 * Hk1^2) - Hk1^2 / (1.0 + 0.02 * Hk1^2)^2 * 0.02 * 2.0 * Hk1)
    end

    nu = param.mu0 / param.rho0
    D = nu * num
    D_U = nu * num_Hk .* Hk_U
    return D, D_U
end


"""
    get_cteq(U, param)

Equilibrium sqrt(ctau). Returns (cteq, cteq_U).
"""
function get_cteq(U::AbstractVector, param)
    CC = 0.5 / (param.GA^2 * param.GB)
    C = param.GC
    Hk, Hk_U = get_Hk(U, param)
    Hs, Hs_U = get_Hs(U, param)
    H, H_U = get_H(U)
    Ret, Ret_U = get_Ret(U, param)
    Us, Us_U = get_Us(U, param)

    if param.wake
        if Hk < 1.00005
            Hk = 1.00005; Hk_U = Hk_U .* 0.0
        end
        Hkc = Hk - 1.0; Hkc_U = copy(Hk_U)
    else
        if Hk < 1.05
            Hk = 1.05; Hk_U = Hk_U .* 0.0
        end
        Hkc = Hk - 1.0 - C / Ret
        Hkc_U = Hk_U .+ C / Ret^2 .* Ret_U
        if Hkc < 0.01
            Hkc = 0.01; Hkc_U = Hkc_U .* 0.0
        end
    end

    num = CC * Hs * (Hk - 1.0) * Hkc^2
    num_U = CC .* (Hs_U .* (Hk - 1.0) * Hkc^2 .+ Hs .* Hk_U .* Hkc^2 .+ Hs * (Hk - 1.0) * 2.0 * Hkc .* Hkc_U)
    den = (1.0 - Us) * H * Hk^2
    den_U = (-Us_U) .* H * Hk^2 .+ (1.0 - Us) .* H_U .* Hk^2 .+ (1.0 - Us) * H * 2.0 * Hk .* Hk_U
    cteq = sqrt(num / den)
    cteq_U = 0.5 / cteq .* (num_U ./ den .- num / den^2 .* den_U)

    return cteq, cteq_U
end


"""
    get_cttr(U, param)

sqrt(shear stress coefficient) at transition. Returns (cttr, cttr_U).
"""
function get_cttr(U::AbstractVector, param)
    param_copy = deepcopy(param)
    param_copy.wake = false
    cteq, cteq_U = get_cteq(U, param_copy)
    Hk, Hk_U = get_Hk(U, param_copy)
    if Hk < 1.05
        Hk = 1.05; Hk_U = Hk_U .* 0.0
    end
    C, E = param.CtauC, param.CtauE
    c = C * exp(-E / (Hk - 1.0))
    c_U = c * E / (Hk - 1.0)^2 .* Hk_U
    cttr = c * cteq
    cttr_U = c_U .* cteq .+ c .* cteq_U
    return cttr, cttr_U
end
