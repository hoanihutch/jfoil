# closures_friction.jl -- friction & amplification closures (Phase 4)
# get_cf, get_cfxt, get_cfutstag, get_Us, get_upw, upwind, get_uq, get_damp

"""
    get_cf(U, param)

Skin friction coefficient. Returns (cf, cf_U).
"""
function get_cf(U::AbstractVector, param)
    if param.wake
        return 0.0, zeros(4)
    end

    Hk, Hk_U = get_Hk(U, param)
    Ret, Ret_U = get_Ret(U, param)

    if param.turb
        M2, M2_U = get_Mach2(U, param)
        Fc = sqrt(1.0 + 0.5 * (param.gam - 1.0) * M2)
        Fc_U = 0.5 / Fc * 0.5 * (param.gam - 1.0) .* M2_U
        aa = -1.33 * Hk; aa_U = -1.33 .* Hk_U
        if aa < -17
            aa = -20.0 + 3.0 * exp((aa + 17.0) / 3.0)
            aa_U = (aa + 20.0) / 3.0 .* aa_U
        end
        bb = log(Ret / Fc); bb_U = Ret_U ./ Ret .- Fc_U ./ Fc
        if bb < 3
            bb = 3.0; bb_U = bb_U .* 0.0
        end
        bb /= log(10); bb_U = bb_U ./ log(10)
        cc = -1.74 - 0.31 * Hk; cc_U = -0.31 .* Hk_U
        dd = tanh(4.0 - Hk / 0.875); dd_U = (1.0 - dd^2) * (-1.0 / 0.875) .* Hk_U
        cf0 = 0.3 * exp(aa) * bb^cc
        cf0_U = cf0 .* aa_U .+ 0.3 * exp(aa) * cc * bb^(cc - 1.0) .* bb_U .+ cf0 * log(bb) .* cc_U
        cf = (cf0 + 1.1e-4 * (dd - 1.0)) / Fc
        cf_U = (cf0_U .+ 1.1e-4 .* dd_U) ./ Fc .- cf / Fc .* Fc_U
    else
        if Hk < 5.5
            num = 0.0727 * (5.5 - Hk)^3 / (Hk + 1.0) - 0.07
            num_Hk = 0.0727 * (3.0 * (5.5 - Hk)^2 / (Hk + 1.0) * (-1.0) - (5.5 - Hk)^3 / (Hk + 1.0)^2)
        else
            num = 0.015 * (1.0 - 1.0 / (Hk - 4.5))^2 - 0.07
            num_Hk = 0.015 * 2.0 * (1.0 - 1.0 / (Hk - 4.5)) / (Hk - 4.5)^2
        end
        cf = num / Ret
        cf_U = num_Hk / Ret .* Hk_U .- num / Ret^2 .* Ret_U
    end

    return cf, cf_U
end


"""
    get_cfxt(U, x, param)

cf*x/theta from state. Returns (cfxt, cfxt_U, cfxt_x).
"""
function get_cfxt(U::AbstractVector, x::Real, param)
    cf, cf_U = get_cf(U, param)
    cfxt = cf * x / U[1]
    cfxt_U = cf_U .* (x / U[1])
    cfxt_U[1] -= cfxt / U[1]
    cfxt_x = cf / U[1]
    return cfxt, cfxt_U, cfxt_x
end


"""
    get_cfutstag(Uin, param)

cf*ue*theta at stagnation (laminar only). Returns (F, F_U).
"""
function get_cfutstag(Uin::AbstractVector, param)
    U = copy(Uin)
    U[4] = 0.0
    Hk, Hk_U = get_Hk(U, param)

    if Hk < 5.5
        num = 0.0727 * (5.5 - Hk)^3 / (Hk + 1.0) - 0.07
        num_Hk = 0.0727 * (3.0 * (5.5 - Hk)^2 / (Hk + 1.0) * (-1.0) - (5.5 - Hk)^3 / (Hk + 1.0)^2)
    else
        num = 0.015 * (1.0 - 1.0 / (Hk - 4.5))^2 - 0.07
        num_Hk = 0.015 * 2.0 * (1.0 - 1.0 / (Hk - 4.5)) / (Hk - 4.5)^2
    end
    nu = param.mu0 / param.rho0
    F = nu * num
    F_U = nu * num_Hk .* Hk_U
    return F, F_U
end


"""
    get_Us(U, param)

Normalized wall slip velocity.
"""
function get_Us(U::AbstractVector, param)
    Hs, Hs_U = get_Hs(U, param)
    Hk, Hk_U = get_Hk(U, param)
    H, H_U = get_H(U)

    if param.wake && Hk < 1.00005
        Hk = 1.00005; Hk_U = Hk_U .* 0.0
    end
    if !param.wake && Hk < 1.05
        Hk = 1.05; Hk_U = Hk_U .* 0.0
    end

    beta = param.GB; bi = 1.0 / beta
    Us = 0.5 * Hs * (1.0 - bi * (Hk - 1.0) / H)
    Us_U = 0.5 .* Hs_U .* (1.0 - bi * (Hk - 1.0) / H) .+ 0.5 * Hs .* (-bi .* Hk_U ./ H .+ bi * (Hk - 1.0) / H^2 .* H_U)

    if !param.wake && Us > 0.95
        Us = 0.98; Us_U = Us_U .* 0.0
    end
    if !param.wake && Us > 0.99995
        Us = 0.99995; Us_U = Us_U .* 0.0
    end

    return Us, Us_U
end


"""
    get_upw(U1, U2, param)

Local upwind factor (0.5=trap, 1=BE). Returns (upw, upw_U) where upw_U is 1x8.
"""
function get_upw(U1::AbstractVector, U2::AbstractVector, param)
    Hk1, Hk1_U1 = get_Hk(U1, param)
    Hk2, Hk2_U2 = get_Hk(U2, param)
    Z = zeros(length(Hk1_U1))
    Hut = 1.0
    C = param.wake ? 1.0 : 5.0
    Huc = C * Hut / Hk2^2
    Huc_U = vcat(Z, -2.0 * Huc / Hk2 .* Hk2_U2)
    aa = (Hk2 - 1.0) / (Hk1 - 1.0)
    sga = sign(aa)
    la = log(sga * aa)
    la_U = vcat(-1.0 / (Hk1 - 1.0) .* Hk1_U1, 1.0 / (Hk2 - 1.0) .* Hk2_U2)
    Hls = la^2; Hls_U = 2.0 * la .* la_U
    if Hls > 15
        Hls = 15.0; Hls_U = Hls_U .* 0.0
    end
    upw = 1.0 - 0.5 * exp(-Hls * Huc)
    upw_U = -0.5 * exp(-Hls * Huc) .* (-Hls_U .* Huc .- Hls .* Huc_U)
    return upw, upw_U
end


"""
    upwind(upw, upw_U, f1, f1_U1, f2, f2_U2)

Upwind average of two scalars. Returns (f, f_U) where f_U is 1x8.
"""
function upwind(upw::Real, upw_U, f1::Real, f1_U1::AbstractVector, f2::Real, f2_U2::AbstractVector)
    f = (1.0 - upw) * f1 + upw * f2
    f_U = (-upw_U) .* f1 .+ upw_U .* f2 .+ vcat((1.0 - upw) .* f1_U1, upw .* f2_U2)
    return f, f_U
end


"""
    get_uq(ds, ds_U, cf, cf_U, Hk, Hk_U, Ret, Ret_U, param)

Equilibrium 1/ue*due/dx. Returns (uq, uq_U).
"""
function get_uq(ds::Real, ds_U, cf::Real, cf_U, Hk::Real, Hk_U, Ret::Real, Ret_U, param)
    beta, A, C = param.GB, param.GA, param.GC
    if param.wake
        A = A * param.Dlr; C = 0.0
    end
    if param.wake && Hk < 1.00005
        Hk = 1.00005; Hk_U = Hk_U .* 0.0
    end
    if !param.wake && Hk < 1.05
        Hk = 1.05; Hk_U = Hk_U .* 0.0
    end
    Hkc = Hk - 1.0 - C / Ret
    Hkc_U = Hk_U .+ C / Ret^2 .* Ret_U
    if Hkc < 0.01
        Hkc = 0.01; Hkc_U = Hkc_U .* 0.0
    end
    ut = 0.5 * cf - (Hkc / (A * Hk))^2
    ut_U = 0.5 .* cf_U .- 2.0 * (Hkc / (A * Hk)) .* (Hkc_U ./ (A * Hk) .- Hkc / (A * Hk^2) .* Hk_U)
    uq = ut / (beta * ds)
    uq_U = ut_U ./ (beta * ds) .- uq / ds .* ds_U
    return uq, uq_U
end


"""
    get_damp(U, param)

Amplification rate dn/dx for transition prediction. Returns (damp, damp_U).
"""
function get_damp(U::AbstractVector, param)
    Hk, Hk_U = get_Hk(U, param)
    Ret, Ret_U = get_Ret(U, param)
    th = U[1]

    if Hk < 1.05
        Hk = 1.05; Hk_U = Hk_U .* 0.0
    end

    Hmi = 1.0 / (Hk - 1.0)
    Hmi_U = -Hmi^2 .* Hk_U
    aa = 2.492 * Hmi^0.43
    aa_U = 0.43 * aa / Hmi .* Hmi_U
    bb = tanh(14.0 * Hmi - 9.24)
    bb_U = (1.0 - bb^2) * 14.0 .* Hmi_U
    lrc = aa + 0.7 * (bb + 1.0)
    lrc_U = aa_U .+ 0.7 .* bb_U
    lten = log(10)
    lr = log(Ret) / lten
    lr_U = (1.0 / Ret) .* Ret_U ./ lten
    dl = 0.1

    damp = 0.0
    damp_U = zeros(length(U))

    if lr >= lrc - dl
        rn = (lr - (lrc - dl)) / (2.0 * dl)
        rn_U = (lr_U .- lrc_U) ./ (2.0 * dl)
        if rn >= 1.0
            rf = 1.0; rf_U = zeros(length(U))
        else
            rf = 3.0 * rn^2 - 2.0 * rn^3
            rf_U = (6.0 * rn - 6.0 * rn^2) .* rn_U
        end
        ar = 3.87 * Hmi - 2.52
        ar_U = 3.87 .* Hmi_U
        ex = exp(-ar^2)
        ex_U = ex * (-2.0 * ar) .* ar_U
        da = 0.028 * (Hk - 1.0) - 0.0345 * ex
        da_U = 0.028 .* Hk_U .- 0.0345 .* ex_U
        af = -0.05 + 2.7 * Hmi - 5.5 * Hmi^2 + 3.0 * Hmi^3 + 0.1 * exp(-20.0 * Hmi)
        af_U = (2.7 - 11.0 * Hmi + 9.0 * Hmi^2 - 1.0 * exp(-20.0 * Hmi)) .* Hmi_U
        damp = rf * af * da / th
        damp_U = (rf_U .* af .* da .+ rf .* af_U .* da .+ rf .* af .* da_U) ./ th .- damp / th .* [1.0, 0.0, 0.0, 0.0]
    end

    # extra amplification near ncrit
    ncrit = param.ncrit
    Cea = 5.0
    nx = Cea * (U[3] - ncrit)
    nx_U = Cea .* [0.0, 0.0, 1.0, 0.0]
    eex = 1.0 + tanh(nx)
    eex_U = (1.0 - tanh(nx)^2) .* nx_U
    ed = eex * 0.001 / th
    ed_U = eex_U .* 0.001 / th .- ed / th .* [1.0, 0.0, 0.0, 0.0]
    damp += ed
    damp_U = damp_U .+ ed_U

    return damp, damp_U
end
