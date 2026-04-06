# boundary_layer.jl -- BL residuals, initialization, and transition (Phase 5)
# Translated from mfoil.py boundary layer functions

#-------------------------------------------------------------------------------
function build_param(M::Mfoil, si::Int)
    # Builds a parameter structure for side si
    # INPUT
    #   si : side number, 1 = lower, 2 = upper, 3 = wake (Julia 1-based)
    # OUTPUT
    #   param : copy of M.param with side information
    param = deepcopy(M.param)
    param.wake = (si == 3)
    param.turb = param.wake  # the wake is fully turbulent
    param.simi = false       # true for similarity station
    return param
end


#-------------------------------------------------------------------------------
function station_param!(M::Mfoil, param::Param, i::Int)
    # Modifies parameter structure to be specific for station i
    # INPUT
    #   i     : station number (node index along the surface, 1-based)
    #   param : parameter structure to modify in place
    param.turb = M.vsol.turb[i]          # turbulent flag
    param.simi = i in M.isol.Istag       # similarity (stagnation) flag
end


#-------------------------------------------------------------------------------
function stagnation_state(U::AbstractMatrix, x::AbstractVector)
    # Extrapolates two states in U to stagnation (x=0)
    # INPUT
    #   U  : 4x2 matrix [U1 U2], states at first two nodes
    #   x  : 2-vector [x1, x2], x-locations of first two nodes
    # OUTPUT
    #   Ust   : 4-vector, stagnation state
    #   Ust_U : 4x8 Jacobian of Ust w.r.t. [U1; U2]
    #   Ust_x : 4x2 Jacobian of Ust w.r.t. [x1, x2]
    #   xst   : stagnation point location (small, near 0)

    @views U1 = U[:, 1]; U2 = U[:, 2]
    x1, x2 = x[1], x[2]
    dx = x2 - x1; dx_x = [-1.0, 1.0]
    rx = x2 / x1; rx_x = [-rx, 1.0] / x1

    # linear extrapolation weights and stagnation state
    w1 = x2 / dx; w1_x = -w1 / dx * dx_x + [0.0, 1.0] / dx
    w2 = -x1 / dx; w2_x = -w2 / dx * dx_x + [-1.0, 0.0] / dx
    Ust = U1 * w1 + U2 * w2

    # quadratic extrapolation of edge velocity for better slope, ue=K*x
    wk1 = rx / dx; wk1_x = rx_x / dx - wk1 / dx * dx_x
    wk2 = -1.0 / (rx * dx); wk2_x = -wk2 * (rx_x / rx + dx_x / dx)
    K = wk1 * U1[4] + wk2 * U2[4]
    K_U = [0.0, 0.0, 0.0, wk1, 0.0, 0.0, 0.0, wk2]
    K_x = U1[4] * wk1_x + U2[4] * wk2_x

    # stagnation coord cannot be zero, but must be small
    xst = 1e-6
    Ust[4] = K * xst  # linear dep of ue on x near stagnation

    # Ust_U: 4x8 Jacobian
    # rows 1-3: linear extrapolation of th, ds, sa
    # row 4: K*xst with K depending on U1[4] and U2[4]
    Ust_U = zeros(4, 8)
    for j in 1:3
        Ust_U[j, j] = w1       # d(Ust_j)/d(U1_j)
        Ust_U[j, j + 4] = w2   # d(Ust_j)/d(U2_j)
    end
    Ust_U[4, :] = K_U * xst

    # Ust_x: 4x2 Jacobian
    Ust_x = zeros(4, 2)
    for j in 1:3
        Ust_x[j, :] = U1[j] * w1_x + U2[j] * w2_x
    end
    Ust_x[4, :] = K_x * xst

    return Ust, Ust_U, Ust_x, xst
end


#-------------------------------------------------------------------------------
function thwaites_init(K::Float64, nu::Float64)
    # Uses Thwaites correlation to initialize first node in stagnation point flow
    # INPUT
    #   K  : stagnation point constant (ue/x)
    #   nu : kinematic viscosity
    # OUTPUT
    #   th : momentum thickness
    #   ds : displacement thickness
    th = sqrt(0.45 * nu / (6.0 * K))
    ds = 2.2 * th
    return th, ds
end


#-------------------------------------------------------------------------------
function residual_station(param, x::AbstractVector, Uin::AbstractMatrix, Aux::AbstractVector)
    # Calculates the viscous residual at one non-transition station
    # INPUT
    #   param : parameter structure
    #   x     : 2-vector [x1, x2], xi values at the two points
    #   Uin   : 4x2 matrix [U1 U2], states at the two points
    #   Aux   : 2-vector [Aux1, Aux2], auxiliary data (wake gap) at the points
    # OUTPUT
    #   R     : 3-vector, residual (momentum, shape-param, amp/lag)
    #   R_U   : 3x8 Jacobian w.r.t. [U1; U2]
    #   R_x   : 3x2 Jacobian w.r.t. [x1, x2]

    # copy so we do not overwrite Uin
    U = copy(Uin)

    # modify ds to take out wake gap (in Aux) for all calculations below
    U[2, 1] -= Aux[1]
    U[2, 2] -= Aux[2]

    # states
    @views U1 = U[:, 1]; U2 = U[:, 2]
    Um = 0.5 * (U1 + U2)
    th = [U[1, 1], U[1, 2]]
    ds = [U[2, 1], U[2, 2]]
    sa = [U[3, 1], U[3, 2]]

    # speed needs compressibility correction
    uk1, uk1_u = get_uk(U1[4], param)
    uk2, uk2_u = get_uk(U2[4], param)

    # log changes
    thlog = log(th[2] / th[1])
    thlog_U = [-1.0/th[1], 0, 0, 0, 1.0/th[2], 0, 0, 0]
    uelog = log(uk2 / uk1)
    uelog_U = [0, 0, 0, -uk1_u/uk1, 0, 0, 0, uk2_u/uk2]
    xlog = log(x[2] / x[1]); xlog_x = [-1.0/x[1], 1.0/x[2]]
    dx = x[2] - x[1]; dx_x = [-1.0, 1.0]

    # upwinding factor
    upw, upw_U = get_upw(U1, U2, param)

    # shape parameter
    H1, H1_U1 = get_H(U1)
    H2, H2_U2 = get_H(U2)
    H = 0.5 * (H1 + H2)
    H_U = 0.5 * vcat(H1_U1, H2_U2)

    # Hstar = KE shape parameter, averaged
    Hs1, Hs1_U1 = get_Hs(U1, param)
    Hs2, Hs2_U2 = get_Hs(U2, param)
    Hs, Hs_U = upwind(0.5, 0, Hs1, Hs1_U1, Hs2, Hs2_U2)

    # log change in Hstar
    Hslog = log(Hs2 / Hs1)
    Hslog_U = vcat(-1.0/Hs1 * Hs1_U1, 1.0/Hs2 * Hs2_U2)

    # similarity station is special: U1 = U2, x1 = x2
    if param.simi
        thlog = 0.0; thlog_U = thlog_U .* 0.0
        Hslog = 0.0; Hslog_U = Hslog_U .* 0.0
        uelog = 1.0; uelog_U = uelog_U .* 0.0
        xlog = 1.0; xlog_x = [0.0, 0.0]
        dx = 0.5 * (x[1] + x[2]); dx_x = [0.5, 0.5]
    end

    # Hw = wake shape parameter
    Hw1, Hw1_U1 = get_Hw(U1, Aux[1])
    Hw2, Hw2_U2 = get_Hw(U2, Aux[2])
    Hw = 0.5 * (Hw1 + Hw2)
    Hw_U = 0.5 * vcat(Hw1_U1, Hw2_U2)

    # set up shear lag or amplification factor equation
    if param.turb
        # log change of root shear stress coeff
        salog = log(sa[2] / sa[1])
        salog_U = [0, 0, -1.0/sa[1], 0, 0, 0, 1.0/sa[2], 0]

        # BL thickness measure, averaged
        de1, de1_U1 = get_de(U1, param)
        de2, de2_U2 = get_de(U2, param)
        de, de_U = upwind(0.5, 0, de1, de1_U1, de2, de2_U2)

        # normalized slip velocity, averaged
        Us1, Us1_U1 = get_Us(U1, param)
        Us2, Us2_U2 = get_Us(U2, param)
        Us, Us_U = upwind(0.5, 0, Us1, Us1_U1, Us2, Us2_U2)

        # Hk, upwinded
        Hk1, Hk1_U1 = get_Hk(U1, param)
        Hk2, Hk2_U2 = get_Hk(U2, param)
        Hk, Hk_U = upwind(upw, upw_U, Hk1, Hk1_U1, Hk2, Hk2_U2)

        # Re_theta, averaged
        Ret1, Ret1_U1 = get_Ret(U1, param)
        Ret2, Ret2_U2 = get_Ret(U2, param)
        Ret, Ret_U = upwind(0.5, 0, Ret1, Ret1_U1, Ret2, Ret2_U2)

        # skin friction, upwinded
        cf1, cf1_U1 = get_cf(U1, param)
        cf2, cf2_U2 = get_cf(U2, param)
        cf, cf_U = upwind(upw, upw_U, cf1, cf1_U1, cf2, cf2_U2)

        # displacement thickness, averaged
        dsa = 0.5 * (ds[1] + ds[2])
        dsa_U = 0.5 * [0, 1, 0, 0, 0, 1, 0, 0]

        # uq = equilibrium 1/ue * due/dx
        uq, uq_U = get_uq(dsa, dsa_U, cf, cf_U, Hk, Hk_U, Ret, Ret_U, param)

        # cteq = root equilibrium wake layer shear coefficient
        cteq1, cteq1_U1 = get_cteq(U1, param)
        cteq2, cteq2_U2 = get_cteq(U2, param)
        cteq, cteq_U = upwind(upw, upw_U, cteq1, cteq1_U1, cteq2, cteq2_U2)

        # root of shear coefficient (a state), upwinded
        saa, saa_U = upwind(upw, upw_U, sa[1], [0, 0, 1, 0], sa[2], [0, 0, 1, 0])

        # lag coefficient
        Klag = param.SlagK
        beta = param.GB
        Clag = Klag / beta * 1.0 / (1.0 + Us)
        Clag_U = -Clag / (1.0 + Us) * Us_U

        # extra dissipation in wake
        ald = 1.0
        if param.wake; ald = param.Dlr; end

        # shear lag equation
        Rlag = Clag * (cteq - ald * saa) * dx - 2 * de * salog + 2 * de * (uq * dx - uelog) * param.Cuq
        Rlag_U = Clag_U * (cteq - ald * saa) * dx + Clag * (cteq_U - ald * saa_U) * dx -
            2 * de_U * salog - 2 * de * salog_U +
            2 * de_U * (uq * dx - uelog) * param.Cuq + 2 * de * (uq_U * dx - uelog_U) * param.Cuq
        Rlag_x = Clag * (cteq - ald * saa) * dx_x + 2 * de * uq * dx_x

    else
        # laminar, amplification factor equation
        if param.simi
            # similarity station
            Rlag = sa[1] + sa[2]  # no amplification
            Rlag_U = Float64[0, 0, 1, 0, 0, 0, 1, 0]
            Rlag_x = Float64[0, 0]
        else
            # amplification factor equation
            damp1, damp1_U1 = get_damp(U1, param)
            damp2, damp2_U2 = get_damp(U2, param)
            damp, damp_U = upwind(0.5, 0, damp1, damp1_U1, damp2, damp2_U2)

            Rlag = sa[2] - sa[1] - damp * dx
            Rlag_U = Float64[0, 0, -1, 0, 0, 0, 1, 0] - damp_U * dx
            Rlag_x = -damp * dx_x
        end
    end

    # squared mach number, symmetrical average
    Ms1, Ms1_U1 = get_Mach2(U1, param)
    Ms2, Ms2_U2 = get_Mach2(U2, param)
    Ms, Ms_U = upwind(0.5, 0, Ms1, Ms1_U1, Ms2, Ms2_U2)

    # skin friction * x/theta, symmetrical average
    cfxt1, cfxt1_U1, cfxt1_x1 = get_cfxt(U1, x[1], param)
    cfxt2, cfxt2_U2, cfxt2_x2 = get_cfxt(U2, x[2], param)
    cfxtm, cfxtm_Um, cfxtm_xm = get_cfxt(Um, 0.5 * (x[1] + x[2]), param)
    cfxt = 0.25 * cfxt1 + 0.5 * cfxtm + 0.25 * cfxt2
    cfxt_U = 0.25 * vcat(cfxt1_U1 + cfxtm_Um, cfxtm_Um + cfxt2_U2)
    cfxt_x = 0.25 * [cfxt1_x1 + cfxtm_xm, cfxtm_xm + cfxt2_x2]

    # momentum equation
    Rmom = thlog + (2 + H + Hw - Ms) * uelog - 0.5 * xlog * cfxt
    Rmom_U = thlog_U + (H_U + Hw_U - Ms_U) * uelog + (2 + H + Hw - Ms) * uelog_U - 0.5 * xlog * cfxt_U
    Rmom_x = -0.5 * xlog_x * cfxt - 0.5 * xlog * cfxt_x

    # dissipation function times x/theta: cDi = (2*cD/H*)*x/theta, upwinded
    cDixt1, cDixt1_U1, cDixt1_x1 = get_cDixt(U1, x[1], param)
    cDixt2, cDixt2_U2, cDixt2_x2 = get_cDixt(U2, x[2], param)
    cDixt, cDixt_U = upwind(upw, upw_U, cDixt1, cDixt1_U1, cDixt2, cDixt2_U2)
    cDixt_x = [(1.0 - upw) * cDixt1_x1, upw * cDixt2_x2]

    # cf*x/theta, upwinded
    cfxtu, cfxtu_U = upwind(upw, upw_U, cfxt1, cfxt1_U1, cfxt2, cfxt2_U2)
    cfxtu_x = [(1.0 - upw) * cfxt1_x1, upw * cfxt2_x2]

    # Hss = density shape parameter, averaged
    Hss1, Hss1_U1 = get_Hss(U1, param)
    Hss2, Hss2_U2 = get_Hss(U2, param)
    Hss, Hss_U = upwind(0.5, 0, Hss1, Hss1_U1, Hss2, Hss2_U2)

    # shape parameter equation
    Rshape = Hslog + (2 * Hss / Hs + 1 - H - Hw) * uelog + xlog * (0.5 * cfxtu - cDixt)
    Rshape_U = Hslog_U + (2 * Hss_U / Hs - 2 * Hss / Hs^2 * Hs_U - H_U - Hw_U) * uelog +
        (2 * Hss / Hs + 1 - H - Hw) * uelog_U + xlog * (0.5 * cfxtu_U - cDixt_U)
    Rshape_x = xlog_x * (0.5 * cfxtu - cDixt) + xlog * (0.5 * cfxtu_x - cDixt_x)

    # put everything together
    R = [Rmom, Rshape, Rlag]
    R_U = vcat(Rmom_U', Rshape_U', Rlag_U')
    R_x = vcat(Rmom_x', Rshape_x', Rlag_x')

    return R, R_U, R_x
end
