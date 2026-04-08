# inviscid.jl -- inviscid solver functions (Phase 3)

"""
    build_gamma!(M, alpha)

Build and solve the inviscid linear system for the given angle of attack.
Sets M.isol.gamref (0,90 deg reference), M.isol.gam, M.isol.AIC.
"""
function build_gamma!(M, alpha::Real)
    N = M.foil.N
    A = zeros(N + 1, N + 1)  # NOTE: sparse candidate
    rhs = zeros(N + 1, 2)
    t, hTE, dtdx, tcp, tdp = TE_info(M.foil.x)
    nogap = abs(hTE) < 1e-10 * M.geom.chord

    vprint(M.param, 1, "\n <<< Solving the inviscid problem >>> \n")

    @inbounds Threads.@threads for i in 1:N
        xi = @view M.foil.x[:, i]
        for j in 1:N-1
            aij, bij = panel_linvortex_stream(M.foil.x[:, j:j+1], xi)
            A[i, j] += aij
            A[i, j+1] += bij
        end
        A[i, N+1] = -1.0  # last unknown = streamfunction value on surface

        # right-hand sides
        rhs[i, :] = [-xi[2], xi[1]]

        # TE source influence
        a = panel_constsource_stream(M.foil.x[:, [N, 1]], xi)
        A[i, 1] += -a * (0.5 * tcp)
        A[i, N] += a * (0.5 * tcp)

        # TE vortex panel
        a_v, b_v = panel_linvortex_stream(M.foil.x[:, [N, 1]], xi)
        A[i, 1] += -(a_v + b_v) * (-0.5 * tdp)
        A[i, N] += (a_v + b_v) * (-0.5 * tdp)
    end

    # special Nth equation (extrapolation of gamma differences) if no gap
    if nogap
        A[N, :] .= 0.0
        A[N, 1] = 1.0
        A[N, 2] = -2.0
        A[N, 3] = 1.0
        A[N, N-2] = -1.0
        A[N, N-1] = 2.0
        A[N, N] = -1.0
    end

    # Kutta condition
    A[N+1, 1] = 1.0
    A[N+1, N] = 1.0

    # Solve
    M.isol.AIC = A
    g = A \ rhs

    M.isol.gamref = g[1:N, :]  # last value is surface streamfunction
    M.isol.gam = M.isol.gamref[:, 1] .* cosd(alpha) .+ M.isol.gamref[:, 2] .* sind(alpha)
end


"""
    get_ueinv(M)

Compute inviscid tangential velocity at every node (airfoil + wake).
"""
function get_ueinv(M)
    @assert length(M.isol.gam) > 0 "No inviscid solution"
    alpha = M.oper.alpha
    cs = [cosd(alpha), sind(alpha)]
    uea = M.isol.sgnue .* (M.isol.gamref * cs)

    if M.oper.viscous && M.wake.N > 0
        uew = M.isol.uewiref * cs
        uew[1] = uea[end]  # continuity of upper surface and wake ue
    else
        uew = Float64[]
    end

    return vcat(uea, uew)
end


"""
    get_ueinvref(M)

Compute 0,90 deg inviscid tangential velocity at all points. Returns (N+Nw) x 2.
"""
function get_ueinvref(M)
    @assert length(M.isol.gam) > 0 "No inviscid solution"
    uearef = M.isol.gamref .* M.isol.sgnue
    if M.oper.viscous && M.wake.N > 0
        uewref = copy(M.isol.uewiref)
        uewref[1, :] = uearef[end, :]
    else
        uewref = zeros(0, 2)
    end
    return vcat(uearef, uewref)
end


"""
    stagpoint_find!(M)

Find the LE stagnation point on the airfoil from the inviscid solution.
"""
function stagpoint_find!(M)
    @assert length(M.isol.gam) > 0 "No inviscid solution"
    N = M.foil.N

    j = 0
    for k in 1:N
        if M.isol.gam[k] > 0
            j = k
            break
        end
    end
    @assert j > 1 && j <= N "no stagnation point"

    I = [j - 1, j]
    G = M.isol.gam[I]
    S = M.foil.s[I]
    M.isol.Istag = I
    den = G[2] - G[1]
    w1 = G[2] / den
    w2 = -G[1] / den
    M.isol.sstag = w1 * S[1] + w2 * S[2]
    M.isol.xstag = M.foil.x[:, j-1] .* w1 .+ M.foil.x[:, j] .* w2
    st_g1 = G[2] * (S[1] - S[2]) / (den^2)
    M.isol.sstag_g = [st_g1, -st_g1]

    sgnue = -ones(N)
    for i in j:N
        sgnue[i] = 1.0
    end
    M.isol.sgnue = sgnue
    M.isol.xi = vcat(abs.(M.foil.s .- M.isol.sstag), M.wake.s .- M.isol.sstag)
end


"""
    inviscid_velocity(X, G, Vinf, alpha, x, dolin)

Inviscid velocity at point x due to gamma distribution G on panels X.
If dolin, also returns linearization V_G [2 x N].
"""
function inviscid_velocity(X::AbstractMatrix, G::AbstractVector, Vinf::Real, alpha::Real, x::AbstractVector, dolin::Bool)
    N = size(X, 2)
    v1 = 0.0; v2 = 0.0
    if dolin
        V_G = zeros(2, N)
    end
    t, hTE, dtdx, tcp, tdp = TE_info(X)

    @inbounds for j in 1:N-1
        a, b = panel_linvortex_velocity(X[:, j:j+1], x, nothing, false)
        v1 += a[1] * G[j] + b[1] * G[j+1]
        v2 += a[2] * G[j] + b[2] * G[j+1]
        if dolin
            V_G[1, j] += a[1]; V_G[2, j] += a[2]
            V_G[1, j+1] += b[1]; V_G[2, j+1] += b[2]
        end
    end

    # TE source influence
    a = panel_constsource_velocity(X[:, [N, 1]], x, nothing)
    f1 = a .* (-0.5 * tcp)
    f2 = a .* (0.5 * tcp)
    v1 += f1[1] * G[1] + f2[1] * G[N]
    v2 += f1[2] * G[1] + f2[2] * G[N]
    if dolin
        V_G[1, 1] += f1[1]; V_G[2, 1] += f1[2]
        V_G[1, N] += f2[1]; V_G[2, N] += f2[2]
    end

    # TE vortex influence
    a_v, b_v = panel_linvortex_velocity(X[:, [N, 1]], x, nothing, false)
    f1 = (a_v .+ b_v) .* (0.5 * tdp)
    f2 = (a_v .+ b_v) .* (-0.5 * tdp)
    v1 += f1[1] * G[1] + f2[1] * G[N]
    v2 += f1[2] * G[1] + f2[2] * G[N]
    if dolin
        V_G[1, 1] += f1[1]; V_G[2, 1] += f1[2]
        V_G[1, N] += f2[1]; V_G[2, N] += f2[2]
    end

    # freestream
    V = SVector(v1 + Vinf * cosd(alpha), v2 + Vinf * sind(alpha))

    if dolin
        return V, V_G
    else
        return V
    end
end


"""
    build_wake!(M)

Build wake panels from the inviscid solution using predictor-corrector method.
"""
function build_wake!(M)
    @assert length(M.isol.gam) > 0 "No inviscid solution"
    N = M.foil.N
    Vinf = M.oper.Vinf
    Nw = Int(ceil(N / 10 + 10 * M.geom.wakelen))
    S = M.foil.s
    ds1 = 0.5 * (S[2] - S[1] + S[N] - S[N-1])
    sv = space_geom(ds1, M.geom.wakelen * M.geom.chord, Nw)
    xyw = zeros(2, Nw)
    tw = zeros(2, Nw)
    xy1x, xy1z = M.foil.x[1, 1], M.foil.x[2, 1]
    xyNx, xyNz = M.foil.x[1, N], M.foil.x[2, N]
    xyte_x = 0.5 * (xy1x + xyNx)
    xyte_z = 0.5 * (xy1z + xyNz)
    nx, nz = xyNx - xy1x, xyNz - xy1z
    tx, tz = nz, -nx
    @assert tx > 0 "Wrong wake direction; ensure airfoil points are CCW"
    xyw[1, 1] = xyte_x + 1e-5 * tx * M.geom.chord
    xyw[2, 1] = xyte_z + 1e-5 * tz * M.geom.chord
    sw = S[N] .+ sv

    for i in 1:Nw-1
        v1 = inviscid_velocity(M.foil.x, M.isol.gam, Vinf, M.oper.alpha, xyw[:, i], false)
        v1 = v1 ./ norm2(v1)
        tw[:, i] = v1
        xyw[:, i+1] = xyw[:, i] .+ (sv[i+1] - sv[i]) .* v1
        v2 = inviscid_velocity(M.foil.x, M.isol.gam, Vinf, M.oper.alpha, xyw[:, i+1], false)
        v2 = v2 ./ norm2(v2)
        tw[:, i+1] = v2
        xyw[:, i+1] = xyw[:, i] .+ (sv[i+1] - sv[i]) .* 0.5 .* (v1 .+ v2)
    end

    # inviscid ue in wake
    uewi = zeros(Nw)
    uewiref = zeros(Nw, 2)
    for i in 1:Nw
        v = inviscid_velocity(M.foil.x, M.isol.gam, Vinf, M.oper.alpha, xyw[:, i], false)
        uewi[i] = dot(v, tw[:, i])
        v0 = inviscid_velocity(M.foil.x, M.isol.gamref[:, 1], Vinf, 0.0, xyw[:, i], false)
        uewiref[i, 1] = dot(v0, tw[:, i])
        v90 = inviscid_velocity(M.foil.x, M.isol.gamref[:, 2], Vinf, 90.0, xyw[:, i], false)
        uewiref[i, 2] = dot(v90, tw[:, i])
    end

    M.wake.N = Nw
    M.wake.x = xyw
    M.wake.s = sw
    M.wake.t = tw
    M.isol.uewi = uewi
    M.isol.uewiref = uewiref
end


"""
    calc_force!(M)

Calculate force and moment coefficients from pressure integration.
"""
function calc_force!(M)
    chord = M.geom.chord
    xref = M.geom.xref
    Vinf = M.param.Vinf
    rho = M.oper.rho
    alpha = M.oper.alpha
    qinf = 0.5 * rho * Vinf^2
    N = M.foil.N

    # pressure coefficient at each node
    ue = M.oper.viscous ? M.glob.U[4, :] : get_ueinv(M)
    cp, cp_ue = get_cp(ue, M.param)
    M.post.cp = cp
    cpi, _ = get_cp(get_ueinv(M), M.param)
    M.post.cpi = cpi

    cl = 0.0
    cl_ue = zeros(N)
    cl_alpha = 0.0
    cm = 0.0
    cdpi = 0.0

    for i0 in 1:N
        if i0 == N
            i = 1
            ip = N
        else
            i = i0 + 1
            ip = i0
        end

        x1 = @view M.foil.x[:, ip]
        x2 = @view M.foil.x[:, i]
        dxv = x2 .- x1
        dx1 = x1 .- xref
        dx2 = x2 .- xref
        dx1nds = dxv[1] * dx1[1] + dxv[2] * dx1[2]
        dx2nds = dxv[1] * dx2[1] + dxv[2] * dx2[2]
        dx = -dxv[1] * cosd(alpha) - dxv[2] * sind(alpha)
        dz = dxv[2] * cosd(alpha) - dxv[1] * sind(alpha)
        cp1, cp2 = cp[ip], cp[i]
        cpbar = 0.5 * (cp1 + cp2)
        cl += dx * cpbar
        cl_ue[ip] += dx * 0.5 * cp_ue[ip]
        cl_ue[i] += dx * 0.5 * cp_ue[i]
        cl_alpha += cpbar * (sind(alpha) * dxv[1] - cosd(alpha) * dxv[2]) * π / 180.0
        cm += cp1 * dx1nds / 3.0 + cp1 * dx2nds / 6.0 + cp2 * dx1nds / 6.0 + cp2 * dx2nds / 3.0
        cdpi += dz * cpbar

        cl /= chord
        cm /= chord^2
        cdpi /= chord
        M.post.cl = cl
        M.post.cl_ue = cl_ue
        M.post.cl_alpha = cl_alpha
        M.post.cm = cm
        M.post.cdpi = cdpi
    end

    # viscous contributions
    cd = 0.0
    cdf = 0.0
    if M.oper.viscous
        # Squire-Young relation for total drag (extrapolates theta from end of wake)
        iw = M.vsol.Is[3][end]  # station at the end of the wake
        Uw = M.glob.U[:, iw]
        Hw, Hw_U = get_H(Uw)
        uk, uk_ue = get_uk(Uw[4], M.param)
        cd = 2.0 * Uw[1] * (uk / Vinf)^((5 + Hw) / 2.0)
        M.post.cd_U = 2.0 * Uw[1] * (uk / Vinf)^((5 + Hw) / 2.0) * log(uk / Vinf) * 0.5 * Hw_U
        M.post.cd_U[1] += 2.0 * (uk / Vinf)^((5 + Hw) / 2.0)
        M.post.cd_U[4] += 2.0 * Uw[1] * (5 + Hw) / 2.0 * (uk / Vinf)^((3 + Hw) / 2.0) * (1.0 / Vinf) * uk_ue

        # skin friction drag
        Df = 0.0
        for si in 1:2
            Is = M.vsol.Is[si]
            param = build_param(M, si)
            station_param!(M, param, Is[1])
            cf1 = 0.0
            ue1 = 0.0
            rho1 = rho
            x1 = M.isol.xstag
            for ii in 1:length(Is)
                station_param!(M, param, Is[ii])
                cf2, cf2_U = get_cf(M.glob.U[:, Is[ii]], param)
                ue2, ue2_ue = get_uk(M.glob.U[4, Is[ii]], param)
                rho2, rho2_U = get_rho(M.glob.U[:, Is[ii]], param)
                x2 = M.foil.x[:, Is[ii]]
                dxv = x2 .- x1
                dx_cf = dxv[1] * cosd(alpha) + dxv[2] * sind(alpha)
                Df += 0.25 * (rho1 * cf1 * ue1^2 + rho2 * cf2 * ue2^2) * dx_cf
                cf1 = cf2; ue1 = ue2; x1 = x2; rho1 = rho2
            end
        end
        cdf = Df / (qinf * chord)
    end

    M.post.cd = cd
    M.post.cdf = cdf
    M.post.cdp = cd - cdf

    s = @sprintf("  alpha=%.2fdeg, cl=%.6f, cm=%.6f, cdpi=%.6f, cd=%.6f, cdf=%.6f, cdp=%.6f",
        M.oper.alpha, M.post.cl, M.post.cm, M.post.cdpi, M.post.cd, M.post.cdf, M.post.cdp)
    vprint(M.param, 1, s)
end


"""
    solve_inviscid!(M)

Top-level inviscid solve.
"""
function solve_inviscid!(M)
    @assert M.foil.N > 0 "No panels"
    M.oper.viscous = false
    init_thermo!(M)
    M.isol.sgnue = ones(M.foil.N)
    build_gamma!(M, M.oper.alpha)
    calc_force!(M)
    M.glob.conv = true
end


"""
    rebuild_isol!(M)

Rebuild inviscid solution after angle of attack change.
"""
function rebuild_isol!(M)
    @assert length(M.isol.gam) > 0 "No inviscid solution"
    vprint(M.param, 2, "\n  Rebuilding the inviscid solution.")
    alpha = M.oper.alpha
    M.isol.gam = M.isol.gamref[:, 1] .* cosd(alpha) .+ M.isol.gamref[:, 2] .* sind(alpha)
    if !M.oper.viscous
        stagpoint_find!(M)
    elseif M.oper.redowake
        build_wake!(M)
        identify_surfaces!(M)
        calc_ue_m!(M)
    end
end
