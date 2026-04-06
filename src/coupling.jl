# coupling.jl -- viscous matrices and coupling (Phase 6)
# Translated from mfoil.py coupling functions

#-------------------------------------------------------------------------------
"""
    identify_surfaces!(M)

Identify lower/upper/wake surfaces from the stagnation point indices.
Sets M.vsol.Is as a vector of 3 integer vectors (lower, upper, wake).
"""
function identify_surfaces!(M)
    # Julia 1-based: Istag[1] = lower stag index, Istag[2] = upper stag index
    # Lower surface: from Istag[1] down to 1
    # Upper surface: from Istag[2] up to N
    # Wake: from N+1 to N+Nw
    M.vsol.Is = [
        collect(M.isol.Istag[1]:-1:1),
        collect(M.isol.Istag[2]:M.foil.N),
        collect(M.foil.N+1:M.foil.N+M.wake.N)
    ]
    return nothing
end


#-------------------------------------------------------------------------------
"""
    set_wake_gap!(M)

Set height (delta*) of dead air in wake using cubic extrapolation of TE gap.
"""
function set_wake_gap!(M)
    t, hTE, dtdx, tcp, tdp = TE_info(M.foil.x)
    flen = 2.5
    dtdx = clamp(dtdx, -3.0 / flen, 3.0 / flen)
    Lw = flen * hTE
    Nw = M.wake.N
    wgap = zeros(Nw)
    xi_N = M.isol.xi[M.foil.N + 1]  # xi at first wake point (Julia 1-based)
    for i in 1:Nw
        xib = (M.isol.xi[M.foil.N + i] - xi_N) / Lw
        if xib <= 1.0
            wgap[i] = hTE * (1.0 + (2.0 + flen * dtdx) * xib) * (1.0 - xib)^2
        end
    end
    M.vsol.wgap = wgap
    return nothing
end


#-------------------------------------------------------------------------------
"""
    calc_ue_m!(M)

Build the ue-mass sensitivity matrix (D matrix).
Sets M.vsol.ue_sigma, M.vsol.sigma_m, M.vsol.ue_m.
"""
function calc_ue_m!(M)
    @assert length(M.isol.gam) > 0 "No inviscid solution"
    N = M.foil.N
    Nw = M.wake.N
    @assert Nw > 0 "No wake"

    # Cgam = d(wake uei)/d(gamma) [Nw x N]
    Cgam = zeros(Nw, N)
    for i in 1:Nw
        V, V_G = inviscid_velocity(M.foil.x, M.isol.gam, 0.0, 0.0,
                                    M.wake.x[:, i], true)
        @views for j in 1:N
            Cgam[i, j] = V_G[1, j] * M.wake.t[1, i] + V_G[2, j] * M.wake.t[2, i]
        end
    end

    # B = d(airfoil surf streamfunction)/d(source) [(N+1) x (N+Nw-2)]
    # NOTE: sparse candidate
    Npan = N + Nw - 2  # total number of panels
    B = zeros(N + 1, Npan)

    for i in 1:N
        xi = @view M.foil.x[:, i]

        # airfoil panels (constant source)
        for j in 1:N-1
            B[i, j] = panel_constsource_stream(M.foil.x[:, j:j+1], xi)
        end

        # wake panels (piecewise linear source)
        for j in 1:Nw-1
            # Panel endpoint coordinates (Julia 1-based: wake points j, j+1)
            Xj = hcat(M.wake.x[:, j], M.wake.x[:, j+1])
            Xm = 0.5 .* (Xj[:, 1] .+ Xj[:, 2])  # panel midpoint
            Xj3 = hcat(Xj[:, 1], Xm, Xj[:, 2])   # left, mid, right

            if j == Nw - 1
                Xj3[:, 3] = 2.0 .* Xj3[:, 3] .- Xj3[:, 2]  # ghost extension
            end

            # Left half panel
            a, b_val = panel_linsource_stream(Xj3[:, 1:2], xi)
            # Wake columns: airfoil panels cols 1..N-1, wake points cols N..N+Nw-2
            wcol_j = N + j - 1  # column for wake point j (1-based)

            if j > 1
                B[i, wcol_j] += 0.5 * a + b_val
                B[i, wcol_j - 1] += 0.5 * a
            else
                B[i, wcol_j] += b_val
            end

            # Right half panel
            a, b_val = panel_linsource_stream(Xj3[:, 2:3], xi)
            B[i, wcol_j] += a + 0.5 * b_val
            if j < Nw - 1
                B[i, wcol_j + 1] += 0.5 * b_val
            else
                B[i, wcol_j] += 0.5 * b_val
            end
        end
    end

    # Bp = -inv(AIC) * B  -> d(airfoil gamma)/d(source) [N x Npan]
    Bp = -(M.isol.AIC \ B)
    Bp = Bp[1:N, :]  # trim last row (streamfunction equation)

    # Csig = d(wake uei)/d(source) [Nw x Npan]
    # NOTE: sparse candidate
    Csig = zeros(Nw, Npan)

    for i in 1:Nw
        xi = M.wake.x[:, i]
        ti = M.wake.t[:, i]

        # Airfoil constant source panels
        # First/last airfoil panel special for i==1 (first wake point)
        jstart = i == 1 ? 2 : 1
        jend = i == 1 ? N - 2 : N - 1
        for j in jstart:jend
            Csig[i, j] = panel_constsource_velocity(M.foil.x[:, j:j+1], xi, ti)
        end

        # Piecewise linear sources across wake panel halves
        for j in 1:Nw
            # Left, self, right indices (clamped)
            jl = max(j - 1, 1)
            jr = min(j + 1, Nw)

            Xj3 = hcat(M.wake.x[:, jl], M.wake.x[:, j], M.wake.x[:, jr])
            Xj3[:, 1] = 0.5 .* (Xj3[:, 1] .+ Xj3[:, 2])  # left midpoint
            Xj3[:, 3] = 0.5 .* (Xj3[:, 2] .+ Xj3[:, 3])  # right midpoint

            if j == Nw
                Xj3[:, 3] = 2.0 .* Xj3[:, 2] .- Xj3[:, 1]  # ghost extension
            end
            d1 = norm2(Xj3[:, 2] .- Xj3[:, 1])  # left half-panel length
            d2 = norm2(Xj3[:, 3] .- Xj3[:, 2])  # right half-panel length

            # Wake columns: wake point j (1-based) -> column N + j - 1
            wcol_j = N + j - 1
            wcol_jm1 = N + j - 2

            if i == j
                if j == 1  # first point: special TE system
                    dl = norm2(M.foil.x[:, 2] .- M.foil.x[:, 1])
                    du = norm2(M.foil.x[:, N] .- M.foil.x[:, N-1])
                    Csig[i, 1] += (0.5 / π) * (log(dl / d2) + 1)      # lower panel
                    Csig[i, N-1] += (0.5 / π) * (log(du / d2) + 1)    # upper panel
                    Csig[i, wcol_j] += -0.5 / π                        # self
                elseif j == Nw  # last point: no self effect
                    Csig[i, wcol_jm1] += 0.0
                else  # all other points
                    aa = (0.25 / π) * log(d1 / d2)
                    Csig[i, wcol_jm1] += aa + 0.5 / π
                    Csig[i, wcol_j] += aa - 0.5 / π
                end
            else
                if j == 1  # first point: only right half panel
                    a, b_val = panel_linsource_velocity(Xj3[:, 2:3], xi, ti)
                    Csig[i, wcol_j] += b_val       # right half panel
                    Csig[i, 1] += a                 # lower airfoil panel
                    Csig[i, N-1] += a               # upper airfoil panel
                elseif j == Nw  # last point: constant source ghost extension
                    a = panel_constsource_velocity(Xj3[:, [1, 3]], xi, ti)
                    Csig[i, Npan] += a              # last wake column (N+Nw-2)
                else  # left and right half panels
                    a1, b1 = panel_linsource_velocity(Xj3[:, 1:2], xi, ti)
                    a2, b2 = panel_linsource_velocity(Xj3[:, 2:3], xi, ti)
                    Csig[i, wcol_jm1] += a1 + 0.5 * b1
                    Csig[i, wcol_j] += 0.5 * a2 + b2
                end
            end
        end
    end

    # Combine: ue_sigma = d(unsigned ue)/d(source) [(N+Nw) x Npan]
    # Dw = Cgam * Bp + Csig
    Dw = Cgam * Bp .+ Csig
    Dw[1, :] = Bp[N, :]  # first wake point has same ue as TE
    M.vsol.ue_sigma = vcat(Bp, Dw)  # [(N+Nw) x Npan]

    # Build ue_m from ue_sigma using sgnue
    rebuild_ue_m!(M)
    return nothing
end


#-------------------------------------------------------------------------------
"""
    rebuild_ue_m!(M)

Rebuild ue_m matrix after stagnation panel change (new sgnue).
Uses existing ue_sigma.
"""
function rebuild_ue_m!(M)
    @assert size(M.vsol.ue_sigma, 1) > 0 "Need ue_sigma to build ue_m"

    N = M.foil.N
    Nw = M.wake.N
    Npan = N + Nw - 2
    Ntot = N + Nw

    # sigma_m = d(source)/d(mass) [Npan x Ntot]
    # NOTE: sparse candidate
    sigma_m = zeros(Npan, Ntot)

    # Airfoil panels
    for i in 1:N-1
        ds = M.foil.s[i+1] - M.foil.s[i]
        sigma_m[i, i] = M.isol.sgnue[i] * (-1.0) / ds
        sigma_m[i, i+1] = M.isol.sgnue[i+1] * (1.0) / ds
    end

    # Wake panels
    for i in 1:Nw-1
        ds = M.wake.s[i+1] - M.wake.s[i]
        sigma_m[N-1+i, N+i] = -1.0 / ds
        sigma_m[N-1+i, N+i+1] = 1.0 / ds
    end

    M.vsol.sigma_m = sigma_m

    # Sign of ue at all points
    sgue = vcat(M.isol.sgnue, ones(Nw))

    # ue_m = diag(sgue) * ue_sigma * sigma_m [(N+Nw) x (N+Nw)]
    # NOTE: sparse candidate
    M.vsol.ue_m = Diagonal(sgue) * M.vsol.ue_sigma * sigma_m
    return nothing
end


#-------------------------------------------------------------------------------
"""
    wake_sys(M, param)

Construct residual system for wake initialization.
Returns (R, R_U, J) where R is 3x1, R_U is 3x12, J is [il, iu, iw].
"""
function wake_sys(M, param::Param)
    il = M.vsol.Is[1][end]   # lower surface TE index (last of lower)
    Ul = @view M.glob.U[:, il]
    iu = M.vsol.Is[2][end]   # upper surface TE index (last of upper)
    Uu = @view M.glob.U[:, iu]
    iw = M.vsol.Is[3][1]     # first wake index
    Uw = @view M.glob.U[:, iw]

    t, hTE, dtdx, tcp, tdp = TE_info(M.foil.x)

    # Obtain wake shear stress from upper/lower
    p = deepcopy(param)
    p.turb = true
    p.wake = false

    if M.vsol.turb[il]
        ctl = Ul[3]
        ctl_Ul = [0.0, 0.0, 1.0, 0.0]
    else
        ctl, ctl_Ul = get_cttr(collect(Ul), p)
    end

    if M.vsol.turb[iu]
        ctu = Uu[3]
        ctu_Uu = [0.0, 0.0, 1.0, 0.0]
    else
        ctu, ctu_Uu = get_cttr(collect(Uu), p)
    end

    thsum = Ul[1] + Uu[1]
    ctw = (ctl * Ul[1] + ctu * Uu[1]) / thsum
    ctw_Ul = (ctl_Ul .* Ul[1] .+ (ctl - ctw) .* [1.0, 0.0, 0.0, 0.0]) ./ thsum
    ctw_Uu = (ctu_Uu .* Uu[1] .+ (ctu - ctw) .* [1.0, 0.0, 0.0, 0.0]) ./ thsum

    # Residual: delta star in wake includes the TE gap hTE
    R = [Uw[1] - (Ul[1] + Uu[1]),
         Uw[2] - (Ul[2] + Uu[2] + hTE),
         Uw[3] - ctw]

    J = [il, iu, iw]

    # R_U = [R_Ul | R_Uu | R_Uw]  3x12
    R_Ul = zeros(3, 4)
    R_Ul[1, 1] = -1.0
    R_Ul[2, 2] = -1.0
    R_Ul[3, :] = -ctw_Ul

    R_Uu = zeros(3, 4)
    R_Uu[1, 1] = -1.0
    R_Uu[2, 2] = -1.0
    R_Uu[3, :] = -ctw_Uu

    R_Uw = zeros(3, 4)
    R_Uw[1, 1] = 1.0
    R_Uw[2, 2] = 1.0
    R_Uw[3, 3] = 1.0

    R_U = hcat(R_Ul, R_Uu, R_Uw)

    return R, R_U, J
end


#-------------------------------------------------------------------------------
"""
    wake_init(M, ue)

Initialize first wake point state using TE data.
Returns Uw (4-vector).
"""
function wake_init(M, ue::Real)
    iw = M.vsol.Is[3][1]
    Uw = copy(M.glob.U[:, iw])
    R, R_U, J = wake_sys(M, M.param)
    Uw[1:3] .-= R
    Uw[4] = ue
    return Uw
end
