# solver.jl -- coupled Newton solver (Phase 7)
# Translated from mfoil.py solver functions

#-------------------------------------------------------------------------------
"""
    init_boundary_layer!(M)

Initialize BL solution on foil and wake by marching with given edge velocity.
"""
function init_boundary_layer!(M)
    Hmaxl = 3.8  # above this, laminar separation
    Hmaxt = 2.5  # above this, turbulent separation

    ueinv = get_ueinv(M)  # inviscid velocity

    M.glob.Nsys = M.foil.N + M.wake.N

    # do we need to initialize?
    if !M.oper.initbl && size(M.glob.U, 2) == M.glob.Nsys
        vprint(M.param, 1, "\n <<< Starting with current boundary layer >>> \n")
        M.glob.U[4, :] = ueinv  # set new edge velocity
        return
    end

    vprint(M.param, 1, "\n <<< Initializing the boundary layer >>> \n")

    M.glob.U = zeros(4, M.glob.Nsys)
    M.vsol.turb = falses(M.glob.Nsys)

    for si in 1:3  # loop over surfaces
        vprint(M.param, 3, @sprintf("\nSide is = %d:\n", si))

        Is = M.vsol.Is[si]
        xi = M.isol.xi[Is]
        ue = ueinv[Is] |> copy
        N = length(Is)
        U = zeros(4, N)
        Aux = zeros(N)

        # ensure edge velocities are not tiny
        uemax = maximum(abs.(ue))
        for i in 1:N
            ue[i] = max(ue[i], 1e-8 * uemax)
        end

        # get parameter structure
        param = build_param(M, si)

        # set auxiliary data
        if si == 3
            Aux[:] = M.vsol.wgap
        end

        # initialize state at first point
        i0 = 1
        if si < 3
            # Solve for the stagnation state (Thwaites initialization + Newton)
            if xi[1] < 1e-8 * xi[end]
                K = ue[2] / xi[2]; hitstag = true
            else
                K = ue[1] / xi[1]; hitstag = false
            end
            th, ds = thwaites_init(K, param.mu0 / param.rho0)
            xst = 1.0e-6
            Ust = [th, ds, 0.0, K * xst]
            nNewton = 20
            for iNewton in 0:nNewton-1
                param.turb = false; param.simi = true
                R, R_U, R_x = residual_station(param, SVector(xst, xst), hcat(SVector{4}(Ust), SVector{4}(Ust)), SVector(0.0, 0.0))
                param.simi = false
                if norm(R) < 1e-10; break; end
                A = R_U[:, 5:7] + R_U[:, 1:3]
                b = -R
                sol = A \ b
                dU = SVector(sol[1], sol[2], sol[3], 0.0)
                dm = max(abs(dU[1] / Ust[1]), abs(dU[2] / Ust[2]))
                omega = dm < 0.2 ? 1.0 : 0.2 / dm
                dU = dU .* omega
                Ust .+= dU
            end

            if hitstag
                U[:, 1] = Ust; U[4, 1] = ue[1]; i0 = 2
            end
            U[:, i0] = Ust; U[4, i0] = ue[i0]

        else  # wake
            U[:, 1] = wake_init(M, ue[1])
            param.turb = true
            M.vsol.turb[Is[1]] = true
        end

        # march over rest of points
        tran = false
        i = i0 + 1
        while i <= N
            Ip = SVector(i - 1, i)
            U[:, i] = U[:, i-1]; U[4, i] = ue[i]  # guess
            if tran
                ct, ct_U = get_cttr(U[:, i], param); U[3, i] = ct
            end
            M.vsol.turb[Is[i]] = tran || param.turb
            direct = true
            nNewton = 30; iNswitch = 12
            Hktgt = 0.0  # will be set if switching to inverse
            iNewton = 0
            for iN in 0:nNewton-1
                iNewton = iN
                # call residual at this station
                if tran
                    vprint(param, 4, @sprintf("i=%d, residual_transition (iNewton = %d)", i, iNewton))
                    try
                        R, R_U, R_x = residual_transition(M, param, xi[Ip], U[:, Ip], Aux[Ip])
                    catch
                        vprint(param, 1, "Transition calculation failed in BL init. Continuing.")
                        M.vsol.xt = 0.5 * sum(xi[Ip])
                        U[:, i] = U[:, i-1]; U[4, i] = ue[i]; U[3, i] = ct
                        R = zero(SVector{3,Float64})  # so we move on
                    end
                else
                    vprint(param, 4, @sprintf("i=%d, residual_station (iNewton = %d)", i, iNewton))
                    R, R_U, R_x = residual_station(param, xi[Ip], U[:, Ip], Aux[Ip])
                end
                if norm(R) < 1e-10; break; end

                if direct  # direct mode => ue prescribed => solve for th, ds, sa
                    A = R_U[:, 5:7]; b = -R; dU = vcat(A \ b, 0.0)
                else  # inverse mode => Hk prescribed
                    Hk, Hk_U = get_Hk(U[:, i], param)
                    A = vcat(R_U[:, 5:8], Hk_U')
                    b = vcat(-R, Hktgt - Hk)
                    dU = A \ b
                end

                # under-relaxation
                dm = max(abs(dU[1] / U[1, i-1]), abs(dU[2] / U[2, i-1]))
                if !direct; dm = max(dm, abs(dU[4] / U[4, i-1])); end
                if param.turb; dm = max(dm, abs(dU[3] / U[3, i-1]))
                elseif direct; dm = max(dm, abs(dU[3] / 10)); end

                omega = dm > 0.3 ? 0.3 / dm : 1.0
                dU = dU .* omega

                Ui = Vector(U[:, i] + dU)

                # clip extreme values
                if param.turb; Ui[3] = clamp(Ui[3], 1e-7, 0.3); end

                # check if about to separate
                Hmax = param.turb ? Hmaxt : Hmaxl
                Hk, Hk_U = get_Hk(Ui, param)

                if direct && (Hk > Hmax || iNewton > iNswitch)
                    direct = false
                    vprint(param, 2, @sprintf("** switching to inverse: i=%d, iNewton=%d", i, iNewton))
                    Hk, Hk_U = get_Hk(U[:, i-1], param)
                    Hkr = (xi[i] - xi[i-1]) / U[1, i-1]
                    if param.wake
                        H2 = Hk
                        for k in 1:6
                            H2 -= (H2 + 0.03 * Hkr * (H2 - 1)^3 - Hk) / (1 + 0.09 * Hkr * (H2 - 1)^2)
                        end
                        Hktgt = max(H2, 1.01)
                    elseif param.turb
                        Hktgt = Hk - 0.15 * Hkr
                    else
                        Hktgt = Hk + 0.03 * Hkr
                    end
                    if !param.wake; Hktgt = max(Hktgt, Hmax); end
                    if iNewton > iNswitch
                        U[:, i] = U[:, i-1]; U[4, i] = ue[i]
                    end
                else
                    U[:, i] = Ui
                end
            end

            if iNewton >= nNewton - 1
                vprint(param, 1, @sprintf("** BL init not converged: si=%d, i=%d **\n", si, i))
                U[:, i] = U[:, i-1]; U[4, i] = ue[i]
                if si < 3
                    U[1, i] = U[1, i-1] * sqrt(xi[i] / xi[i-1])
                    U[2, i] = U[2, i-1] * sqrt(xi[i] / xi[i-1])
                else
                    rlen = (xi[i] - xi[i-1]) / (10.0 * U[2, i-1])
                    U[2, i] = (U[2, i-1] + U[1, i-1] * rlen) / (1.0 + rlen)
                end
            end

            # check for transition
            if !param.turb && !tran && U[3, i] > param.ncrit
                vprint(param, 2, @sprintf("Identified transition at (si=%d, i=%d): n=%.5f, ncrit=%.5f", si, i, U[3, i], param.ncrit))
                tran = true
                continue  # redo station with transition
            end

            if tran
                store_transition!(M, si, i)
                param.turb = true; tran = false
            end

            i += 1
        end

        # store states
        M.glob.U[:, Is] = U
    end
    return nothing
end


#-------------------------------------------------------------------------------
"""
    stagpoint_move!(M)

Move the LE stagnation point on the airfoil using the global solution ue.
"""
function stagpoint_move!(M)
    N = M.foil.N
    I = copy(M.isol.Istag)
    ue = M.glob.U[4, :]
    sstag0 = M.isol.sstag

    newpanel = true
    if ue[I[2]] < 0
        # move stagnation point up (larger s, new panel)
        vprint(M.param, 2, "  Moving stagnation point up")
        j = I[2]
        while j <= N
            if ue[j] > 0; break; end
            j += 1
        end
        @assert j <= N "no stagnation point"
        I1 = j
        for jj in I[2]:I1-1
            ue[jj] *= -1.0
        end
        I[1] = I1 - 1; I[2] = I1
    elseif ue[I[1]] < 0
        # move stagnation point down (smaller s, new panel)
        vprint(M.param, 2, "  Moving stagnation point down")
        j = I[1]
        while j >= 1
            if ue[j] > 0; break; end
            j -= 1
        end
        @assert j >= 1 "no stagnation point"
        I0 = j
        for jj in I0+1:I[1]
            ue[jj] *= -1.0
        end
        I[1] = I0; I[2] = I0 + 1
    else
        newpanel = false
    end

    # move point along panel
    ues = ue[I]; S = M.foil.s[I]
    @assert ues[1] > 0 && ues[2] > 0 "stagpoint_move: velocity error"
    den = ues[1] + ues[2]
    w1 = ues[2] / den; w2 = ues[1] / den
    M.isol.sstag = w1 * S[1] + w2 * S[2]
    M.isol.xstag = M.foil.x[:, I] * [w1, w2]
    M.isol.sstag_ue = [ues[2], -ues[1]] .* (S[2] - S[1]) / (den * den)
    vprint(M.param, 2, @sprintf("  Moving stagnation point: s=%.15e -> s=%.15e", sstag0, M.isol.sstag))

    # set new xi coordinates for every point
    M.isol.xi = vcat(abs.(M.foil.s .- M.isol.sstag), M.wake.s .- M.isol.sstag)

    # matrices need to be recalculated if on a new panel
    if newpanel
        vprint(M.param, 2, @sprintf("  New stagnation panel = %d %d", I[1], I[2]))
        M.isol.Istag = I
        for i in 1:I[1]
            M.isol.sgnue[i] = -1.0
        end
        for i in I[1]+1:N
            M.isol.sgnue[i] = 1.0
        end
        identify_surfaces!(M)
        M.glob.U[4, :] = ue  # sign of ue changed on some points
        M.glob.realloc = true
        rebuild_ue_m!(M)
    end
    return nothing
end


#-------------------------------------------------------------------------------
"""
    build_glob_sys!(M)

Assemble global BL residual and Jacobian for the coupled problem.
"""
function build_glob_sys!(M)
    Nsys = M.glob.Nsys

    # allocate or zero arrays
    do_alloc = M.glob.realloc || length(M.glob.R) != 3 * Nsys
    if do_alloc
        M.glob.R = zeros(3 * Nsys)
    else
        fill!(M.glob.R, 0.0)
    end

    alloc_R_U = M.glob.realloc || size(M.glob.R_U) != (3 * Nsys, 4 * Nsys)
    if alloc_R_U
        M.glob.R_U = zeros(3 * Nsys, 4 * Nsys)  # NOTE: sparse candidate
    else
        fill!(M.glob.R_U, 0.0)
    end

    alloc_R_x = M.glob.realloc || size(M.glob.R_x) != (3 * Nsys, Nsys)
    if alloc_R_x
        M.glob.R_x = zeros(3 * Nsys, Nsys)  # NOTE: sparse candidate
    else
        fill!(M.glob.R_x, 0.0)
    end

    M.glob.realloc = false

    for si in 1:3  # loop over surfaces
        Is = M.vsol.Is[si]
        xi = M.isol.xi[Is]
        N = length(Is)
        U = M.glob.U[:, Is]

        param = build_param(M, si)

        # Aux: wake gap for wake surface, zeros for airfoil surfaces
        Aux = si == 3 ? M.vsol.wgap : zeros(N)

        # special case of tiny first xi
        i0 = (si < 3 && xi[1] < 1e-8 * xi[end]) ? 2 : 1

        # first point system
        if si < 3
            Ip = SVector(i0, i0 + 1)
            Ust, Ust_U, Ust_x, xst = stagnation_state(U[:, Ip], xi[Ip])
            param.turb = false; param.simi = true
            R1, R1_Ut, R1_x = residual_station(param, SVector(xst, xst), hcat(Ust, Ust), SVector(Aux[i0], Aux[i0]))
            param.simi = false
            R1_Ust = R1_Ut[:, 1:4] + R1_Ut[:, 5:8]
            R1_U = R1_Ust * Ust_U
            R1_x = R1_Ust * Ust_x
            J = SVector(Is[i0], Is[i0 + 1])

            if i0 == 2
                # i0=1 landed on stagnation: set value to Ust
                vprint(param, 2, "hit stagnation!")
                Ig = (3*(Is[1]-1)+1):(3*(Is[1]-1)+3)
                M.glob.R[Ig] = U[1:3, 1] - Ust[1:3]
                M.glob.R_U[Ig, (4*(Is[1]-1)+1):(4*(Is[1]-1)+4)] += SMatrix{3,4}(1,0,0, 0,1,0, 0,0,1, 0,0,0)
                M.glob.R_U[Ig, (4*(J[1]-1)+1):(4*(J[1]-1)+4)] -= Ust_U[1:3, 1:4]
                M.glob.R_U[Ig, (4*(J[2]-1)+1):(4*(J[2]-1)+4)] -= Ust_U[1:3, 5:8]
                M.glob.R_x[Ig, J] = -Ust_x[1:3, :]
            end
        else
            # wake initialization
            R1, R1_U, J = wake_sys(M, param)
            R1_x = zeros(0, 0)
            param.turb = true; param.wake = true
        end

        # store first point system in global residual, Jacobian
        Ig = (3*(Is[i0]-1)+1):(3*(Is[i0]-1)+3)
        M.glob.R[Ig] = R1
        for j in 1:length(J)
            jcols = (4*(J[j]-1)+1):(4*(J[j]-1)+4)
            M.glob.R_U[Ig, jcols] += R1_U[:, (4*(j-1)+1):(4*(j-1)+4)]
            if size(R1_x, 1) > 0
                M.glob.R_x[Ig, J[j]] += R1_x[:, j]
            end
        end

        # march over rest of points
        @inbounds for i in (i0 + 1):N
            Ip = SVector(i - 1, i)
            tran = M.vsol.turb[Is[i-1]] ⊻ M.vsol.turb[Is[i]]

            if tran
                Ri, Ri_U, Ri_x = residual_transition(M, param, xi[Ip], U[:, Ip], Aux[Ip])
                store_transition!(M, si, i)
            else
                Ri, Ri_U, Ri_x = residual_station(param, xi[Ip], U[:, Ip], Aux[Ip])
            end

            Ig = (3*(Is[i]-1)+1):(3*(Is[i]-1)+3)
            M.glob.R[Ig] += Ri
            M.glob.R_U[Ig, (4*(Is[i-1]-1)+1):(4*(Is[i-1]-1)+4)] += Ri_U[:, 1:4]
            M.glob.R_U[Ig, (4*(Is[i]-1)+1):(4*(Is[i]-1)+4)] += Ri_U[:, 5:8]
            M.glob.R_x[Ig, SVector(Is[i-1], Is[i])] += Ri_x

            if tran; param.turb = true; end
        end
    end

    # include effects of R_x into R_U: R_ue += R_x*x_st*st_ue
    Nsys = M.glob.Nsys
    Iue = 4:4:4*Nsys  # ue indices in U (Julia: every 4th starting from 4)
    x_st = -M.isol.sgnue
    x_st = vcat(x_st, -ones(M.wake.N))  # wake same sens as upper surface
    R_st = M.glob.R_x * x_st  # [3*Nsys]
    Ist = M.isol.Istag
    st_ue = M.isol.sstag_ue
    M.glob.R_U[:, Iue[Ist[1]]] += R_st .* st_ue[1]
    M.glob.R_U[:, Iue[Ist[2]]] += R_st .* st_ue[2]

    return nothing
end


#-------------------------------------------------------------------------------
"""
    clalpha_residual(M)

Compute cl constraint (or alpha-prescribed) residual and Jacobian.
Returns (Rcla, Ru_alpha, Rcla_U).
"""
function clalpha_residual(M)
    Nsys = M.glob.Nsys
    N = M.foil.N
    alpha = M.oper.alpha

    if M.oper.givencl
        Rcla = M.post.cl - M.oper.cltgt
        Rcla_U = zeros(4 * Nsys + 1)
        Rcla_U[end] = M.post.cl_alpha
        # cl_ue only affects airfoil nodes
        for i in 1:N
            Rcla_U[4*(i-1)+4] = M.post.cl_ue[i]
        end
        # Ru_alpha: ue residual sensitivity to alpha
        uref = get_ueinvref(M)
        Ru_alpha = -(uref * [-sind(alpha), cosd(alpha)]) .* (π / 180)
    else
        Rcla = 0.0
        Ru_alpha = zeros(Nsys)
        Rcla_U = zeros(4 * Nsys + 1)
        Rcla_U[end] = 1.0
    end

    return Rcla, Ru_alpha, Rcla_U
end


#-------------------------------------------------------------------------------
"""
    solve_glob!(M)

Solve global system for the primary variable update dU.
"""
function solve_glob!(M)
    Nsys = M.glob.Nsys
    docl = M.oper.givencl ? 1 : 0

    ue = copy(M.glob.U[4, :])
    ds = M.glob.U[2, :]
    uemax = maximum(abs.(ue))
    for i in 1:length(ue)
        ue[i] = max(ue[i], 1e-10 * uemax)
    end

    # inviscid edge velocity
    ueinv = get_ueinv(M)

    # initialize global variable Jacobian
    NN = 4 * Nsys + docl
    alloc_R_V = size(M.glob.R_V) != (NN, NN)
    if alloc_R_V
        M.glob.R_V = zeros(NN, NN)  # NOTE: sparse candidate
    else
        fill!(M.glob.R_V, 0.0)
    end

    # state indices in global system
    Ids = 2:4:4*Nsys  # delta star indices (Julia: every 4th starting at 2)
    Iue = 4:4:4*Nsys  # ue indices

    # assemble residual: BL residual + ue equation
    R = vcat(M.glob.R, ue .- (ueinv .+ M.vsol.ue_m * (ds .* ue)))

    # assemble Jacobian
    M.glob.R_V[1:3*Nsys, 1:4*Nsys] = M.glob.R_U
    I = (3*Nsys+1):4*Nsys
    # ue equation: d/d(ue) = I - ue_m * diag(ds)
    for k in 1:Nsys
        M.glob.R_V[3*Nsys+k, Iue[k]] = 1.0
    end
    M.glob.R_V[I, Iue] .-= M.vsol.ue_m * Diagonal(ds)
    # ue equation: d/d(ds) = -ue_m * diag(ue)
    M.glob.R_V[I, Ids] = -M.vsol.ue_m * Diagonal(ue)

    if docl == 1
        Rcla, Ru_alpha, Rcla_U = clalpha_residual(M)
        R = vcat(R, Rcla)
        M.glob.R_V[I, 4*Nsys+1] = Ru_alpha
        M.glob.R_V[4*Nsys+1, :] = Rcla_U
    end

    # solve system for dV
    dV = -(M.glob.R_V \ R)

    # store dU, reshaped, in M
    M.glob.dU = reshape(dV[1:4*Nsys], 4, Nsys)
    M.glob.dalpha = docl == 1 ? dV[end] : 0.0
    return nothing
end


#-------------------------------------------------------------------------------
"""
    update_state!(M)

Update state with under-relaxation, respecting physical constraints.
"""
function update_state!(M)
    any(imag.(M.glob.U[3, :]) .!= 0) && error("imaginary amp in U")
    any(imag.(M.glob.dU[3, :]) .!= 0) && error("imaginary amp in dU")

    # max ctau
    It = findall(M.vsol.turb)
    ctmax = length(It) > 0 ? maximum(M.glob.U[3, It]) : 0.0

    omega = 1.0

    # limit theta and delta*
    for k in 1:2
        Uk = @view M.glob.U[k, :]
        dUk = @view M.glob.dU[k, :]
        fmin = minimum(dUk ./ Uk)
        if fmin < -0.5
            om = abs(0.5 / fmin)
            if om < omega
                omega = om
                vprint(M.param, 3, @sprintf("  th/ds decrease: omega = %.5f", omega))
            end
        end
    end

    # limit negative amp/ctau
    Uk = M.glob.U[3, :]
    dUk = M.glob.dU[3, :]
    for i in 1:length(Uk)
        if !M.vsol.turb[i] && Uk[i] < 0.2; continue; end
        if M.vsol.turb[i] && Uk[i] < 0.1 * ctmax; continue; end
        if Uk[i] == 0.0 || dUk[i] == 0.0; continue; end
        if Uk[i] + dUk[i] < 0
            om = 0.8 * abs(Uk[i] / dUk[i])
            if om < omega
                omega = om
                vprint(M.param, 3, @sprintf("  neg sa: omega = %.5f", omega))
            end
        end
    end

    # prevent big changes in amp
    Ilam = findall(.!M.vsol.turb)
    if length(Ilam) > 0
        dumax = maximum(abs.(dUk[Ilam]))
        om = dumax > 0 ? abs(2.0 / dumax) : 1.0
        if om < omega
            omega = om
            vprint(M.param, 3, @sprintf("  amp: omega = %.5f", omega))
        end
    end

    # prevent big changes in ctau
    if length(It) > 0
        dumax = maximum(abs.(dUk[It]))
        om = dumax > 0 ? abs(0.05 / dumax) : 1.0
        if om < omega
            omega = om
            vprint(M.param, 3, @sprintf("  ctau: omega = %.5f", omega))
        end
    end

    # prevent large ue changes
    dUk_ue = M.glob.dU[4, :]
    fmax = maximum(abs.(dUk_ue) ./ M.oper.Vinf)
    om = fmax > 0 ? 0.2 / fmax : 1.0
    if om < omega
        omega = om
        vprint(M.param, 3, @sprintf("  ue: omega = %.5f", omega))
    end

    # prevent large alpha changes
    if abs(M.glob.dalpha) > 2
        omega = min(omega, abs(2 / M.glob.dalpha))
    end

    # take the update
    vprint(M.param, 2, @sprintf("  state update: under-relaxation = %.5f", omega))
    M.glob.U .+= omega .* M.glob.dU
    M.oper.alpha += omega * M.glob.dalpha

    # fix bad Hk after the update
    for si in 1:3
        Hkmin = si == 3 ? 1.00005 : 1.02
        Is = M.vsol.Is[si]
        param = build_param(M, si)
        for i in 1:length(Is)
            j = Is[i]
            Uj = M.glob.U[:, j]
            station_param!(M, param, j)
            Hk, Hk_U = get_Hk(Uj, param)
            if Hk < Hkmin
                M.glob.U[2, j] += 2 * (Hkmin - Hk) * M.glob.U[2, j]
            end
        end
    end

    # fix negative ctau after the update
    for i in It
        if M.glob.U[3, i] < 0
            M.glob.U[3, i] = 0.1 * ctmax
        end
    end

    # rebuild inviscid solution if angle of attack changed
    if abs(omega * M.glob.dalpha) > 1e-10
        rebuild_isol!(M)
    end
    return nothing
end


#-------------------------------------------------------------------------------
"""
    solve_coupled!(M)

Solve the coupled inviscid and viscous system using Newton iteration.
"""
function solve_coupled!(M)
    nNewton = M.param.niglob
    M.glob.conv = false
    M.glob.realloc = true
    vprint(M.param, 1, "\n <<< Beginning coupled solver iterations >>>")

    for iNewton in 0:nNewton-1
        vprint(M.param, 2, "Building global system")
        build_glob_sys!(M)

        vprint(M.param, 2, "Calculating force")
        calc_force!(M)

        Rnorm = norm(M.glob.R)
        vprint(M.param, 1, @sprintf("\nNewton iteration %d, Rnorm = %.10e", iNewton, Rnorm))
        if Rnorm < M.param.rtol
            M.glob.conv = true
            break
        end

        vprint(M.param, 2, "Solving global system")
        solve_glob!(M)

        vprint(M.param, 2, "Updating the state")
        update_state!(M)

        M.glob.realloc = false

        vprint(M.param, 2, "Moving stagnation point")
        stagpoint_move!(M)

        vprint(M.param, 2, "Updating transition")
        update_transition!(M)
    end

    if !M.glob.conv
        vprint(M.param, 1, "\n** Global Newton NOT CONVERGED **\n")
    end
    return nothing
end


#-------------------------------------------------------------------------------
"""
    solve_viscous!(M)

Top-level viscous solve: inviscid + BL init + coupled Newton.
"""
function solve_viscous!(M)
    solve_inviscid!(M)
    M.oper.viscous = true
    init_thermo!(M)
    build_wake!(M)
    stagpoint_find!(M)
    identify_surfaces!(M)
    set_wake_gap!(M)
    calc_ue_m!(M)
    init_boundary_layer!(M)
    stagpoint_move!(M)
    solve_coupled!(M)
    calc_force!(M)
    get_distributions!(M)
    return nothing
end


"""
    aseq(M, a1, a2, da)

Alpha sweep: solve viscous at each alpha in a1:da:a2, returning a polar table.
Geometry setup is done once; each successive point warm-starts from the previous.
On convergence failure, records NaN and reverts to the last converged BL state.
"""
function aseq(M::Mfoil, a1::Real, a2::Real, da::Real)
    @assert M.foil.N > 0 "No panels"

    # One-time geometry setup
    solve_inviscid!(M)
    M.oper.viscous = true
    init_thermo!(M)
    build_wake!(M)
    stagpoint_find!(M)
    identify_surfaces!(M)
    set_wake_gap!(M)
    calc_ue_m!(M)

    alphas = collect(a1:da:a2)
    npts = length(alphas)

    # Polar accumulation
    pol_cl = fill(NaN, npts)
    pol_cd = fill(NaN, npts)
    pol_cm = fill(NaN, npts)
    pol_cdf = fill(NaN, npts)
    pol_cdp = fill(NaN, npts)
    pol_xtr_upper = fill(NaN, npts)
    pol_xtr_lower = fill(NaN, npts)

    # Last converged state for revert on failure
    U_last = Matrix{Float64}(undef, 0, 0)
    turb_last = Bool[]

    for (k, alpha) in enumerate(alphas)
        M.oper.alpha = alpha
        M.oper.givencl = false

        # Update gam from gamref for this alpha
        M.isol.gam = M.isol.gamref[:, 1] .* cosd(alpha) .+ M.isol.gamref[:, 2] .* sind(alpha)

        # First point: cold init; subsequent: warm-start
        M.oper.initbl = (k == 1)

        try
            init_boundary_layer!(M)
            stagpoint_move!(M)
            solve_coupled!(M)
            calc_force!(M)
            get_distributions!(M)

            if !M.glob.conv
                error("not converged")
            end

            pol_cl[k] = M.post.cl
            pol_cd[k] = M.post.cd
            pol_cm[k] = M.post.cm
            pol_cdf[k] = M.post.cdf
            pol_cdp[k] = M.post.cdp

            # Transition locations from distributions
            Is_upper = M.vsol.Is[2]
            Is_lower = M.vsol.Is[1]
            pol_xtr_upper[k] = M.vsol.turb[Is_upper[1]] ? 0.0 : M.vsol.xt
            pol_xtr_lower[k] = M.vsol.turb[Is_lower[1]] ? 0.0 : M.vsol.xt

            # Save converged state
            U_last = copy(M.glob.U)
            turb_last = copy(M.vsol.turb)
        catch
            vprint(M.param, 1, "  aseq: alpha=$(alpha) failed, recording NaN")
            # Revert to last converged state if available
            if length(turb_last) > 0
                M.glob.U .= U_last
                M.vsol.turb .= turb_last
            end
        end
    end

    return (alpha=alphas, cl=pol_cl, cd=pol_cd, cm=pol_cm,
            cdf=pol_cdf, cdp=pol_cdp, xtr_upper=pol_xtr_upper, xtr_lower=pol_xtr_lower)
end


"""
    cseq(M, cl1, cl2, dcl)

CL sweep: solve viscous at each cl in cl1:dcl:cl2, returning a polar table.
Geometry setup is done once; each successive point warm-starts from the previous.
On convergence failure, records NaN and reverts to the last converged BL state.
"""
function cseq(M::Mfoil, cl1::Real, cl2::Real, dcl::Real)
    @assert M.foil.N > 0 "No panels"

    # One-time geometry setup
    solve_inviscid!(M)
    M.oper.viscous = true
    init_thermo!(M)
    build_wake!(M)
    stagpoint_find!(M)
    identify_surfaces!(M)
    set_wake_gap!(M)
    calc_ue_m!(M)

    cls = collect(cl1:dcl:cl2)
    npts = length(cls)

    # Polar accumulation
    pol_alpha = fill(NaN, npts)
    pol_cl = fill(NaN, npts)
    pol_cd = fill(NaN, npts)
    pol_cm = fill(NaN, npts)
    pol_cdf = fill(NaN, npts)
    pol_cdp = fill(NaN, npts)
    pol_xtr_upper = fill(NaN, npts)
    pol_xtr_lower = fill(NaN, npts)

    # Last converged state for revert on failure
    U_last = Matrix{Float64}(undef, 0, 0)
    turb_last = Bool[]
    alpha_last = M.oper.alpha

    for (k, cltgt) in enumerate(cls)
        M.oper.givencl = true
        M.oper.cltgt = cltgt

        # Update gam from gamref for current alpha
        alpha = M.oper.alpha
        M.isol.gam = M.isol.gamref[:, 1] .* cosd(alpha) .+ M.isol.gamref[:, 2] .* sind(alpha)

        # First point: cold init; subsequent: warm-start
        M.oper.initbl = (k == 1)

        try
            init_boundary_layer!(M)
            stagpoint_move!(M)
            solve_coupled!(M)
            calc_force!(M)
            get_distributions!(M)

            if !M.glob.conv
                error("not converged")
            end

            pol_alpha[k] = M.oper.alpha
            pol_cl[k] = M.post.cl
            pol_cd[k] = M.post.cd
            pol_cm[k] = M.post.cm
            pol_cdf[k] = M.post.cdf
            pol_cdp[k] = M.post.cdp

            Is_upper = M.vsol.Is[2]
            Is_lower = M.vsol.Is[1]
            pol_xtr_upper[k] = M.vsol.turb[Is_upper[1]] ? 0.0 : M.vsol.xt
            pol_xtr_lower[k] = M.vsol.turb[Is_lower[1]] ? 0.0 : M.vsol.xt

            # Save converged state
            U_last = copy(M.glob.U)
            turb_last = copy(M.vsol.turb)
            alpha_last = M.oper.alpha
        catch
            vprint(M.param, 1, "  cseq: cl=$(cltgt) failed, recording NaN")
            if length(turb_last) > 0
                M.glob.U .= U_last
                M.vsol.turb .= turb_last
                M.oper.alpha = alpha_last
            end
        end
    end

    M.oper.givencl = false  # restore
    return (alpha=pol_alpha, cl=pol_cl, cd=pol_cd, cm=pol_cm,
            cdf=pol_cdf, cdp=pol_cdp, xtr_upper=pol_xtr_upper, xtr_lower=pol_xtr_lower)
end
