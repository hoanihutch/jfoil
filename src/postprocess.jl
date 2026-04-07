# postprocess.jl -- post-processing and geometry modifications (Phase 8)
# Translated from mfoil.py post-processing functions

#-------------------------------------------------------------------------------
"""
    get_distributions!(M)

Compute various BL distributions and store in M.post.
"""
function get_distributions!(M)
    @assert size(M.glob.U, 2) > 0 "no global solution"

    # quantities already in global state
    M.post.th = copy(M.glob.U[1, :])   # theta
    M.post.ds = copy(M.glob.U[2, :])   # delta*
    M.post.sa = copy(M.glob.U[3, :])   # amp or ctau
    M.post.ue = [get_uk(M.glob.U[4, i], M.param)[1] for i in 1:M.glob.Nsys]
    M.post.uei = get_ueinv(M)

    # derived viscous quantities
    N = M.glob.Nsys
    cf = zeros(N); Ret = zeros(N); Hk = zeros(N)
    for si in 1:3
        Is = M.vsol.Is[si]
        param = build_param(M, si)
        for ii in 1:length(Is)
            j = Is[ii]
            Uj = M.glob.U[:, j]
            station_param!(M, param, j)
            uk, _ = get_uk(Uj[4], param)
            cfloc, _ = get_cf(Uj, param)
            cf[j] = cfloc * uk^2 / param.Vinf^2
            Ret[j], _ = get_Ret(Uj, param)
            Hk[j], _ = get_Hk(Uj, param)
        end
    end
    M.post.cf = cf
    M.post.Ret = Ret
    M.post.Hk = Hk
    return nothing
end


#-------------------------------------------------------------------------------
"""
    mgeom_flap!(M, xzhinge, eta)

Deploy a flap at hinge location xzhinge with flap angle eta (degrees, positive=down).
"""
function mgeom_flap!(M, xzhinge::AbstractVector, eta::Real)
    X = M.geom.xpoint
    N = size(X, 2)
    xh = xzhinge[1]

    # identify points on flap
    If = findall(X[1, :] .> xh)

    # rotate flap points
    R = [cosd(eta) sind(eta); -sind(eta) cosd(eta)]
    for i in If
        X[:, i] = xzhinge .+ R * (X[:, i] .- xzhinge)
    end

    # remove flap points that ended up to the left of hinge
    Iremove = If[X[1, If] .< xh]
    Ikeep = setdiff(1:N, Iremove)

    M.geom.xpoint = X[:, Ikeep]
    M.geom.npoint = size(M.geom.xpoint, 2)

    # repanel if panels exist
    if M.foil.N > 0
        make_panels!(M, M.foil.N - 1)
    end

    clear_solution!(M)
    return nothing
end


#-------------------------------------------------------------------------------
"""
    mgeom_addcamber!(M, xzcamb)

Add camber to airfoil from given camberline increment coordinates (2 x Nc).
"""
function mgeom_addcamber!(M, xzcamb::AbstractMatrix)
    # ensure xzcamb is 2 x Nc
    if size(xzcamb, 1) > size(xzcamb, 2)
        xzcamb = transpose(xzcamb) |> collect
    end

    X = M.geom.xpoint

    # interpolate camber delta using cubic spline, add to z-coordinates
    spl = CubicSpline1D(xzcamb[1, :], xzcamb[2, :])
    for i in 1:size(X, 2)
        dz = spl(X[1, i])
        X[2, i] += dz
    end

    M.geom.xpoint = X
    M.geom.npoint = size(M.geom.xpoint, 2)

    if M.foil.N > 0
        make_panels!(M, M.foil.N - 1)
    end

    clear_solution!(M)
    return nothing
end


#-------------------------------------------------------------------------------
"""
    mgeom_derotate!(M)

Derotate airfoil about leading edge to make chordline horizontal.
"""
function mgeom_derotate!(M)
    X = M.geom.xpoint
    N = size(X, 2)

    # LE = point with minimum x
    iLE = argmin(X[1, :])
    xLE = X[:, iLE]

    # TE = midpoint of first and last points
    xTE = 0.5 .* (X[:, 1] .+ X[:, N])

    # rotation angle
    theta = atan(xTE[2] - xLE[2], xTE[1] - xLE[1])
    R = [cos(theta) sin(theta); -sin(theta) cos(theta)]
    for i in 1:N
        X[:, i] = xLE .+ R * (X[:, i] .- xLE)
    end

    M.geom.xpoint = X
    M.geom.npoint = size(M.geom.xpoint, 2)

    if M.foil.N > 0
        make_panels!(M, M.foil.N - 1)
    end

    clear_solution!(M)
    return nothing
end


#-------------------------------------------------------------------------------
"""
    plot_cpplus(p, M)

Add Cp plot with annotation text to subplot p.
"""
function plot_cpplus(p, M)
    chord = M.geom.chord
    x = copy(M.foil.x[1, :])
    N = M.foil.N

    if M.oper.viscous
        x = vcat(x, M.wake.x[1, :])
        colors = [:red, :blue, :black]
        for si in 1:3
            Is = M.vsol.Is[si]
            Plots.plot!(p, x[Is], M.post.cp[Is]; color=colors[si], lw=2, label="")
            Plots.plot!(p, x[Is], M.post.cpi[Is]; color=colors[si], ls=:dash, lw=2, label="")
        end
    else
        Plots.plot!(p, x, M.post.cp; color=:blue, lw=2, label="")
    end

    Plots.yflip!(p, true)
    Plots.ylabel!(p, "cp")

    # annotation
    txt = @sprintf("%s\nMa = %.4f\nalpha = %.2f deg\ncl = %.4f\ncm = %.4f\ncd = %.6f",
        M.geom.name, M.oper.Ma, M.oper.alpha, M.post.cl, M.post.cm, M.post.cd)
    if M.oper.viscous
        txt *= @sprintf("\nRe = %.1e\ncdf = %.5f\ncdp = %.5f",
            M.oper.Re, M.post.cdf, M.post.cdp)
    end
    Plots.annotate!(p, 0.85 * chord, minimum(M.post.cp[1:N]) + 0.3,
        Plots.text(txt, 8, :left))
    return p
end


#-------------------------------------------------------------------------------
"""
    plot_airfoil(p, M)

Add airfoil (and wake) geometry to subplot p.
"""
function plot_airfoil(p, M)
    xz = copy(M.foil.x)
    if M.oper.viscous
        xz = hcat(xz, M.wake.x)
    end
    Plots.plot!(p, xz[1, :], xz[2, :]; color=:black, lw=1, label="", aspect_ratio=:equal)
    return p
end


#-------------------------------------------------------------------------------
"""
    plot_boundary_layer(p, M)

Add boundary layer displacement thickness to subplot p.
"""
function plot_boundary_layer(p, M)
    if !M.oper.viscous; return p; end
    xz = hcat(M.foil.x, M.wake.x)
    N = M.foil.N
    ds = M.post.ds
    rl = 0.5 * (1.0 + (ds[1] - ds[N]) / ds[N + 1])
    ru = 1.0 - rl
    t = hcat(M.foil.t, M.wake.t)
    n = vcat(-t[2:2, :], t[1:1, :])  # outward normals
    for i in 1:size(n, 2)
        nn = sqrt(n[1, i]^2 + n[2, i]^2)
        if nn > 0; n[:, i] ./= nn; end
    end
    xzd = xz .+ n .* ds'
    colors = [:red, :blue, :black]
    for i in 0:3
        si = i
        if si == 2; xzd = xz .+ n .* (ds .* ru)'; end
        if si == 3; xzd = xz .- n .* (ds .* rl)'; si = 2; end
        Is = M.vsol.Is[si + 1]  # Julia 1-based
        Plots.plot!(p, xzd[1, Is], xzd[2, Is]; color=colors[si + 1], lw=2, label="")
    end
    return p
end


#-------------------------------------------------------------------------------
"""
    plot_results(M)

Make a summary results plot with Cp, airfoil shape, and BL thickness.
Requires Plots.jl.
"""
function plot_results(M)
    @assert length(M.post.cp) > 0 "no cp for results plot"

    p1 = Plots.plot(; title="", legend=false)
    plot_cpplus(p1, M)

    p2 = Plots.plot(; title="", legend=false)
    plot_airfoil(p2, M)
    plot_boundary_layer(p2, M)

    fig = Plots.plot(p1, p2; layout=(2, 1), size=(800, 700))
    Plots.display(fig)
    return fig
end


#-------------------------------------------------------------------------------
"""
    check_ping(ep, v, v_u, sname) -> (E, rate)

Check convergence of finite-difference derivative pinging.
- `v`: vector of 3 function evaluations at 0, +ep, +2*ep
- `v_u`: vector of 3 directional derivative evaluations
- Returns error values E and convergence rate (expect ~2 for 2nd order).
"""
function check_ping(ep::Real, v::Vector, v_u::Vector, sname::String)
    E = zeros(2)
    for i in 1:2
        E[i] = norm((v[1+i] .- v[1]) ./ (ep * i) .- 0.5 .* (v_u[1] .+ v_u[1+i]))
    end
    rate = log2(E[2] / E[1])
    @printf("%s ping error convergence rate = %.4f\n", sname, rate)
    return E, rate
end


#-------------------------------------------------------------------------------
"""
    ping_test!(M)

Check derivatives of various functions via finite-difference pinging.
Prints convergence rates (2 = second-order expected).
"""
function ping_test!(M)
    rng = Random.MersenneTwister(17)

    M.oper.alpha = 3.0
    M.oper.Ma = 0.4
    M.oper.viscous = true
    M.param.verb = 2

    # freestream Reynolds numbers
    Rev = [2e3, 1e5]

    # laminar/turbulent test states: th, ds, sa, ue
    Uv = [[0.01, 0.02, 8.4, 0.9], [0.023, 0.05, 0.031, 1.1]]

    # functions to test (scalar closure functions)
    fv = [get_Hk, get_Ret, get_cf, get_cDi, get_Hs, get_Us,
          get_cDi_turbwall, get_cDi_lam, get_cDi_lamwake, get_cDi_outer,
          get_cDi_lamstress, get_cteq, get_cttr, get_de, get_damp,
          get_Mach2, get_Hss, residual_station]

    sturb = ["lam", "turb", "wake"]
    for iRe in 1:length(Rev)
        M.oper.Re = Rev[iRe]
        init_thermo!(M)
        param = build_param(M, 1)
        for it in 1:3  # lam, turb, wake
            param.turb = (it > 1)
            param.wake = (it == 3)
            for ih in 1:length(fv)
                f = fv[ih]
                U = copy(Uv[min(it, 2)])
                srates = ""; smark = ""; serr = ""
                if f === residual_station
                    U = vcat(U, U .* [1.1, 0.8, 0.9, 1.2])
                end
                for k in 1:length(U)
                    ep = 1e-2 * U[k]
                    E = zeros(2)
                    if f === residual_station
                        xi = [0.7, 0.8]; Aux = [0.002, 0.0018]; dx = [-0.2, 0.3]
                        Umat = hcat(U[1:4], U[5:8])
                        v0, v_U0, v_x0 = f(param, xi, Umat, Aux)
                        for iep in 1:2
                            U[k] += ep; xi .+= ep .* dx
                            Umat = hcat(U[1:4], U[5:8])
                            v1, v_U1, v_x1 = f(param, xi, Umat, Aux)
                            U[k] -= ep; xi .-= ep .* dx
                            E[iep] = norm((v1 .- v0) ./ ep .- 0.5 .* (v_U1[:, k] .+ v_U0[:, k] .+ (v_x0 .+ v_x1) * dx))
                            ep /= 2
                        end
                    else
                        v0, v_U0 = f(U, param)
                        for iep in 1:2
                            U[k] += ep
                            v1, v_U1 = f(U, param)
                            U[k] -= ep
                            E[iep] = abs((v1 - v0) / ep - 0.5 * (v_U1[k] + v_U0[k]))
                            ep /= 2
                        end
                    end
                    srate = " N/A"
                    if E[1] > 5e-11 && E[2] > 5e-11
                        m = log2(E[1] / E[2])
                        srate = @sprintf(" %4.1f", m)
                        if m < 1.5; smark = "<==="; end
                    end
                    srates *= " " * srate
                    serr *= @sprintf(" %.2e->%.2e", E[1], E[2])
                end
                fname = f === residual_station ? "residual_station" : string(nameof(f))
                vprint(param, 0, @sprintf("%-18s %-5s err=[%s]  rates=[%s] %s", fname, sturb[it], serr, srates, smark))
            end
        end
    end

    # transition residual ping
    M.oper.Re = 2e6; init_thermo!(M); param = build_param(M, 1)
    U = [0.01 0.013; 0.02 0.023; 8.95 0.028; 0.9 0.85]
    x = [0.7, 0.8]; Aux = [0.0, 0.0]
    dU = rand(rng, 4, 2); dx = rand(rng, 2); ep = 1e-4
    v = Vector[]; v_u = Vector[]
    for ie in 1:3
        R, R_U, R_x = residual_transition(M, param, x, U, Aux)
        push!(v, R); push!(v_u, R_U * vec(dU) .+ R_x * dx)
        U .+= ep .* dU; x .+= ep .* dx
    end
    check_ping(ep, v, v_u, "transition residual")

    # stagnation residual ping
    M.oper.Re = 1e6; M.oper.alpha = 1.0; init_thermo!(M); param = build_param(M, 1)
    U_stag = [0.00004673616, 0.000104289, 0.0, 0.11977917547]
    x_stag = 4.590816441485401e-05
    dU_stag = rand(rng, 4); dx_stag = rand(rng, 1); ep = 1e-6
    v = Vector[]; v_u = Vector[]
    for ie in 1:3
        param.simi = true
        R, R_U, R_x = residual_station(param, [x_stag, x_stag],
            hcat(U_stag, U_stag), [0.0, 0.0])
        param.simi = false
        push!(v, R)
        push!(v_u, (R_U[:, 1:4] .+ R_U[:, 5:8]) * dU_stag .+ (R_x[:, 1] .+ R_x[:, 2]) .* dx_stag[1])
        U_stag .+= ep .* dU_stag; x_stag += ep * dx_stag[1]
    end
    check_ping(ep, v, v_u, "stagnation residual")

    # need a viscous solution for the next tests
    solve_viscous!(M)

    # entire system ping
    Nsys = size(M.glob.U, 2)
    dU_glob = rand(rng, 4, Nsys); dx_glob = 0.1 .* rand(rng, Nsys); ep = 1e-6
    for ix in 0:1
        if ix == 1
            dx_glob .= 0.0
            stagpoint_move!(M)
        end
        v = Vector[]; v_u = Vector[]
        for ie in 1:3
            build_glob_sys!(M)
            push!(v, copy(M.glob.R))
            push!(v_u, M.glob.R_U * vec(dU_glob) .+ M.glob.R_x * dx_glob)
            M.glob.U .+= ep .* dU_glob; M.isol.xi .+= ep .* dx_glob
            if ix == 1; stagpoint_move!(M); end
        end
        M.glob.U .-= 3 * ep .* dU_glob; M.isol.xi .-= 3 * ep .* dx_glob
        check_ping(ep, v, v_u, @sprintf("global system, ix=%d", ix))
    end

    return nothing
end
