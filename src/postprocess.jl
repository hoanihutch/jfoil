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
