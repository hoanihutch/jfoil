# geometry.jl -- geometry and paneling functions (Phase 1)

"""
    set_coords!(M, X)

Set geometry from coordinate matrix X [2 x npoint] or [npoint x 2].
Coordinates should start/end at trailing edge, CCW ordering.
"""
function set_coords!(M, X::AbstractMatrix)
    if size(X, 1) > size(X, 2)
        X = collect(X')
    end

    # ensure CCW
    A = 0.0
    for i in 1:size(X, 2)-1
        A += (X[1, i+1] - X[1, i]) * (X[2, i+1] + X[2, i])
    end
    if A < 0
        X = X[:, end:-1:1]
    end

    M.geom.npoint = size(X, 2)
    M.geom.xpoint = X
    M.geom.chord = maximum(X[1, :]) - minimum(X[1, :])
end


"""
    space_geom(dx0, L, Np)

Space Np points geometrically from [0, L] with dx0 as first interval.
Returns vector of point locations.
"""
function space_geom(dx0::Real, L::Real, Np::Int)
    @assert Np > 1 "Need at least two points for spacing."
    N = Np - 1
    d = L / dx0
    a = N * (N - 1.0) * (N - 2.0) / 6.0
    b = N * (N - 1.0) / 2.0
    c_val = N - d
    disc = max(b^2 - 4.0 * a * c_val, 0.0)
    r = 1.0 + (-b + sqrt(disc)) / (2.0 * a)
    for k in 1:10
        R = r^N - 1.0 - d * (r - 1.0)
        R_r = N * r^(N - 1) - d
        dr = -R / R_r
        if abs(dr) < 1e-6
            break
        end
        r -= R / R_r
    end
    intervals = dx0 .* r .^ (0:N-1)
    return vcat(0.0, cumsum(intervals))
end


"""
    TE_info(X)

Returns trailing-edge information for airfoil with node coordinates X [2 x N].
Returns: (t, hTE, dtdx, tcp, tdp)
"""
function TE_info(X::AbstractMatrix)
    t1 = X[:, 1] .- X[:, 2]
    t1 = t1 ./ norm2(t1)    # lower tangent vector
    t2 = X[:, end] .- X[:, end-1]
    t2 = t2 ./ norm2(t2)    # upper tangent vector
    t = 0.5 .* (t1 .+ t2)
    t = t ./ norm2(t)         # average tangent; gap bisector
    s = X[:, end] .- X[:, 1]  # lower to upper connector
    hTE = -s[1] * t[2] + s[2] * t[1]  # TE gap
    dtdx = t1[1] * t2[2] - t2[1] * t1[2]  # thickness slope
    p = s ./ norm2(s)  # unit vector along TE panel
    tcp = abs(t[1] * p[2] - t[2] * p[1])
    tdp = dot(t, p)

    return t, hTE, dtdx, tcp, tdp
end


"""
    panel_info(Xj, xi)

Calculate common panel properties (distance, angles).
INPUT:
  Xj : panel endpoint coordinates [2 x 2]
  xi : control point coordinates [2]
OUTPUT:
  t, n, x, z, d, r1, r2, theta1, theta2
"""
function panel_info(Xj::AbstractMatrix, xi::AbstractVector)
    xj1, zj1 = Xj[1, 1], Xj[2, 1]
    xj2, zj2 = Xj[1, 2], Xj[2, 2]

    # panel-aligned tangent and normal vectors
    t = [xj2 - xj1, zj2 - zj1]
    t = t ./ norm2(t)
    n = [-t[2], t[1]]

    # control point relative to (xj1, zj1)
    xz = [xi[1] - xj1, xi[2] - zj1]
    x = dot(xz, t)
    z = dot(xz, n)

    # distances and angles
    d = dist(xj2 - xj1, zj2 - zj1)
    r1 = dist(x, z)
    r2 = dist(x - d, z)
    theta1 = atan(z, x)
    theta2 = atan(z, x - d)

    return t, n, x, z, d, r1, r2, theta1, theta2
end


"""
    clear_solution!(M)

Clear inviscid/viscous solutions by re-initializing structures.
"""
function clear_solution!(M)
    M.isol = Isol()
    M.vsol = Vsol()
    M.glob = Glob()
    M.post = Post()
    M.wake.N = 0
    M.wake.x = zeros(2, 0)
    M.wake.s = Float64[]
    M.wake.t = zeros(2, 0)
end


"""
    make_panels!(M, npanel, stgt)

Place panels on the current airfoil using curvature-based point distribution.
"""
function make_panels!(M, npanel::Int, stgt)
    clear_solution!(M)
    Ufac = 2      # uniformity factor
    TEfac = 0.1   # trailing-edge factor
    M.foil.x, M.foil.s, M.foil.t = spline_curvature(M.geom.xpoint, npanel + 1, Ufac, TEfac, stgt)
    M.foil.N = size(M.foil.x, 2)
end
