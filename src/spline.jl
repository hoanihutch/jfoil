# spline.jl -- spline utilities (Phase 1)

"""
    CubicSpline1D

Piecewise cubic spline with not-a-knot boundary conditions.
Stores knots `x` and coefficients `c` (4 x n_intervals).
p_i(dt) = c[1,i]*dt^3 + c[2,i]*dt^2 + c[3,i]*dt + c[4,i]
where dt = t - x[i].
"""
struct CubicSpline1D
    x::Vector{Float64}        # knot positions [n]
    c::Matrix{Float64}        # coefficients [4 x (n-1)]
end

"""
    CubicSpline1D(x, y)

Construct a not-a-knot cubic spline through points (x[i], y[i]).
Matches scipy.interpolate.CubicSpline default behavior.
"""
function CubicSpline1D(x::AbstractVector, y::AbstractVector)
    n = length(x)
    @assert n >= 2 "Need at least 2 points"
    @assert n == length(y) "x and y must have same length"

    if n == 2
        # Linear case
        dx = x[2] - x[1]
        slope = (y[2] - y[1]) / dx
        c = zeros(4, 1)
        c[3, 1] = slope
        c[4, 1] = y[1]
        return CubicSpline1D(collect(Float64, x), c)
    end

    if n == 3
        # Quadratic (not-a-knot with 3 points = single cubic through all)
        # With not-a-knot: d[1] = d[2], which for 3 points means one polynomial
        dx = diff(collect(Float64, x))
        dy = diff(collect(Float64, y))
        slopes = dy ./ dx

        # For 3 points with not-a-knot, we solve for a single quadratic (d=0 forced)
        # Actually not-a-knot means third derivative is continuous at internal knots.
        # For n=3, that means c[1,1] = c[1,2] (same cubic coeff).
        # Let's just build the full tridiagonal system.
    end

    # General case: n >= 3
    xf = collect(Float64, x)
    yf = collect(Float64, y)
    dx = diff(xf)
    dy = diff(yf)
    slopes = dy ./ dx
    ni = n - 1  # number of intervals

    # Solve for second derivatives (m = s'') at each knot
    # Natural spline system: dx[i]*m[i] + 2*(dx[i]+dx[i+1])*m[i+1] + dx[i+1]*m[i+2] = 6*(slopes[i+1]-slopes[i])
    # Not-a-knot: third derivative continuous at x[2] and x[n-1]
    # i.e. c[1,1]=c[1,2] at left, c[1,n-2]=c[1,n-1] at right

    # We solve for m[1..n] (second derivatives at knots)
    # Interior equations (i=2..n-1):
    # dx[i-1]*m[i-1] + 2*(dx[i-1]+dx[i])*m[i] + dx[i]*m[i+1] = 6*(slopes[i] - slopes[i-1])

    # Set up tridiagonal system for m[1..n]
    # Not-a-knot BCs:
    # Left: d3 continuous at x[2] => (m[2]-m[1])/dx[1] = (m[3]-m[2])/dx[2]
    #   => dx[2]*m[1] - (dx[1]+dx[2])*m[2] + dx[1]*m[3] = 0
    # Right: d3 continuous at x[n-1] => (m[n-1]-m[n-2])/dx[n-2] = (m[n]-m[n-1])/dx[n-1]
    #   => dx[n-1]*m[n-2] - (dx[n-2]+dx[n-1])*m[n-1] + dx[n-2]*m[n] = 0

    # Build tridiagonal system Am = rhs
    A = zeros(n, n)  # NOTE: sparse candidate
    rhs = zeros(n)

    # Row 1: not-a-knot left BC
    A[1, 1] = dx[2]
    A[1, 2] = -(dx[1] + dx[2])
    A[1, 3] = dx[1]
    rhs[1] = 0.0

    # Interior rows
    for i in 2:n-1
        A[i, i-1] = dx[i-1]
        A[i, i] = 2.0 * (dx[i-1] + dx[i])
        A[i, i+1] = dx[i]
        rhs[i] = 6.0 * (slopes[i] - slopes[i-1])
    end

    # Row n: not-a-knot right BC
    A[n, n-2] = dx[ni]
    A[n, n-1] = -(dx[ni-1] + dx[ni])
    A[n, n] = dx[ni-1]
    rhs[n] = 0.0

    m = A \ rhs

    # Build piecewise polynomial coefficients
    # On interval i: p(dt) = a*dt^3 + b*dt^2 + c_coef*dt + d
    # where dt = t - x[i]
    # a = (m[i+1] - m[i]) / (6*dx[i])
    # b = m[i] / 2
    # c_coef = slopes[i] - dx[i]*(2*m[i] + m[i+1])/6
    # d = y[i]
    c = zeros(4, ni)
    for i in 1:ni
        c[1, i] = (m[i+1] - m[i]) / (6.0 * dx[i])
        c[2, i] = m[i] / 2.0
        c[3, i] = slopes[i] - dx[i] * (2.0 * m[i] + m[i+1]) / 6.0
        c[4, i] = yf[i]
    end

    return CubicSpline1D(xf, c)
end

"""
    (cs::CubicSpline1D)(t)

Evaluate the cubic spline at value(s) t.
"""
function (cs::CubicSpline1D)(t::Real)
    # Find interval
    i = _find_interval(cs.x, t)
    dt = t - cs.x[i]
    return cs.c[1, i] * dt^3 + cs.c[2, i] * dt^2 + cs.c[3, i] * dt + cs.c[4, i]
end

function (cs::CubicSpline1D)(t::AbstractVector)
    return [cs(ti) for ti in t]
end

"""
    _find_interval(x, t)

Find interval index i such that x[i] <= t <= x[i+1].
Clamps to valid range.
"""
function _find_interval(x::Vector{Float64}, t::Real)
    n = length(x)
    if t <= x[1]
        return 1
    elseif t >= x[n]
        return n - 1
    end
    # Binary search
    lo, hi = 1, n
    while lo < hi - 1
        mid = (lo + hi) >> 1
        if x[mid] <= t
            lo = mid
        else
            hi = mid
        end
    end
    return lo
end

"""
    derivative(cs::CubicSpline1D)

Return a new CubicSpline1D representing the derivative.
"""
function derivative(cs::CubicSpline1D)
    ni = size(cs.c, 2)
    dc = zeros(4, ni)
    for i in 1:ni
        dc[1, i] = 0.0
        dc[2, i] = 3.0 * cs.c[1, i]
        dc[3, i] = 2.0 * cs.c[2, i]
        dc[4, i] = cs.c[3, i]
    end
    return CubicSpline1D(copy(cs.x), dc)
end


"""
    Spline2D

Two-dimensional parametric spline (x(s), y(s)).
"""
struct Spline2D
    X::CubicSpline1D
    Y::CubicSpline1D
end


"""
    quadseg()

Returns quadrature points and weights for a [0,1] line segment (5-point Gauss-Legendre).
"""
function quadseg()
    x = [0.046910077030668, 0.230765344947158, 0.500000000000000,
         0.769234655052842, 0.953089922969332]
    w = [0.118463442528095, 0.239314335249683, 0.284444444444444,
         0.239314335249683, 0.118463442528095]
    return x, w
end


"""
    spline2d(X)

Spline 2D points with arclength parameterization.
INPUT: X [2 x N] point coordinates
OUTPUT: Spline2D struct
"""
function spline2d(X::AbstractMatrix)
    N = size(X, 2)
    S = zeros(N)
    Snew = zeros(N)

    # estimate arclength
    for i in 2:N
        S[i] = S[i-1] + norm2(@view(X[:, i]) .- @view(X[:, i-1]))
    end
    PPX = CubicSpline1D(S, X[1, :])
    PPY = CubicSpline1D(S, X[2, :])

    # re-integrate to true arclength via several passes
    xq, wq = quadseg()
    for ipass in 1:10
        Snew[1] = S[1]
        for i in 1:N-1
            ds = S[i+1] - S[i]
            st = xq .* ds
            px = PPX.c[:, i]
            xs = 3.0 .* px[1] .* st .* st .+ 2.0 .* px[2] .* st .+ px[3]
            py = PPY.c[:, i]
            ys = 3.0 .* py[1] .* st .* st .+ 2.0 .* py[2] .* st .+ py[3]
            sint = dot(wq, sqrt.(xs .* xs .+ ys .* ys)) * ds
            Snew[i+1] = Snew[i] + sint
        end
        S .= Snew
        PPX = CubicSpline1D(S, X[1, :])
        PPY = CubicSpline1D(S, X[2, :])
    end

    return Spline2D(PPX, PPY)
end


"""
    splineval(PP, S)

Evaluate 2D spline at given S values. Returns [2 x length(S)] matrix.
"""
function splineval(PP::Spline2D, S::AbstractVector)
    return hcat([[PP.X(s), PP.Y(s)] for s in S]...)::Matrix{Float64}
end


"""
    splinetan(PP, S)

Evaluate 2D spline tangent at given S values. Returns [2 x length(S)] matrix.
"""
function splinetan(PP::Spline2D, S::AbstractVector)
    DPX = derivative(PP.X)
    DPY = derivative(PP.Y)
    return hcat([[DPX(s), DPY(s)] for s in S]...)::Matrix{Float64}
end


"""
    spline_curvature(Xin, N, Ufac, TEfac, stgt)

Splines 2D points and samples using curvature-based spacing.
INPUT:
  Xin  : points to spline [2 x npt]
  N    : number of output points
  Ufac : uniformity factor (higher = more uniform)
  TEfac: trailing-edge resolution factor
  stgt : optional target s values, or nothing
OUTPUT:
  X  : new points [2 x N]
  S  : spline s values [N]
  XS : spline tangents [2 x N]
"""
function spline_curvature(Xin::AbstractMatrix, N::Int, Ufac::Real, TEfac::Real, stgt)
    xmin = minimum(@view(Xin[1, :]))
    xmax = maximum(@view(Xin[1, :]))

    # spline given points
    PP = spline2d(Xin)

    # curvature-based spacing on geom
    nfine = 501
    smax = PP.X.x[end]
    s = collect(range(0.0, smax, length=nfine))
    xyfine = splineval(PP, s)
    PPfine = spline2d(xyfine)

    if stgt === nothing
        s = PPfine.X.x
        sk = zeros(nfine)
        xq, wq = quadseg()
        for i in 1:nfine-1
            ds = s[i+1] - s[i]
            st = xq .* ds
            px = PPfine.X.c[:, i]
            xss = 6.0 .* px[1] .* st .+ 2.0 .* px[2]
            py = PPfine.Y.c[:, i]
            yss = 6.0 .* py[1] .* st .+ 2.0 .* py[2]
            skint = 0.01 * Ufac + 0.5 * dot(wq, sqrt.(xss .* xss .+ yss .* yss)) * ds

            # force TE resolution
            xx = (0.5 * (xyfine[1, i] + xyfine[1, i+1]) - xmin) / (xmax - xmin)
            skint = skint + TEfac * 0.5 * exp(-100.0 * (1.0 - xx))

            sk[i+1] = sk[i] + skint
        end

        # offset by fraction of average
        sk .= sk .+ 2.0 * sum(sk) / nfine

        # arclength values at points
        skl = collect(range(minimum(sk), maximum(sk), length=N))
        # cubic interpolation: sk -> s
        cs_interp = CubicSpline1D(sk, s)
        s = [cs_interp(v) for v in skl]
    else
        s = collect(Float64, stgt)
    end

    X = splineval(PPfine, s)
    XS = splinetan(PPfine, s)

    return X, s, XS
end
