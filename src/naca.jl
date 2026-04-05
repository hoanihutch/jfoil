# naca.jl -- NACA coordinate generation (Phase 1)

"""
    naca_points!(M, digits)

Calculate coordinates of a NACA 4-digit airfoil, store in M.geom.
"""
function naca_points!(M, digits::String)
    M.geom.name = "NACA " * digits
    N_pts = 100  # points per side
    te = 1.5     # trailing-edge bunching factor
    f = collect(range(0.0, 1.0, length=N_pts + 1))
    x = @. 1.0 - (te + 1.0) * f * (1.0 - f)^te - (1.0 - f)^(te + 1.0)

    # normalized thickness (finite TE gap with -.1015*x^4)
    t = @. 0.2969 * sqrt(x) - 0.126 * x - 0.3516 * x^2 + 0.2843 * x^3 - 0.1015 * x^4
    tmax = parse(Float64, digits[end-1:end]) * 0.01
    t = t .* (tmax / 0.2)

    ndigits = length(digits)
    if ndigits == 4
        # 4-digit series
        m_val = parse(Float64, string(digits[1])) * 0.01
        p_val = parse(Float64, string(digits[2])) * 0.1
        c = @. m_val / (1.0 - p_val)^2 * ((1.0 - 2.0 * p_val) + 2.0 * p_val * x - x^2)
        for i in eachindex(x)
            if x[i] < p_val
                c[i] = m_val / p_val^2 * (2.0 * p_val * x[i] - x[i]^2)
            end
        end
    elseif ndigits == 5
        n = parse(Int, string(digits[2]))
        valid = digits[1] == '2' && digits[3] == '0' && n > 0 && n < 6
        @assert valid "5-digit NACA must begin with 2X0, X in 1-5"
        mv = [0.058, 0.126, 0.2025, 0.29, 0.391]
        m_val = mv[n]
        cv = [361.4, 51.64, 15.957, 6.643, 3.23]
        cc = cv[n]
        c = @. (cc / 6.0) * (x^3 - 3.0 * m_val * x^2 + m_val^2 * (3.0 - m_val) * x)
        for i in eachindex(x)
            if x[i] > m_val
                c[i] = (cc / 6.0) * m_val^3 * (1.0 - x[i])
            end
        end
    else
        error("Provide 4 or 5 NACA digits")
    end

    zu = c .+ t
    zl = c .- t
    xs = vcat(reverse(x), x[2:end])
    zs = vcat(reverse(zl), zu[2:end])

    M.geom.npoint = length(xs)
    M.geom.xpoint = vcat(xs', zs')
    M.geom.chord = maximum(xs) - minimum(xs)
end
