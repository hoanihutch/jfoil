# panels.jl -- panel influence functions (Phase 2)

"""
    panel_linvortex_velocity(Xj, xi, vdir, onmid)

Linear vortex panel velocity influence coefficients.
Returns (a, b) where velocity = a*g1 + b*g2.
If vdir is nothing, a,b are 2-vectors; otherwise scalars (dotted with vdir).
"""
@inline function panel_linvortex_velocity(Xj::AbstractMatrix, xi::AbstractVector, vdir, onmid::Bool)
    t, n, x, z, d, r1, r2, theta1, theta2 = panel_info(Xj, xi)

    if onmid
        ug1 = 0.25
        ug2 = 0.25
        wg1 = -1.0 / (2.0 * π)
        wg2 = 1.0 / (2.0 * π)
    else
        temp1 = (theta2 - theta1) / (2.0 * π)
        temp2 = (2.0 * z * log(r1 / r2) - 2.0 * x * (theta2 - theta1)) / (4.0 * π * d)
        ug1 = temp1 + temp2
        ug2 = -temp2
        temp1 = log(r2 / r1) / (2.0 * π)
        temp2 = (x * log(r1 / r2) - d + z * (theta2 - theta1)) / (2.0 * π * d)
        wg1 = temp1 + temp2
        wg2 = -temp2
    end

    a = SVector(ug1 * t[1] + wg1 * n[1], ug1 * t[2] + wg1 * n[2])
    b = SVector(ug2 * t[1] + wg2 * n[1], ug2 * t[2] + wg2 * n[2])

    if vdir !== nothing
        a = dot(a, vdir)
        b = dot(b, vdir)
    end

    return a, b
end


"""
    panel_linvortex_stream(Xj, xi)

Linear vortex panel streamfunction influence coefficients.
Returns (a, b) where psi = a*g1 + b*g2.
"""
@inline function panel_linvortex_stream(Xj::AbstractMatrix, xi::AbstractVector)
    t, n, x, z, d, r1, r2, theta1, theta2 = panel_info(Xj, xi)

    ep = 1e-9
    logr1 = r1 < ep ? 0.0 : log(r1)
    logr2 = r2 < ep ? 0.0 : log(r2)

    P1 = (0.5 / π) * (z * (theta2 - theta1) - d + x * logr1 - (x - d) * logr2)
    P2 = x * P1 + (0.5 / π) * (0.5 * r2^2 * logr2 - 0.5 * r1^2 * logr1 - r2^2 / 4.0 + r1^2 / 4.0)

    a = P1 - P2 / d
    b = P2 / d

    return a, b
end


"""
    panel_constsource_velocity(Xj, xi, vdir)

Constant source panel velocity influence coefficient.
Returns a where velocity = a*sigma.
"""
@inline function panel_constsource_velocity(Xj::AbstractMatrix, xi::AbstractVector, vdir)
    t, n, x, z, d, r1, r2, theta1, theta2 = panel_info(Xj, xi)

    ep = 1e-9
    if r1 < ep
        logr1 = 0.0; theta1 = π; theta2 = π
    else
        logr1 = log(r1)
    end
    if r2 < ep
        logr2 = 0.0; theta1 = 0.0; theta2 = 0.0
    else
        logr2 = log(r2)
    end

    u = (0.5 / π) * (logr1 - logr2)
    w = (0.5 / π) * (theta2 - theta1)

    a = SVector(u * t[1] + w * n[1], u * t[2] + w * n[2])
    if vdir !== nothing
        a = dot(a, vdir)
    end

    return a
end


"""
    panel_constsource_stream(Xj, xi)

Constant source panel streamfunction influence coefficient.
Returns a where psi = a*sigma.
"""
@inline function panel_constsource_stream(Xj::AbstractMatrix, xi::AbstractVector)
    t, n, x, z, d, r1, r2, theta1, theta2 = panel_info(Xj, xi)

    ep = 1e-9
    if r1 < ep
        logr1 = 0.0; theta1 = π; theta2 = π
    else
        logr1 = log(r1)
    end
    if r2 < ep
        logr2 = 0.0; theta1 = 0.0; theta2 = 0.0
    else
        logr2 = log(r2)
    end

    P = (x * (theta1 - theta2) + d * theta2 + z * logr1 - z * logr2) / (2.0 * π)

    dP = d
    P = (theta1 + theta2) > π ? (P - 0.25 * dP) : (P + 0.75 * dP)

    return P
end


"""
    panel_linsource_velocity(Xj, xi, vdir)

Linear source panel velocity influence coefficients.
Returns (a, b) where velocity = a*s1 + b*s2.
"""
@inline function panel_linsource_velocity(Xj::AbstractMatrix, xi::AbstractVector, vdir)
    t, n, x, z, d, r1, r2, theta1, theta2 = panel_info(Xj, xi)

    temp1 = log(r1 / r2) / (2.0 * π)
    temp2 = (x * log(r1 / r2) - d + z * (theta2 - theta1)) / (2.0 * π * d)
    ug1 = temp1 - temp2
    ug2 = temp2
    temp1 = (theta2 - theta1) / (2.0 * π)
    temp2 = (-z * log(r1 / r2) + x * (theta2 - theta1)) / (2.0 * π * d)
    wg1 = temp1 - temp2
    wg2 = temp2

    a = SVector(ug1 * t[1] + wg1 * n[1], ug1 * t[2] + wg1 * n[2])
    b = SVector(ug2 * t[1] + wg2 * n[1], ug2 * t[2] + wg2 * n[2])

    if vdir !== nothing
        a = dot(a, vdir)
        b = dot(b, vdir)
    end

    return a, b
end


"""
    panel_linsource_stream(Xj, xi)

Linear source panel streamfunction influence coefficients.
Returns (a, b) where psi = a*s1 + b*s2.
"""
@inline function panel_linsource_stream(Xj::AbstractMatrix, xi::AbstractVector)
    t, n, x, z, d, r1, r2, theta1, theta2 = panel_info(Xj, xi)

    # make branch cut at theta = 0
    if theta1 < 0
        theta1 += 2.0 * π
    end
    if theta2 < 0
        theta2 += 2.0 * π
    end

    ep = 1e-9
    if r1 < ep
        logr1 = 0.0; theta1 = π; theta2 = π
    else
        logr1 = log(r1)
    end
    if r2 < ep
        logr2 = 0.0; theta1 = 0.0; theta2 = 0.0
    else
        logr2 = log(r2)
    end

    P1 = (0.5 / π) * (x * (theta1 - theta2) + theta2 * d + z * logr1 - z * logr2)
    P2 = x * P1 + (0.5 / π) * (0.5 * r2^2 * theta2 - 0.5 * r1^2 * theta1 - 0.5 * z * d)

    a = P1 - P2 / d
    b = P2 / d

    return a, b
end
