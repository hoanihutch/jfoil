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
