# types.jl -- mutable structs for JFoil (Phase 0)
# Translated from mfoil.py data structures

#-------------------------------------------------------------------------------
mutable struct Geom
    chord::Float64
    wakelen::Float64
    npoint::Int
    name::String
    xpoint::Matrix{Float64}   # [2 x npoint]
    xref::Vector{Float64}     # moment reference point [2]
end

function Geom()
    Geom(1.0, 1.0, 1, "noname", zeros(2, 0), [0.25, 0.0])
end

#-------------------------------------------------------------------------------
mutable struct Panel
    N::Int
    x::Matrix{Float64}    # [2 x N]
    s::Vector{Float64}    # arclength at nodes
    t::Matrix{Float64}    # dx/ds, dy/ds tangents [2 x N]
end

function Panel()
    Panel(0, zeros(2, 0), Float64[], zeros(2, 0))
end

#-------------------------------------------------------------------------------
mutable struct Oper
    Vinf::Float64
    alpha::Float64
    rho::Float64
    cltgt::Float64
    givencl::Bool
    initbl::Bool
    viscous::Bool
    redowake::Bool
    Re::Float64
    Ma::Float64
end

function Oper()
    Oper(1.0, 0.0, 1.0, 0.0, false, true, false, false, 1e5, 0.0)
end

#-------------------------------------------------------------------------------
mutable struct Isol
    AIC::Matrix{Float64}
    gamref::Matrix{Float64}    # [N x 2]
    gam::Vector{Float64}
    sstag::Float64
    sstag_g::Vector{Float64}   # [2]
    sstag_ue::Vector{Float64}  # [2]
    Istag::Vector{Int}         # [2]
    sgnue::Vector{Float64}
    xi::Vector{Float64}
    uewi::Vector{Float64}
    uewiref::Matrix{Float64}   # wake ue at 0,90 deg
    xstag::Vector{Float64}     # [2] physical stagnation location
end

function Isol()
    Isol(
        zeros(0, 0), zeros(0, 0), Float64[],
        0.0, [0.0, 0.0], [0.0, 0.0], [0, 0],
        Float64[], Float64[], Float64[], zeros(0, 0), zeros(2)
    )
end

#-------------------------------------------------------------------------------
mutable struct Vsol
    th::Vector{Float64}
    ds::Vector{Float64}
    Is::Vector{Vector{Int}}    # 3 arrays of surface indices
    wgap::Vector{Float64}
    ue_m::Matrix{Float64}      # NOTE: sparse candidate
    sigma_m::Matrix{Float64}   # NOTE: sparse candidate
    ue_sigma::Matrix{Float64}  # NOTE: sparse candidate
    turb::Vector{Bool}
    xt::Float64
    Xt::Matrix{Float64}        # [2 x 2]
end

function Vsol()
    Vsol(
        Float64[], Float64[], Vector{Int}[], Float64[],
        zeros(0, 0), zeros(0, 0), zeros(0, 0),
        Bool[], 0.0, zeros(2, 2)
    )
end

#-------------------------------------------------------------------------------
mutable struct Glob
    Nsys::Int
    U::Matrix{Float64}      # [4 x Nsys]
    dU::Matrix{Float64}
    dalpha::Float64
    conv::Bool
    R::Vector{Float64}       # [3*Nsys]
    R_U::Matrix{Float64}     # [3*Nsys x 4*Nsys]  NOTE: sparse candidate
    R_x::Matrix{Float64}     # [3*Nsys x Nsys]    NOTE: sparse candidate
    R_V::Matrix{Float64}     # [4*Nsys x 4*Nsys]  NOTE: sparse candidate
    realloc::Bool
end

function Glob()
    Glob(
        0, zeros(4, 0), zeros(4, 0), 0.0, true,
        Float64[], zeros(0, 0), zeros(0, 0), zeros(0, 0), false
    )
end

#-------------------------------------------------------------------------------
mutable struct Post
    cp::Vector{Float64}
    cpi::Vector{Float64}
    cl::Float64
    cl_ue::Vector{Float64}
    cl_alpha::Float64
    cm::Float64
    cdpi::Float64
    cd::Float64
    cd_U::Vector{Float64}    # [4]
    cdf::Float64
    cdp::Float64

    # distributions
    th::Vector{Float64}
    ds::Vector{Float64}
    sa::Vector{Float64}
    ue::Vector{Float64}
    uei::Vector{Float64}
    cf::Vector{Float64}
    Ret::Vector{Float64}
    Hk::Vector{Float64}
end

function Post()
    Post(
        Float64[], Float64[], 0.0, Float64[], 0.0, 0.0, 0.0, 0.0, Float64[], 0.0, 0.0,
        Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]
    )
end

#-------------------------------------------------------------------------------
mutable struct Param
    verb::Int
    rtol::Float64
    niglob::Int
    doplot::Bool

    # viscous parameters
    ncrit::Float64
    Cuq::Float64
    Dlr::Float64
    SlagK::Float64

    # initial Ctau after transition
    CtauC::Float64
    CtauE::Float64

    # G-beta locus
    GA::Float64
    GB::Float64
    GC::Float64

    # operating conditions and thermodynamics
    Minf::Float64
    Vinf::Float64
    muinf::Float64
    mu0::Float64
    rho0::Float64
    H0::Float64
    Tsrat::Float64
    gam::Float64
    KTb::Float64
    KTl::Float64
    cps::Float64

    # station information
    simi::Bool
    turb::Bool
    wake::Bool
end

function Param()
    Param(
        1, 1e-10, 50, true,
        9.0, 1.0, 0.9, 5.6,
        1.8, 3.3,
        6.7, 0.75, 18.0,
        0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.35, 1.4, 1.0, 0.0, 0.0,
        false, false, false
    )
end

"""
    copy_param(p::Param)

Lightweight field-by-field copy of Param (all scalar fields, no heap refs).
Replaces deepcopy(param) which is expensive for a 28-field mutable struct.
"""
function copy_param(p::Param)
    new_p = Param()
    for f in fieldnames(Param)
        setfield!(new_p, f, getfield(p, f))
    end
    return new_p
end

#-------------------------------------------------------------------------------
mutable struct Mfoil
    version::String
    geom::Geom
    foil::Panel
    wake::Panel
    oper::Oper
    isol::Isol
    vsol::Vsol
    glob::Glob
    post::Post
    param::Param
end

function Mfoil()
    Mfoil(
        "2022-02-22",
        Geom(), Panel(), Panel(), Oper(), Isol(), Vsol(), Glob(), Post(), Param()
    )
end
