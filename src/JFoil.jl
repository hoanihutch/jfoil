module JFoil

using LinearAlgebra
using Printf

# Phase 0: Types and utilities
include("types.jl")
include("utils.jl")

# Phase 1: Spline and geometry
include("spline.jl")
include("naca.jl")
include("geometry.jl")

# Phase 2: Panel influence functions
include("panels.jl")

# Phase 3: Inviscid solver
include("thermo.jl")
include("inviscid.jl")

# Exports
export Geom, Panel, Oper, Isol, Vsol, Glob, Post, Param, Mfoil
export vprint, norm2, dist
export CubicSpline1D, Spline2D, quadseg, spline2d, splineval, splinetan, spline_curvature
export naca_points!
export set_coords!, space_geom, TE_info, panel_info, clear_solution!, make_panels!
export panel_linvortex_velocity, panel_linvortex_stream
export panel_constsource_velocity, panel_constsource_stream
export panel_linsource_velocity, panel_linsource_stream
export init_thermo!, get_cp, get_uk
export build_gamma!, get_ueinv, get_ueinvref, stagpoint_find!
export inviscid_velocity, build_wake!, calc_force!, solve_inviscid!, rebuild_isol!

end # module JFoil
