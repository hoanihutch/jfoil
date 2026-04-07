module JFoil

using LinearAlgebra
using Printf
using Random
using Plots

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

# Phase 4: Closure relations
include("closures_shape.jl")
include("closures_friction.jl")
include("closures_dissipation.jl")

# Phase 5: Boundary layer
include("boundary_layer.jl")

# Phase 6: Coupling
include("coupling.jl")

# Phase 7: Coupled solver
include("solver.jl")

# Phase 8: Post-processing
include("postprocess.jl")

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
export get_H, get_Hw, get_Mach2, get_Hk, get_Hs, get_Hss, get_Ret, get_rho, get_de
export get_cf, get_cfxt, get_cfutstag, get_Us, get_upw, upwind, get_uq, get_damp
export get_cDi, get_cDi_turbwall, get_cDi_lam, get_cDi_lamwake, get_cDi_outer
export get_cDi_lamstress, get_cDixt, get_cdutstag, get_cteq, get_cttr
export build_param, station_param!
export stagnation_state, thwaites_init
export residual_station, residual_transition
export store_transition!, march_amplification!, update_transition!
export identify_surfaces!, set_wake_gap!, calc_ue_m!, rebuild_ue_m!, wake_sys, wake_init
export init_boundary_layer!, stagpoint_move!, build_glob_sys!, clalpha_residual
export solve_glob!, update_state!, solve_coupled!, solve_viscous!
export get_distributions!, mgeom_flap!, mgeom_addcamber!, mgeom_derotate!
export check_ping, ping_test!
export plot_results, plot_cpplus, plot_airfoil, plot_boundary_layer

end # module JFoil
