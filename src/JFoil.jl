module JFoil

using LinearAlgebra

# Phase 0: Types and utilities
include("types.jl")
include("utils.jl")

# Phase 1: Spline and geometry
include("spline.jl")
include("naca.jl")
include("geometry.jl")

# Exports
export Geom, Panel, Oper, Isol, Vsol, Glob, Post, Param, Mfoil
export vprint, norm2, dist
export CubicSpline1D, Spline2D, quadseg, spline2d, splineval, splinetan, spline_curvature
export naca_points!
export set_coords!, space_geom, TE_info, panel_info, clear_solution!, make_panels!

end # module JFoil
