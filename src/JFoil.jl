module JFoil

using LinearAlgebra

# Phase 0: Types and utilities
include("types.jl")
include("utils.jl")

# Phase 1: Spline and geometry
include("spline.jl")

# Exports
export Geom, Panel, Oper, Isol, Vsol, Glob, Post, Param, Mfoil
export vprint, norm2, dist
export CubicSpline1D, Spline2D, quadseg, spline2d, splineval, splinetan, spline_curvature

end # module JFoil
