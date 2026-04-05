module JFoil

using LinearAlgebra

# Phase 0: Types and utilities
include("types.jl")
include("utils.jl")

# Exports
export Geom, Panel, Oper, Isol, Vsol, Glob, Post, Param, Mfoil
export vprint, norm2, dist

end # module JFoil
