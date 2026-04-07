using JFoil
using Test

@testset "JFoil.jl" begin
    include("test_utils.jl")
    include("test_spline.jl")
    include("test_naca.jl")
    include("test_geometry.jl")
    include("test_panels.jl")
    include("test_thermo.jl")
    include("test_inviscid.jl")
    include("test_closures.jl")
    include("test_boundary_layer.jl")
    include("test_coupling.jl")
    include("test_solver.jl")
    include("test_postprocess.jl")
end
