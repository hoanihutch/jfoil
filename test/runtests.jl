using JFoil
using Test

@testset "JFoil.jl" begin
    include("test_utils.jl")
    include("test_spline.jl")
    include("test_naca.jl")
end
