@testset "init_thermo!" begin
    M = JFoil.Mfoil()
    JFoil.naca_points!(M, "0012")
    JFoil.make_panels!(M, 199, nothing)
    M.oper.Ma = 0.0
    M.oper.Re = 1e5
    JFoil.init_thermo!(M)
    @test M.param.Vinf == 1.0
    @test M.param.Minf == 0.0
    @test M.param.KTb == 1.0  # default, not modified for Ma=0
    @test M.param.mu0 == M.param.muinf  # finf=1 for Ma=0
end

@testset "get_cp" begin
    param = JFoil.Param()
    param.Vinf = 1.0
    param.Minf = 0.0

    # Incompressible
    cp, cp_u = JFoil.get_cp(0.5, param)
    @test cp ≈ 1.0 - 0.25 atol=1e-12
    @test cp_u ≈ -1.0 atol=1e-12

    cp0, cp0_u = JFoil.get_cp(1.0, param)
    @test cp0 ≈ 0.0 atol=1e-12

    # Vectorized
    cpv, cpv_u = JFoil.get_cp([0.5, 1.0], param)
    @test cpv[1] ≈ 0.75 atol=1e-12
    @test cpv[2] ≈ 0.0 atol=1e-12
end

@testset "get_uk" begin
    param = JFoil.Param()
    param.Vinf = 1.0
    param.Minf = 0.0

    # Incompressible: uk = u
    uk, uk_u = JFoil.get_uk(0.5, param)
    @test uk ≈ 0.5 atol=1e-12
    @test uk_u ≈ 1.0 atol=1e-12
end
