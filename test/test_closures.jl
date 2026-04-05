@testset "Closures - Laminar" begin
    param = JFoil.Param()
    param.Vinf = 1.0; param.Minf = 0.0; param.mu0 = 1e-5; param.rho0 = 1.0
    param.turb = false; param.wake = false; param.simi = false

    U = [0.001, 0.003, 5.0, 0.8]

    H, H_U = JFoil.get_H(U)
    @test H ≈ 3.0 atol=1e-12
    @test H_U ≈ [-3000.0, 1000.0, 0.0, 0.0] atol=1e-6

    Hk, Hk_U = JFoil.get_Hk(U, param)
    @test Hk ≈ 3.0 atol=1e-12

    Ret, Ret_U = JFoil.get_Ret(U, param)
    @test Ret ≈ 80.0 atol=1e-6
    @test Ret_U ≈ [80000.0, 0.0, 0.0, 100.0] atol=1e-2

    Hs, _ = JFoil.get_Hs(U, param)
    @test Hs ≈ 1.54687654375 atol=1e-6

    cf, cf_U = JFoil.get_cf(U, param)
    @test cf ≈ 0.0026748046875 atol=1e-8

    cDi, _ = JFoil.get_cDi(U, param)
    @test cDi ≈ 0.002613125 atol=1e-6

    de, _ = JFoil.get_de(U, param)
    @test de ≈ 0.00701 atol=1e-6

    Us, _ = JFoil.get_Us(U, param)
    @test Us ≈ 0.08593758576388894 atol=1e-6

    damp, _ = JFoil.get_damp(U, param)
    @test damp ≈ 7.269229041499859 atol=1e-3
end

@testset "Closures - Turbulent" begin
    param = JFoil.Param()
    param.Vinf = 1.0; param.Minf = 0.0; param.mu0 = 1e-5; param.rho0 = 1.0
    param.turb = true; param.wake = false; param.simi = false

    U = [0.005, 0.010, 0.1, 1.2]

    H, _ = JFoil.get_H(U)
    @test H ≈ 2.0 atol=1e-12

    Hk, _ = JFoil.get_Hk(U, param)
    @test Hk ≈ 2.0 atol=1e-12

    Ret, _ = JFoil.get_Ret(U, param)
    @test Ret ≈ 600.0 atol=1e-6

    Hs, _ = JFoil.get_Hs(U, param)
    @test Hs ≈ 1.6222916666666667 atol=1e-6

    cf, _ = JFoil.get_cf(U, param)
    @test cf ≈ 0.0018751513374696976 atol=1e-8

    cDi, _ = JFoil.get_cDi(U, param)
    @test cDi ≈ 0.00940691542288937 atol=1e-6

    Us, _ = JFoil.get_Us(U, param)
    @test Us ≈ 0.27038194444444447 atol=1e-6

    cteq, _ = JFoil.get_cteq(U, param)
    @test cteq ≈ 0.06231936088984647 atol=1e-6
end

@testset "Closures - Wake" begin
    param = JFoil.Param()
    param.Vinf = 1.0; param.Minf = 0.0; param.mu0 = 1e-5; param.rho0 = 1.0
    param.turb = true; param.wake = true; param.simi = false

    U = [0.01, 0.02, 0.05, 0.9]

    Hw, _ = JFoil.get_Hw(U, 0.005)
    @test Hw ≈ 0.5 atol=1e-12

    cDi, _ = JFoil.get_cDi(U, param)
    @test cDi ≈ 0.004739185200418441 atol=1e-6

    # cf should be 0 in wake
    cf, _ = JFoil.get_cf(U, param)
    @test cf == 0.0
end

@testset "Closures - Misc" begin
    param = JFoil.Param()
    param.Vinf = 1.0; param.Minf = 0.0; param.mu0 = 1e-5; param.rho0 = 1.0
    param.turb = false; param.wake = false

    U = [0.001, 0.003, 5.0, 0.8]

    # upwind with identical states should give 0.5
    upw, _ = JFoil.get_upw(U, U, param)
    @test upw ≈ 0.5 atol=1e-10

    # Mach2 should be 0 for incompressible
    M2, M2_U = JFoil.get_Mach2(U, param)
    @test M2 == 0.0
    @test M2_U == zeros(4)

    # rho should be rho0 for incompressible
    rho, rho_U = JFoil.get_rho(U, param)
    @test rho == 1.0
    @test rho_U == zeros(4)

    # Hss should be 0 for incompressible
    Hss, _ = JFoil.get_Hss(U, param)
    @test Hss == 0.0
end
