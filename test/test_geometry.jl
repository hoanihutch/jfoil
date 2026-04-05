@testset "set_coords!" begin
    M = JFoil.Mfoil()
    X = [0.0 0.5 1.0 0.5 0.0;
         0.0 0.1 0.0 -0.1 0.0]
    JFoil.set_coords!(M, X)
    @test M.geom.npoint == 5
    @test M.geom.chord ≈ 1.0
    @test size(M.geom.xpoint) == (2, 5)
end

@testset "space_geom" begin
    sg = JFoil.space_geom(0.01, 1.0, 20)
    @test length(sg) == 20
    @test sg[1] ≈ 0.0
    @test sg[2] ≈ 0.01 atol=1e-8
    @test sg[end] ≈ 1.0 atol=0.01  # Python gives 1.00029
    # Monotonically increasing
    for i in 2:length(sg)
        @test sg[i] > sg[i-1]
    end
end

@testset "TE_info" begin
    M = JFoil.Mfoil()
    JFoil.naca_points!(M, "0012")
    t, hTE, dtdx, tcp, tdp = JFoil.TE_info(M.geom.xpoint)
    @test t ≈ [1.0, 0.0] atol=1e-10
    @test hTE ≈ 0.00252 atol=1e-4
    @test dtdx ≈ -0.27479566541781286 atol=1e-4
    @test tcp ≈ 1.0 atol=1e-10
    @test tdp ≈ 0.0 atol=1e-10
end

@testset "panel_info" begin
    Xj = [0.0 1.0; 0.0 0.0]
    xi = [0.5, 0.1]
    t, n, x, z, d, r1, r2, theta1, theta2 = JFoil.panel_info(Xj, xi)
    @test t ≈ [1.0, 0.0] atol=1e-12
    @test n ≈ [0.0, 1.0] atol=1e-12
    @test x ≈ 0.5 atol=1e-12
    @test z ≈ 0.1 atol=1e-12
    @test d ≈ 1.0 atol=1e-12
    @test r1 ≈ 0.5099019513592785 atol=1e-10
    @test r2 ≈ 0.5099019513592785 atol=1e-10
    @test theta1 ≈ 0.19739555984988078 atol=1e-10
    @test theta2 ≈ 2.9441970937399127 atol=1e-10
end

@testset "clear_solution!" begin
    M = JFoil.Mfoil()
    M.isol.sstag = 5.0
    JFoil.clear_solution!(M)
    @test M.isol.sstag == 0.0
    @test M.wake.N == 0
end

@testset "make_panels!" begin
    M = JFoil.Mfoil()
    JFoil.naca_points!(M, "0012")
    JFoil.make_panels!(M, 199, nothing)
    @test M.foil.N == 200
    @test size(M.foil.x) == (2, 200)
    @test length(M.foil.s) == 200
    @test size(M.foil.t) == (2, 200)
    @test M.foil.s[1] ≈ 0.0 atol=1e-10
    # TE points should be near (1, ±0.00126)
    @test M.foil.x[1, 1] ≈ 1.0 atol=0.01
    @test M.foil.x[1, end] ≈ 1.0 atol=0.01
    # Python reference: s[-1] ≈ 2.039271378057922
    @test M.foil.s[end] ≈ 2.039271378057922 atol=0.01

    # NACA 2412 with 99 panels
    M2 = JFoil.Mfoil()
    JFoil.naca_points!(M2, "2412")
    JFoil.make_panels!(M2, 99, nothing)
    @test M2.foil.N == 100
end
