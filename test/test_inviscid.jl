@testset "solve_inviscid! NACA 0012 alpha=2" begin
    M = JFoil.Mfoil()
    JFoil.naca_points!(M, "0012")
    JFoil.make_panels!(M, 199, nothing)
    M.oper.alpha = 2.0
    M.param.doplot = false
    M.param.verb = 0  # quiet
    JFoil.solve_inviscid!(M)

    # Python reference
    @test M.post.cl ≈ 0.24167575662386903 atol=1e-4
    @test M.post.cm ≈ -0.002816308048432195 atol=1e-4
    @test M.post.cdpi ≈ -0.001105166929456156 atol=1e-4

    # Gamma values
    @test M.isol.gam[1] ≈ -0.75662476 atol=1e-3
    @test M.isol.gam[end] ≈ 0.75662476 atol=1e-3

    # gamref should be N x 2
    @test size(M.isol.gamref) == (200, 2)

    # Convergence flag
    @test M.glob.conv == true
end

@testset "solve_inviscid! NACA 2412 alpha=0" begin
    M = JFoil.Mfoil()
    JFoil.naca_points!(M, "2412")
    JFoil.make_panels!(M, 99, nothing)
    M.oper.alpha = 0.0
    M.param.doplot = false
    M.param.verb = 0
    JFoil.solve_inviscid!(M)

    @test M.post.cl ≈ 0.25532498405914295 atol=1e-3
    @test M.post.cm ≈ -0.055736685131018496 atol=1e-3
end

@testset "get_ueinv" begin
    M = JFoil.Mfoil()
    JFoil.naca_points!(M, "0012")
    JFoil.make_panels!(M, 199, nothing)
    M.oper.alpha = 2.0
    M.param.verb = 0
    JFoil.solve_inviscid!(M)

    ueinv = JFoil.get_ueinv(M)
    @test length(ueinv) == M.foil.N
    # At TE, ue should be moderate
    @test abs(ueinv[1]) > 0
end

@testset "inviscid_velocity" begin
    M = JFoil.Mfoil()
    JFoil.naca_points!(M, "0012")
    JFoil.make_panels!(M, 199, nothing)
    M.oper.alpha = 0.0
    M.param.verb = 0
    JFoil.solve_inviscid!(M)

    # Velocity at a point far from airfoil should be close to freestream
    V = JFoil.inviscid_velocity(M.foil.x, M.isol.gam, 1.0, 0.0, [10.0, 0.0], false)
    @test V[1] ≈ 1.0 atol=0.05
    @test abs(V[2]) < 0.05
end
