@testset "naca_points!" begin
    # NACA 0012
    M = JFoil.Mfoil()
    JFoil.naca_points!(M, "0012")
    @test M.geom.npoint == 201
    @test M.geom.chord ≈ 1.0
    @test M.geom.name == "NACA 0012"
    @test size(M.geom.xpoint) == (2, 201)

    # First 5 x values (from Python reference)
    @test M.geom.xpoint[1, 1] ≈ 1.0 atol=1e-10
    @test M.geom.xpoint[1, 2] ≈ 0.997515 atol=1e-5
    @test M.geom.xpoint[1, 3] ≈ 0.99301379 atol=1e-5

    # First 5 z values
    @test M.geom.xpoint[2, 1] ≈ -0.00126 atol=1e-5
    @test M.geom.xpoint[2, 2] ≈ -0.00160813 atol=1e-5

    # Last 5 values
    @test M.geom.xpoint[1, end] ≈ 1.0 atol=1e-10
    @test M.geom.xpoint[2, end] ≈ 0.00126 atol=1e-5

    # Symmetry: NACA 0012 is symmetric, so z should be antisymmetric
    mid = 101
    @test M.geom.xpoint[1, mid] ≈ 0.0 atol=1e-10  # x at LE is 0
    @test M.geom.xpoint[2, mid] ≈ 0.0 atol=1e-10  # z at LE is 0

    # NACA 2412
    M2 = JFoil.Mfoil()
    JFoil.naca_points!(M2, "2412")
    @test M2.geom.npoint == 201
    @test M2.geom.chord ≈ 1.0
    @test M2.geom.name == "NACA 2412"

    # Reference values from Python
    @test M2.geom.xpoint[1, 1] ≈ 1.0 atol=1e-10
    @test M2.geom.xpoint[2, 1] ≈ -0.00126 atol=1e-5
    @test M2.geom.xpoint[2, end] ≈ 0.00126 atol=1e-5
    @test M2.geom.xpoint[1, end-1] ≈ 0.997515 atol=1e-5
    @test M2.geom.xpoint[2, end-1] ≈ 0.00177346 atol=1e-4
end
