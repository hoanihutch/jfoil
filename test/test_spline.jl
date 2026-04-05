@testset "quadseg" begin
    x, w = JFoil.quadseg()
    @test length(x) == 5
    @test length(w) == 5
    @test x[1] ≈ 0.046910077030668
    @test x[3] ≈ 0.500000000000000
    @test w[3] ≈ 0.284444444444444
    @test sum(w) ≈ 1.0 atol=1e-14
end

@testset "CubicSpline1D" begin
    # Test against scipy CubicSpline (not-a-knot) reference
    x = [0.0, 1.0, 3.0, 5.0, 7.0]
    y = [1.0, 2.0, 0.5, 3.0, 1.5]
    cs = JFoil.CubicSpline1D(x, y)

    # Check interpolation at knots
    for i in 1:5
        @test cs(x[i]) ≈ y[i] atol=1e-12
    end

    # Reference from scipy: cs(0.5) = 1.0, cs(1.5) = 0.5 (approx from scipy)
    # Scipy coefficients (c[0] = cubic term):
    # [[ 0.3377193   0.3377193  -0.26754386 -0.26754386]
    #  [-1.93421053 -0.92105263  1.10526316 -0.5       ]
    #  [ 2.59649123 -0.25877193  0.10964912  1.32017544]
    #  [ 1.          2.          0.5         3.        ]]
    # Our c layout is same: c[1,i]=cubic, c[4,i]=constant
    @test cs.c[1, 1] ≈ 0.3377193 atol=1e-5
    @test cs.c[4, 1] ≈ 1.0 atol=1e-12
    @test cs.c[4, 2] ≈ 2.0 atol=1e-12

    # Derivative test
    dcs = JFoil.derivative(cs)
    @test dcs(0.0) ≈ 2.59649123 atol=1e-5
    @test dcs(1.0) ≈ -0.25877193 atol=1e-5
end

@testset "spline2d" begin
    X = [0.0 0.5 1.0 0.5 0.0;
         0.0 0.1 0.0 -0.1 0.0]
    PP = JFoil.spline2d(X)

    smax = PP.X.x[end]
    @test smax ≈ 2.129754998959201 atol=1e-4

    # Evaluate at endpoints
    xy0 = JFoil.splineval(PP, [0.0])
    @test xy0[1, 1] ≈ 0.0 atol=1e-6
    @test xy0[2, 1] ≈ 0.0 atol=1e-6

    xy_end = JFoil.splineval(PP, [smax])
    @test xy_end[1, 1] ≈ 0.0 atol=1e-4
    @test xy_end[2, 1] ≈ 0.0 atol=1e-4

    # Midpoint
    xy_mid = JFoil.splineval(PP, [smax / 2])
    @test xy_mid[1, 1] ≈ 1.0 atol=1e-3
    @test xy_mid[2, 1] ≈ 0.0 atol=1e-3
end

@testset "splinetan" begin
    X = [0.0 0.5 1.0 0.5 0.0;
         0.0 0.1 0.0 -0.1 0.0]
    PP = JFoil.spline2d(X)
    smax = PP.X.x[end]

    # Tangent at start should point roughly in +x direction with +z component
    xys = JFoil.splinetan(PP, [0.0])
    @test xys[2, 1] > 0  # positive z-tangent at start (going CCW)
end
