using Test
using JFoil
using LinearAlgebra

@testset "Phase 8: Post-Processing" begin

    @testset "8.1: get_distributions!" begin
        M = Mfoil()
        naca_points!(M, "2412")
        make_panels!(M, 199)
        M.oper.alpha = 2.0
        M.oper.viscous = true
        M.oper.Re = 1e6
        M.oper.Ma = 0.0
        M.param.verb = 0

        solve_viscous!(M)

        # get_distributions! is called inside solve_viscous!, so post should be filled
        @test length(M.post.th) == M.glob.Nsys
        @test length(M.post.ds) == M.glob.Nsys
        @test length(M.post.sa) == M.glob.Nsys
        @test length(M.post.ue) == M.glob.Nsys
        @test length(M.post.uei) == M.glob.Nsys
        @test length(M.post.cf) == M.glob.Nsys
        @test length(M.post.Ret) == M.glob.Nsys
        @test length(M.post.Hk) == M.glob.Nsys

        # Norm checks against Python reference
        @test norm(M.post.th) ≈ 2.305542897896558e-02 rtol=1e-2
        @test norm(M.post.ds) ≈ 3.921014618069256e-02 rtol=1e-2
        @test norm(M.post.sa) ≈ 4.191455727787810e+01 rtol=1e-2
        @test norm(M.post.ue) ≈ 1.596621631868030e+01 rtol=1e-3
        @test norm(M.post.uei) ≈ 1.594912648139533e+01 rtol=1e-3
        @test norm(M.post.cf) ≈ 6.092328356465126e-02 rtol=1e-2
        @test norm(M.post.Ret) ≈ 2.192218186095311e+04 rtol=1e-2
        @test norm(M.post.Hk) ≈ 3.697652408162989e+01 rtol=1e-2

        # Spot checks on first/last values
        # Python post.th[:5] first = 0.001082614719283 (0-based index 0 = Julia index 1)
        @test M.post.th[1] ≈ 0.001082614719283 rtol=1e-2
        # Python post.cf[-5:] = [0,0,0,0,0] (wake has zero cf)
        @test all(M.post.cf[end-4:end] .≈ 0.0)
        # Hk in wake should be > 1
        @test all(M.post.Hk[end-4:end] .> 1.0)
    end

    @testset "8.2: mgeom_flap!" begin
        M = Mfoil()
        naca_points!(M, "2412")
        make_panels!(M, 199)

        npoint_before = M.geom.npoint
        mgeom_flap!(M, [0.7, 0.0], 10.0)

        @test M.geom.npoint == npoint_before  # should keep same count
        # Check that trailing edge has been deflected downward
        # Python: xpoint[1,0] ~ -0.053 (first point, z-coord)
        @test M.geom.xpoint[2, 1] ≈ -0.053335311068874 rtol=1e-3
        @test M.geom.chord ≈ 1.0 rtol=1e-10  # chord unchanged
    end

    @testset "8.3: mgeom_addcamber!" begin
        M = Mfoil()
        naca_points!(M, "0012")
        make_panels!(M, 199)

        x_before = copy(M.geom.xpoint)
        xcamb = zeros(2, 11)
        xcamb[1, :] = range(0, 1, length=11)
        xcamb[2, :] = 0.02 .* sin.(π .* range(0, 1, length=11))
        mgeom_addcamber!(M, xcamb)

        # Max dz should be close to 0.02
        max_dz = maximum(abs.(M.geom.xpoint[2, :] .- x_before[2, :]))
        @test max_dz ≈ 1.999737888024614e-02 rtol=1e-3
    end

    @testset "8.4: mgeom_derotate!" begin
        M = Mfoil()
        naca_points!(M, "2412")
        make_panels!(M, 199)

        mgeom_derotate!(M)

        # After derotation, TE y should be close to original
        @test M.geom.xpoint[2, 1] ≈ -1.26e-3 atol=1e-4
        # LE should have min x
        iLE = argmin(M.geom.xpoint[1, :])
        @test M.geom.xpoint[1, iLE] < 0.01  # LE near x=0
    end

    @testset "8.5: check_ping" begin
        # Use f(x)=x^3, f'(x)=3x^2 — FD has 2nd-order truncation error
        ep = 0.01; x0 = 1.0
        f(x) = x^3; fp(x) = 3*x^2
        v = [f(x0), f(x0+ep), f(x0+2*ep)]
        v_u = [fp(x0), fp(x0+ep), fp(x0+2*ep)]
        E, rate = check_ping(ep, v, v_u, "test cubic")
        @test rate ≈ 2.0 atol=0.3
    end

    @testset "8.5: ping_test! runs" begin
        M = Mfoil()
        naca_points!(M, "2412")
        make_panels!(M, 199)
        M.param.verb = 0

        # Verify it runs without errors (full derivative check)
        ping_test!(M)
        @test true
    end

end
