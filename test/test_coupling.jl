using Test
using JFoil

@testset "Phase 6: Coupling" begin

    # Build a NACA 2412 with wake for all coupling tests
    M = Mfoil()
    naca_points!(M, "2412")
    make_panels!(M, 199)
    init_thermo!(M)
    M.isol.sgnue = ones(M.foil.N)
    build_gamma!(M, 2.0)
    M.oper.alpha = 2.0
    M.oper.viscous = true
    M.oper.Re = 1e6
    M.oper.Ma = 0.0
    build_wake!(M)
    stagpoint_find!(M)

    N = M.foil.N
    Nw = M.wake.N

    @testset "identify_surfaces!" begin
        identify_surfaces!(M)
        # Python (0-based): Is[0]=[94..0], Is[1]=[95..199], Is[2]=[200..229]
        # Julia (1-based):  Is[1]=[95..1], Is[2]=[96..200], Is[3]=[201..230]
        @test length(M.vsol.Is) == 3
        @test length(M.vsol.Is[1]) == 95   # lower
        @test M.vsol.Is[1][1] == 95        # first = Istag[1]
        @test M.vsol.Is[1][end] == 1       # last = 1
        @test length(M.vsol.Is[2]) == 105  # upper
        @test M.vsol.Is[2][1] == 96        # first = Istag[2]
        @test M.vsol.Is[2][end] == 200     # last = N
        @test length(M.vsol.Is[3]) == 30   # wake
        @test M.vsol.Is[3][1] == 201       # first = N+1
        @test M.vsol.Is[3][end] == 230     # last = N+Nw
    end

    @testset "set_wake_gap!" begin
        set_wake_gap!(M)
        @test length(M.vsol.wgap) == Nw
        @test M.vsol.wgap[1] ≈ 0.00251466 atol=1e-6
        @test M.vsol.wgap[2] ≈ 0.00061254 atol=1e-6
        @test M.vsol.wgap[3] ≈ 0.0 atol=1e-10
        @test maximum(M.vsol.wgap) ≈ 2.514663879044001e-3 rtol=1e-8
        @test all(M.vsol.wgap .>= 0.0)
    end

    # Set up BL state for wake_sys / wake_init / calc_ue_m tests
    Ntot = N + Nw
    M.glob.U = zeros(4, Ntot)
    M.vsol.turb = falses(Ntot)

    il = M.vsol.Is[1][end]  # lower TE = 1 (Julia)
    iu = M.vsol.Is[2][end]  # upper TE = 200 (Julia)
    iw = M.vsol.Is[3][1]    # first wake = 201 (Julia)

    M.glob.U[:, il] = [0.005, 0.015, 0.0001, 0.8]
    M.vsol.turb[il] = false
    M.glob.U[:, iu] = [0.004, 0.012, 0.0002, 0.9]
    M.vsol.turb[iu] = false
    M.glob.U[:, iw] = [0.01, 0.03, 0.001, 0.85]
    M.vsol.turb[iw] = true

    @testset "calc_ue_m!" begin
        calc_ue_m!(M)

        @test size(M.vsol.ue_sigma) == (N + Nw, N + Nw - 2)
        @test size(M.vsol.sigma_m) == (N + Nw - 2, N + Nw)
        @test size(M.vsol.ue_m) == (N + Nw, N + Nw)

        # Check ue_sigma values (row 1, first 5 cols)
        @test M.vsol.ue_sigma[1, 1] ≈ -0.37321835 rtol=1e-5
        @test M.vsol.ue_sigma[1, 2] ≈ -0.08022909 rtol=1e-5
        @test M.vsol.ue_sigma[1, 3] ≈ -0.07542016 rtol=1e-5

        # ue_sigma at midpoint
        mid = N ÷ 2 + 1  # Julia 1-based for Python N//2
        @test M.vsol.ue_sigma[mid, 1] ≈ 0.25622239 rtol=1e-5
        @test M.vsol.ue_sigma[mid, 2] ≈ 0.12802499 rtol=1e-5

        # Frobenius norm
        ue_sigma_fro = sqrt(sum(M.vsol.ue_sigma .^ 2))
        @test ue_sigma_fro ≈ 2.033284824011517e+01 rtol=1e-6

        # ue_m norms
        ue_m_fro = sqrt(sum(M.vsol.ue_m .^ 2))
        @test ue_m_fro ≈ 5.732093271986834e+03 rtol=1e-5

        # ue_m spot checks
        @test M.vsol.ue_m[1, 1] ≈ 92.87470892 rtol=1e-4
        @test M.vsol.ue_m[1, 2] ≈ -77.85985329 rtol=1e-4
        @test M.vsol.ue_m[mid, 1] ≈ 63.7604773 rtol=1e-4
    end

    @testset "rebuild_ue_m!" begin
        old_ue_m = copy(M.vsol.ue_m)
        rebuild_ue_m!(M)
        # Should be identical since sgnue hasn't changed
        @test M.vsol.ue_m ≈ old_ue_m atol=1e-12
    end

    @testset "wake_sys" begin
        param = build_param(M, 3)  # Julia side 3 = wake
        R, R_U, J = wake_sys(M, param)

        @test length(R) == 3
        @test R[1] ≈ 0.001 atol=1e-10
        @test R[2] ≈ 0.00048534 rtol=1e-4
        @test R[3] ≈ -0.02828143 rtol=1e-4

        @test size(R_U) == (3, 12)
        @test R_U[1, 1] ≈ -1.0
        @test R_U[1, 5] ≈ -1.0
        @test R_U[1, 9] ≈ 1.0
        @test R_U[2, 2] ≈ -1.0
        @test R_U[2, 6] ≈ -1.0
        @test R_U[2, 10] ≈ 1.0
        @test R_U[3, 11] ≈ 1.0

        # ctw linearization entries
        @test R_U[3, 1] ≈ 9.87189065 rtol=1e-4
        @test R_U[3, 2] ≈ -3.29240094 rtol=1e-4
        @test R_U[3, 5] ≈ 9.86792211 rtol=1e-4
        @test R_U[3, 6] ≈ -3.29104276 rtol=1e-4

        # Python J=[0, 199, 200] -> Julia J=[1, 200, 201]
        @test J == [1, 200, 201]
    end

    @testset "wake_init" begin
        Uw = wake_init(M, 0.85)
        @test length(Uw) == 4
        @test Uw[1] ≈ 0.009 atol=1e-10
        @test Uw[2] ≈ 0.02951466 rtol=1e-4
        @test Uw[3] ≈ 0.02928143 rtol=1e-4
        @test Uw[4] ≈ 0.85
    end

end
