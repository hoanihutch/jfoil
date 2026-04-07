using Test
using JFoil
using LinearAlgebra

@testset "Phase 7: Coupled Solver" begin

    # Shared setup: NACA 2412, Ma=0, Re=1e6, alpha=2
    function setup_viscous()
        M = Mfoil()
        naca_points!(M, "2412")
        make_panels!(M, 199)
        init_thermo!(M)
        M.oper.alpha = 2.0
        M.oper.viscous = true
        M.oper.Re = 1e6
        M.oper.Ma = 0.0
        M.param.verb = 0
        build_gamma!(M, 2.0)
        build_wake!(M)
        stagpoint_find!(M)
        identify_surfaces!(M)
        set_wake_gap!(M)
        calc_ue_m!(M)
        return M
    end

    @testset "7.1: init_boundary_layer!" begin
        M = setup_viscous()
        init_boundary_layer!(M)

        @test size(M.glob.U) == (4, 230)
        @test M.glob.Nsys == 230

        # Python ref: norm(U[0,:]) = 6.036250904548351e-02  (row 1 = theta)
        @test norm(M.glob.U[1, :]) ≈ 6.036250904548351e-02 rtol=1e-4
        # norm(U[1,:]) = 1.927191439005563e-01  (row 2 = delta*)
        @test norm(M.glob.U[2, :]) ≈ 1.927191439005563e-01 rtol=1e-4
        # norm(U[2,:]) = 2.982365028343824e+01  (row 3 = sa)
        @test norm(M.glob.U[3, :]) ≈ 2.982365028343824e+01 rtol=1e-3
        # norm(U[3,:]) = 1.608670710255753e+01  (row 4 = ue)
        @test norm(M.glob.U[4, :]) ≈ 1.608670710255753e+01 rtol=1e-4

        # turb sum = 56
        @test sum(M.vsol.turb) == 56

        # Surface 0 (lower), Is[0] in Python = Is[1] in Julia
        # Python Is[0][0] = 94 -> Julia Is[1][1] = 95
        Is_lower = M.vsol.Is[1]
        # U[:,Is[0]] = stag point: th ~ 1.054e-4, ds ~ 2.349e-4, sa ~ 0, ue ~ 0.128
        @test M.glob.U[1, Is_lower[1]] ≈ 1.053819193182027e-04 rtol=1e-3
        @test M.glob.U[2, Is_lower[1]] ≈ 2.349498067840702e-04 rtol=1e-3
        @test M.glob.U[3, Is_lower[1]] ≈ 0.0 atol=1e-10
        @test M.glob.U[4, Is_lower[1]] ≈ 1.277569884358418e-01 rtol=1e-3

        # Surface 1 (upper), last point (TE)
        Is_upper = M.vsol.Is[2]
        @test M.glob.U[1, Is_upper[end]] ≈ 0.012105140096354 rtol=1e-3
        @test M.glob.U[2, Is_upper[end]] ≈ 0.030262850240885 rtol=1e-3

        # Surface 2 (wake), first point
        Is_wake = M.vsol.Is[3]
        @test M.glob.U[1, Is_wake[1]] ≈ 0.014654729364969 rtol=1e-3
        @test M.glob.U[2, Is_wake[1]] ≈ 0.046678038772966 rtol=1e-3

        # Wake last point
        @test M.glob.U[1, Is_wake[end]] ≈ 0.005415986346254 rtol=1e-3
        @test M.glob.U[4, Is_wake[end]] ≈ 0.993192700551308 rtol=1e-3
    end

    @testset "7.2: stagpoint_move!" begin
        M = setup_viscous()
        init_boundary_layer!(M)

        sstag_before = M.isol.sstag
        stagpoint_move!(M)

        # Python: stag point did not move to a new panel
        @test M.isol.Istag == [95, 96]  # Julia 1-based (Python [94, 95])
        # sstag should be very close to before (same panel)
        @test M.isol.sstag ≈ sstag_before rtol=1e-6

        # sstag_ue
        @test M.isol.sstag_ue[1] ≈ 0.000281844100594 rtol=1e-3
        @test M.isol.sstag_ue[2] ≈ -0.012811317208623 rtol=1e-3

        # xstag
        @test M.isol.xstag[1] ≈ 0.000999355638379 rtol=1e-3
        @test M.isol.xstag[2] ≈ -0.005449030581565 rtol=1e-3
    end

    @testset "7.3: build_glob_sys!" begin
        M = setup_viscous()
        init_boundary_layer!(M)
        stagpoint_move!(M)
        M.glob.realloc = true
        build_glob_sys!(M)

        @test length(M.glob.R) == 3 * 230
        @test size(M.glob.R_U) == (3 * 230, 4 * 230)

        # Python: norm(R) = 6.192132146853815e-01
        @test norm(M.glob.R) ≈ 6.192132146853815e-01 rtol=1e-3

        # First 5 residual entries should be near zero (stagnation)
        @test all(abs.(M.glob.R[1:5]) .< 1e-10)
    end

    @testset "7.4: clalpha_residual" begin
        M = setup_viscous()
        init_boundary_layer!(M)
        stagpoint_move!(M)
        M.glob.realloc = true
        build_glob_sys!(M)
        calc_force!(M)

        # cl-constrained mode
        M.oper.givencl = true
        M.oper.cltgt = 0.5
        Rcla, Ru_alpha, Rcla_U = clalpha_residual(M)

        @test Rcla ≈ 2.621158572303939e-02 rtol=1e-3
        @test norm(Ru_alpha) ≈ 1.062502741889674e+00 rtol=1e-3
        @test Rcla_U[end] ≈ M.post.cl_alpha rtol=1e-6
        @test norm(Rcla_U) ≈ 3.312242009129048e-01 rtol=1e-3

        # alpha-prescribed mode
        M.oper.givencl = false
        Rcla2, Ru_alpha2, Rcla_U2 = clalpha_residual(M)
        @test Rcla2 == 0.0
        @test Rcla_U2[end] == 1.0
        @test all(Ru_alpha2 .== 0.0)
    end

    @testset "7.5: solve_glob!" begin
        M = setup_viscous()
        init_boundary_layer!(M)
        stagpoint_move!(M)
        M.glob.realloc = true
        build_glob_sys!(M)
        calc_force!(M)

        # cl-constrained mode
        M.oper.givencl = true
        M.oper.cltgt = 0.5
        solve_glob!(M)

        @test size(M.glob.dU) == (4, 230)
        @test norm(M.glob.dU[1, :]) ≈ 2.263746558673215e-02 rtol=1e-3
        @test norm(M.glob.dU[2, :]) ≈ 2.029657738057825e-01 rtol=1e-3
        @test norm(M.glob.dU[3, :]) ≈ 1.852576869804451e+01 rtol=1e-2
        @test norm(M.glob.dU[4, :]) ≈ 1.025798831498111e+00 rtol=1e-3
        @test M.glob.dalpha ≈ -1.656442071208610e+00 rtol=1e-2

        # alpha-prescribed mode
        M.oper.givencl = false
        M.glob.realloc = true
        solve_glob!(M)
        @test norm(M.glob.dU[1, :]) ≈ 1.981702205778808e-02 rtol=1e-3
        @test norm(M.glob.dU[2, :]) ≈ 1.697776145275515e-01 rtol=1e-3
    end

    @testset "7.6: update_state!" begin
        M = setup_viscous()
        init_boundary_layer!(M)
        stagpoint_move!(M)
        M.glob.realloc = true
        build_glob_sys!(M)
        calc_force!(M)

        # Python script ran solve_glob twice: first cl-constrained, then alpha-prescribed
        # update_state was called with the alpha-prescribed dU
        M.oper.givencl = true
        M.oper.cltgt = 0.5
        solve_glob!(M)  # cl-constrained (7.5 test)
        M.oper.givencl = false
        M.glob.realloc = true
        solve_glob!(M)  # alpha-prescribed (7.5b test) -- this dU is used by update_state

        U_before = copy(M.glob.U)
        alpha_before = M.oper.alpha
        update_state!(M)

        # Python: alpha after = 1.490167627369780 (dalpha=0 since givencl=false)
        @test M.oper.alpha ≈ alpha_before atol=1e-10
        # Python: norm(U_after - U_before) = 2.137164820569818
        @test norm(M.glob.U - U_before) ≈ 2.137164820569818 rtol=1e-2
    end

    @testset "7.7-7.8: solve_viscous! (full, alpha-prescribed)" begin
        M = Mfoil()
        naca_points!(M, "2412")
        make_panels!(M, 199)
        M.oper.alpha = 2.0
        M.oper.viscous = true
        M.oper.Re = 1e6
        M.oper.Ma = 0.0
        M.param.verb = 1

        solve_viscous!(M)

        # Python converged values:
        # cl=0.449351, cd=0.005778, cm=-0.048030
        @test M.glob.conv == true
        @test M.post.cl ≈ 4.493512202037986e-01 rtol=1e-3
        @test M.post.cd ≈ 5.778465143722671e-03 rtol=1e-2
        @test M.post.cm ≈ -4.803034279073116e-02 rtol=1e-2
        @test M.post.cdf ≈ 4.115997520571147e-03 rtol=1e-2
        @test M.post.cdpi ≈ 6.819341035993196e-04 rtol=5e-2

        # State norms at convergence
        @test norm(M.glob.U[1, :]) ≈ 2.305542897896558e-02 rtol=1e-2
        @test norm(M.glob.U[2, :]) ≈ 3.921014618069256e-02 rtol=1e-2
    end

    @testset "7.7-7.8: solve_viscous! (cl-constrained)" begin
        M = Mfoil()
        naca_points!(M, "2412")
        make_panels!(M, 199)
        M.oper.viscous = true
        M.oper.Re = 1e6
        M.oper.Ma = 0.0
        M.oper.givencl = true
        M.oper.cltgt = 0.5
        M.oper.alpha = 2.0
        M.param.verb = 1

        solve_viscous!(M)

        # Python converged:
        # cl=0.5, cd=0.005985, alpha=2.379
        @test M.glob.conv == true
        @test M.post.cl ≈ 0.5 atol=1e-3
        @test M.post.cd ≈ 5.984982982290877e-03 rtol=5e-2
        @test M.oper.alpha ≈ 2.379070166622572 rtol=1e-2
    end

end
