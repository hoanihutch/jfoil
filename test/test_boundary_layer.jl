using Test
using JFoil

@testset "Boundary Layer - Phase 5" begin

    @testset "Chunk A: build_param and station_param!" begin
        M = Mfoil()
        init_thermo!(M)

        @testset "build_param - lower surface (si=1)" begin
            p = build_param(M, 1)
            @test p.wake == false
            @test p.turb == false
            @test p.simi == false
            p.turb = true
            @test M.param.turb == false
        end

        @testset "build_param - upper surface (si=2)" begin
            p = build_param(M, 2)
            @test p.wake == false
            @test p.turb == false
            @test p.simi == false
        end

        @testset "build_param - wake (si=3)" begin
            p = build_param(M, 3)
            @test p.wake == true
            @test p.turb == true
            @test p.simi == false
        end

        @testset "build_param preserves thermo params" begin
            M.param.ncrit = 7.5
            M.param.Cuq = 2.0
            p = build_param(M, 1)
            @test p.ncrit == 7.5
            @test p.Cuq == 2.0
        end

        @testset "station_param!" begin
            M.vsol.turb = vcat(falses(10), trues(11))
            M.isol.Istag = [5, 6]

            param = build_param(M, 2)

            station_param!(M, param, 1)
            @test param.turb == false
            @test param.simi == false

            station_param!(M, param, 15)
            @test param.turb == true
            @test param.simi == false

            station_param!(M, param, 5)
            @test param.turb == false
            @test param.simi == true

            station_param!(M, param, 6)
            @test param.turb == false
            @test param.simi == true

            M.vsol.turb[5] = true
            station_param!(M, param, 5)
            @test param.turb == true
            @test param.simi == true
        end
    end

    @testset "Chunk B: stagnation_state and thwaites_init" begin

        @testset "stagnation_state" begin
            # Two states near stagnation: ue proportional to x
            K = 5.0  # ue/x constant
            x = [0.001, 0.002]
            U = zeros(4, 2)
            U[1, :] = [0.001, 0.0012]   # th
            U[2, :] = [0.0022, 0.0026]  # ds
            U[3, :] = [0.0, 0.0]        # sa
            U[4, :] = K .* x            # ue = K*x

            Ust, Ust_U, Ust_x, xst = stagnation_state(U, x)

            # xst should be small
            @test xst ≈ 1e-6

            # ue at stag should be K*xst (very small)
            @test Ust[4] ≈ K * xst atol=1e-12

            # th, ds extrapolated linearly to x=0
            # w1 = x2/dx = 2.0, w2 = -x1/dx = -1.0
            @test Ust[1] ≈ 2.0 * 0.001 + (-1.0) * 0.0012 atol=1e-12  # 0.0008
            @test Ust[2] ≈ 2.0 * 0.0022 + (-1.0) * 0.0026 atol=1e-12  # 0.0018

            # Jacobian sizes
            @test size(Ust_U) == (4, 8)
            @test size(Ust_x) == (4, 2)

            # Finite-difference check on Ust_U
            ε = 1e-7
            for k in 1:8
                Up = copy(U)
                col = (k - 1) ÷ 4 + 1
                row = (k - 1) % 4 + 1
                Up[row, col] += ε
                Ustp, _, _, _ = stagnation_state(Up, x)
                fd = (Ustp - Ust) / ε
                for j in 1:4
                    @test Ust_U[j, k] ≈ fd[j] atol=1e-5
                end
            end

            # Finite-difference check on Ust_x
            for k in 1:2
                xp = copy(x)
                xp[k] += ε
                Ustp, _, _, _ = stagnation_state(U, xp)
                fd = (Ustp - Ust) / ε
                for j in 1:4
                    @test Ust_x[j, k] ≈ fd[j] atol=1e-5
                end
            end
        end

        @testset "thwaites_init" begin
            # Known test: K=5.0, nu=1.5e-5 (air at standard conditions)
            K = 5.0
            nu = 1.5e-5
            th, ds = thwaites_init(K, nu)

            # th = sqrt(0.45*nu/(6*K))
            @test th ≈ sqrt(0.45 * nu / (6.0 * K))
            @test ds ≈ 2.2 * th

            # Sanity: th should be small and positive
            @test th > 0
            @test ds > th  # ds/th = 2.2 > 1
        end
    end

end
