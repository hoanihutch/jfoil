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
            K = 5.0
            x = [0.001, 0.002]
            U = zeros(4, 2)
            U[1, :] = [0.001, 0.0012]
            U[2, :] = [0.0022, 0.0026]
            U[3, :] = [0.0, 0.0]
            U[4, :] = K .* x

            Ust, Ust_U, Ust_x, xst = stagnation_state(U, x)

            @test xst ≈ 1e-6
            @test Ust[4] ≈ K * xst atol=1e-12
            @test Ust[1] ≈ 2.0 * 0.001 + (-1.0) * 0.0012 atol=1e-12
            @test Ust[2] ≈ 2.0 * 0.0022 + (-1.0) * 0.0026 atol=1e-12
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
            K = 5.0
            nu = 1.5e-5
            th, ds = thwaites_init(K, nu)
            @test th ≈ sqrt(0.45 * nu / (6.0 * K))
            @test ds ≈ 2.2 * th
            @test th > 0
            @test ds > th
        end
    end

    @testset "Chunk C: residual_station" begin
        M = Mfoil()
        init_thermo!(M)

        @testset "laminar station" begin
            param = build_param(M, 2)  # upper, laminar
            param.turb = false
            param.simi = false
            x = [0.01, 0.02]
            U = [0.001  0.0015;
                 0.0022 0.003;
                 2.0    3.0;
                 0.5    0.6]
            Aux = [0.0, 0.0]

            R, R_U, R_x = residual_station(param, x, U, Aux)

            # Python reference values
            @test R[1] ≈ 1.10116912 atol=1e-6
            @test R[2] ≈ -0.15975624 atol=1e-6
            @test R[3] ≈ 1.0 atol=1e-6

            @test size(R_U) == (3, 8)
            @test size(R_x) == (3, 2)

            # Check R_U first row against Python
            @test R_U[1, 1] ≈ -1220.11060965 atol=1e-3
            @test R_U[1, 5] ≈ 533.12962021 atol=1e-3

            # Finite-difference Jacobian check on R_U
            ε = 1e-7
            for k in 1:8
                Up = copy(U)
                col = (k - 1) ÷ 4 + 1
                row = (k - 1) % 4 + 1
                Up[row, col] += ε
                Rp, _, _ = residual_station(param, x, Up, Aux)
                fd = (Rp - R) / ε
                for j in 1:3
                    @test R_U[j, k] ≈ fd[j] rtol=1e-4 atol=1e-8
                end
            end

            # Finite-difference check on R_x
            for k in 1:2
                xp = copy(x)
                xp[k] += ε
                Rp, _, _ = residual_station(param, xp, U, Aux)
                fd = (Rp - R) / ε
                for j in 1:3
                    @test R_x[j, k] ≈ fd[j] rtol=1e-4 atol=1e-8
                end
            end
        end

        @testset "turbulent station" begin
            param = build_param(M, 2)
            param.turb = true
            param.simi = false
            x = [0.01, 0.02]
            U = [0.003  0.004;
                 0.006  0.008;
                 0.05   0.06;
                 0.8    0.9]
            Aux = [0.0, 0.0]

            R, R_U, R_x = residual_station(param, x, U, Aux)

            @test R[1] ≈ 0.75521061 atol=1e-6
            @test R[2] ≈ -0.1268235 atol=1e-5
            @test R[3] ≈ -0.01443033 atol=1e-5

            # Finite-difference Jacobian check
            ε = 1e-7
            for k in 1:8
                Up = copy(U)
                col = (k - 1) ÷ 4 + 1
                row = (k - 1) % 4 + 1
                Up[row, col] += ε
                Rp, _, _ = residual_station(param, x, Up, Aux)
                fd = (Rp - R) / ε
                for j in 1:3
                    @test R_U[j, k] ≈ fd[j] rtol=1e-3 atol=1e-7
                end
            end
        end

        @testset "similarity station" begin
            param = build_param(M, 2)
            param.turb = false
            param.simi = true
            x = [0.001, 0.001]
            U = [0.001  0.001;
                 0.0022 0.0022;
                 0.0    0.0;
                 0.005  0.005]
            Aux = [0.0, 0.0]

            R, R_U, R_x = residual_station(param, x, U, Aux)

            @test R[1] ≈ 3.45355628 atol=1e-5
            @test R[2] ≈ -0.97149626 atol=1e-5
            @test R[3] ≈ 0.0 atol=1e-10

            # Check symmetry: R_U columns for U1 and U2 should be equal at similarity
            @test R_U[1, 1] ≈ R_U[1, 5] atol=1e-8
            @test R_U[3, 3] ≈ 1.0 atol=1e-10  # sa coefficient
            @test R_U[3, 7] ≈ 1.0 atol=1e-10
        end
    end

end
