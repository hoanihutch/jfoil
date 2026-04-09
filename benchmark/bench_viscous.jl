# benchmark/bench_viscous.jl -- Phase 0 baseline benchmark
#
# Measures solve_viscous! performance (single-point and sweep) and
# function-level benchmarks for hot-path functions.
#
# Usage:
#   julia --project benchmark/bench_viscous.jl
#
# No source code changes -- this is the reference baseline.
#
# ============================================================================
# Phase 0 Baseline Results (2026-04-08)
# Machine: macOS Darwin 25.4.0
# Julia:   see Project.toml
# Test:    NACA 2412, Re=1e5, alpha=5, 199 panels
# ============================================================================
#
# Single-point solve_viscous!:
#   Median time:    0.087 s
#   Median bytes:   375.5 MB
#   Allocations:    4.63M allocs, 358 MiB (from @time)
#
# Function-level benchmarks:
#   get_Hk:           median  12 ns,   96 bytes, 2 allocs
#   residual_station:  median  1.98 us, 17.8 KiB, 348 allocs
#   build_gamma!:      median  3.18 ms, 17.25 MiB, 403k allocs
#   panel_info:        median  42 ns,   320 bytes, 8 allocs
#
# 5-point alpha sweep (0, 2, 4, 6, 8):
#   Cold-start total:   0.75 s
#   Warm-start total:   0.87 s  (0.86x -- no speedup, solve_viscous! redoes geometry)
#
# @code_warntype:
#   get_Hk:            return type Tuple{Float64, Vector{Float64}} -- clean
#   residual_station:  return type Tuple{Vector{Float64}, Matrix{Float64}, Matrix{Float64}} -- clean
#   build_gamma!:      return type Nothing -- clean
#   solve_glob!:       return type Nothing -- clean
# ============================================================================
#
# Spike 2.1: panel_info return type (2026-04-08)
#   Original (heap Vector): 40.3 ns, 320 B, 8 allocs
#   SVector{2}:             9.9 ns,  0 B,   0 allocs  <-- WINNER
#   Scalars:                9.8 ns,  0 B,   0 allocs
#   Decision: SVector{2} -- zero allocs, preserves vector abstraction
#
# Spike 3.1: closure Jacobian return type (2026-04-08)
#   Original (heap Vector): 10.6 ns, 96 B, 2 allocs
#   SVector{4}:             1.6 ns,  0 B, 0 allocs  <-- WINNER
#   MVector{4}:             6.9 ns, 48 B, 1 allocs
#   NTuple{4}:              1.6 ns,  0 B, 0 allocs
#   vcat SVector 4+4->8:    1.3 ns,  0 B (vs 15.8 ns, 128 B for Vector)
#   Decision: SVector{4} -- zero allocs, native vcat, broadcasting
#
# Spike 4.1: residual_station output types (2026-04-08)
#   Current (after Phase 3):  0.59 us, 3.0 KiB, 60 allocs
#   R_U already SMatrix{3,8} from SVector{8} row transposes
#   Decision: R->SVector{3}, R_U->SMatrix{3,8}(already), R_x->SMatrix{3,2}
#   Remaining allocs from: copy(Uin), array literals, R/R_x heap construction
#
# Spike 8.1: threading build_gamma! (2026-04-08)
#   1 thread:  current 1.39ms, threaded 1.45ms (0.96x - overhead)
#   4 threads: current 1.42ms, threaded 0.70ms (2.04x - beneficial)
#   Decision: Apply @threads to build_gamma! outer loop
#
# ============================================================================
# Final Results (2026-04-08, after all phases)
# ============================================================================
#
# Single-point solve_viscous!:
#   Median time:    59 ms  (was 87 ms — 1.5x faster)
#   Allocations:    554k   (was 4.63M — 8.4x fewer)
#   Memory:         162 MiB (was 358 MiB — 2.2x less)
#
# Function-level benchmarks (all ZERO allocations):
#   get_Hk:           2.1 ns, 0 B  (was 12 ns, 96 B — 5.7x faster)
#   residual_station:  146 ns, 0 B  (was 1980 ns, 17.8 KiB — 13.6x faster)
#   panel_info:        8.4 ns, 0 B  (was 42 ns, 320 B — 5.0x faster)
#   upwind:            1.4 ns, 0 B
#
# aseq 11-point sweep: ~1.1s, 11/11 converged, 1.6x vs cold-starts
# cseq 7-point sweep:  ~1.1s, 5/7 converged
# Threading: 2x speedup on build_gamma! with 4 threads
# Type stability: all key functions clean (@code_warntype)
#
# ============================================================================

using JFoil
using BenchmarkTools
using LinearAlgebra
using InteractiveUtils

# ============================================================================
# Setup: NACA 2412 at Re=1e5, alpha=5, 199 panels
# ============================================================================

function setup_mfoil(; alpha=5.0, Re=1e5, npanels=199, verb=0)
    M = Mfoil()
    naca_points!(M, "2412")
    make_panels!(M, npanels)
    M.oper.alpha = alpha
    M.oper.Re = Re
    M.param.verb = verb
    return M
end

# ============================================================================
# 1. Single-point solve_viscous! benchmark
# ============================================================================

println("=" ^ 70)
println("Phase 0 Baseline Benchmark -- JFoil Performance")
println("=" ^ 70)

# Warm-up / compile run
println("\n--- Compile run (first call) ---")
M_compile = setup_mfoil()
@time solve_viscous!(M_compile)
println("  cl = $(M_compile.post.cl), cd = $(M_compile.post.cd)")

# Timed runs (10 runs, report median)
println("\n--- Single-point solve_viscous! (10 runs) ---")
times = Float64[]
allocs_count = Int[]
allocs_bytes = Int[]
for i in 1:10
    M = setup_mfoil()
    stats = @timed solve_viscous!(M)
    push!(times, stats.time)
    push!(allocs_count, Int(stats.gcstats.allocd > 0 ? -1 : 0))  # @timed doesn't give alloc count directly
    push!(allocs_bytes, stats.bytes)
end
median_time = sort(times)[5]  # approximate median
median_bytes = sort(allocs_bytes)[5]
println("  Median time:  $(round(median_time, digits=4)) s")
println("  Median bytes: $(round(median_bytes / 1e6, digits=2)) MB")

# Also use @time for a clean single-run report
println("\n--- Single @time report ---")
M_timed = setup_mfoil()
@time solve_viscous!(M_timed)

# ============================================================================
# 2. Function-level benchmarks with @benchmark
# ============================================================================

println("\n--- Function-level benchmarks ---")

# Setup a solved Mfoil for extracting test inputs
M_solved = setup_mfoil()
solve_viscous!(M_solved)

# 2a. get_Hk
println("\n  get_Hk:")
let
    U = M_solved.glob.U[:, 30]  # a mid-station state vector
    param = build_param(M_solved, 2)  # upper surface param
    b = @benchmark get_Hk($U, $param)
    display(b)
    println()
end

# 2b. residual_station
println("\n  residual_station:")
let
    param = build_param(M_solved, 2)
    Is = M_solved.vsol.Is[2]  # upper surface indices
    i = Is[5]  # pick a station a few from the start
    x = [M_solved.isol.xi[Is[4]], M_solved.isol.xi[Is[5]]]
    Uin = M_solved.glob.U[:, [Is[4], Is[5]]]
    Aux = zeros(2)
    station_param!(M_solved, param, i)
    b = @benchmark residual_station($param, $x, $Uin, $Aux)
    display(b)
    println()
end

# 2c. build_gamma!
println("\n  build_gamma!:")
let
    M = setup_mfoil()
    b = @benchmark build_gamma!($M, 5.0)
    display(b)
    println()
end

# 2d. panel_info
println("\n  panel_info:")
let
    Xj = M_solved.foil.x[:, 1:2]
    xi = M_solved.foil.x[:, 50]
    b = @benchmark panel_info($Xj, $xi)
    display(b)
    println()
end

# ============================================================================
# 3. Five-point alpha sweep (warm-start)
# ============================================================================

println("\n--- 5-point alpha sweep (0, 2, 4, 6, 8) ---")
let
    alphas = [0.0, 2.0, 4.0, 6.0, 8.0]

    # Cold-start baseline: 5 independent solve_viscous! calls
    println("\n  Cold-start (5 independent solves):")
    cold_times = Float64[]
    for a in alphas
        M = setup_mfoil(alpha=a)
        t = @elapsed solve_viscous!(M)
        push!(cold_times, t)
        println("    alpha=$a: $(round(t, digits=4)) s, cl=$(round(M.post.cl, digits=4))")
    end
    println("  Total cold-start time: $(round(sum(cold_times), digits=4)) s")

    # Warm-start: reuse previous solution
    println("\n  Warm-start sweep:")
    M = setup_mfoil(alpha=alphas[1])
    warm_times = Float64[]

    # First point: full cold solve
    t = @elapsed solve_viscous!(M)
    push!(warm_times, t)
    println("    alpha=$(alphas[1]): $(round(t, digits=4)) s (cold), cl=$(round(M.post.cl, digits=4))")

    # Subsequent points: warm start
    for a in alphas[2:end]
        M.oper.alpha = a
        M.oper.initbl = false
        t = @elapsed solve_viscous!(M)
        push!(warm_times, t)
        println("    alpha=$a: $(round(t, digits=4)) s (warm), cl=$(round(M.post.cl, digits=4))")
    end
    println("  Total warm-start time: $(round(sum(warm_times), digits=4)) s")
    println("  Speedup: $(round(sum(cold_times) / sum(warm_times), digits=2))x")
end

# ============================================================================
# 3b. aseq 11-point sweep benchmark
# ============================================================================

println("\n--- aseq 11-point alpha sweep (0:1:10) ---")
let
    # aseq sweep (fresh M each time to avoid stale state)
    M = setup_mfoil(alpha=0.0)
    t_aseq = @elapsed polar = aseq(M, 0.0, 10.0, 1.0)
    nconv = sum(.!isnan.(polar.cl))
    println("  aseq time: $(round(t_aseq, digits=4)) s, $(nconv)/11 converged")

    # Cold-start comparison: 11 individual solve_viscous! calls
    cold_total = 0.0
    for a in 0.0:1.0:10.0
        Mc = setup_mfoil(alpha=a)
        cold_total += @elapsed solve_viscous!(Mc)
    end
    println("  11x cold-start time: $(round(cold_total, digits=4)) s")
    println("  Speedup: $(round(cold_total / t_aseq, digits=2))x")
    println("  Target (< 3x single cold): $(round(t_aseq, digits=4)) < $(round(3 * cold_total / 11, digits=4)) ? $(t_aseq < 3 * cold_total / 11)")
end

# ============================================================================
# 3c. cseq sweep benchmark
# ============================================================================

println("\n--- cseq sweep (0.2:0.1:0.8) ---")
let
    M = setup_mfoil(alpha=0.0)
    t_cseq = @elapsed polar = cseq(M, 0.2, 0.8, 0.1)
    nconv = sum(.!isnan.(polar.cl))
    println("  cseq time: $(round(t_cseq, digits=4)) s, $(nconv)/7 converged")
end

# ============================================================================
# 4. @code_warntype checks on key functions
# ============================================================================

println("\n--- @code_warntype analysis ---")

# get_Hk
println("\n  get_Hk @code_warntype:")
let
    U = M_solved.glob.U[:, 30]
    param = build_param(M_solved, 2)
    println("  (Check output for Any-typed variables)")
    @code_warntype get_Hk(U, param)
end

# residual_station
println("\n  residual_station @code_warntype:")
let
    param = build_param(M_solved, 2)
    Is = M_solved.vsol.Is[2]
    i = Is[5]
    x = [M_solved.isol.xi[Is[4]], M_solved.isol.xi[Is[5]]]
    Uin = M_solved.glob.U[:, [Is[4], Is[5]]]
    Aux = zeros(2)
    station_param!(M_solved, param, i)
    println("  (Check output for Any-typed variables)")
    @code_warntype residual_station(param, x, Uin, Aux)
end

# build_gamma!
println("\n  build_gamma! @code_warntype:")
let
    M = setup_mfoil()
    println("  (Check output for Any-typed variables)")
    @code_warntype build_gamma!(M, 5.0)
end

# solve_glob! -- needs a fully set up system
println("\n  solve_glob! @code_warntype:")
let
    M = setup_mfoil()
    solve_viscous!(M)
    # Rebuild the system so R, R_U etc. are populated
    build_glob_sys!(M)
    println("  (Check output for Any-typed variables)")
    @code_warntype solve_glob!(M)
end

println("\n" * "=" ^ 70)
println("Benchmark complete.")
println("=" ^ 70)
