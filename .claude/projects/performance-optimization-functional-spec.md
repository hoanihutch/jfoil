# Functional Spec: JFoil Performance Optimization

**Date:** 2026-04-08
**Branch:** master
**Status:** Draft

## 1. Background

JFoil is a 2D airfoil analysis code translated from Python (`original/mfoil.py`). The translation prioritized correctness -- matching Python output exactly -- over Julia performance idioms. As a result, the code carries pervasive allocation patterns inherited from Python: temporary vectors created on every function call, `deepcopy` of parameter structs, `vcat` for Jacobian assembly, and dense matrix storage for structurally sparse systems.

The solver spends the vast majority of its time in the coupled Newton loop (`solve_coupled!`), which repeatedly calls `build_glob_sys!` -> `calc_force!` -> `solve_glob!` -> `update_state!` -> `stagpoint_move!` -> `update_transition!`. Within `build_glob_sys!`, every station calls `residual_station` (or `residual_transition`), which in turn calls ~15 closure functions, each of which allocates multiple 4-element and 8-element Jacobian vectors.

Beyond single-point performance, XFoil's most common workflow is polar generation: sweeping through a sequence of angles (`ASEQ`) or lift coefficients (`CSEQ`), where each new operating point reuses the previous converged solution as its starting guess. JFoil has the mechanical support for this (`oper.initbl` flag, warm-start path at `solver.jl:19`), but no sweep functions exist, and the per-point overhead makes multi-point runs unnecessarily expensive.

The goal is to reduce allocation pressure and runtime across the solver, add `aseq`/`cseq` sweep functions that exploit warm-start convergence, and preserve all existing test results exactly.

### Relevant code paths

| File | Role | Allocation concern |
|------|------|--------------------|
| `src/closures_shape.jl` | Shape parameter closures (H, Hk, Hs, etc.) | Every function returns `(scalar, Vector{Float64}(4))` -- hundreds of allocations per Newton step |
| `src/closures_friction.jl` | Skin friction, amplification rate | Same pattern; `get_damp` allocates ~10 temporary 4-vectors |
| `src/closures_dissipation.jl` | Dissipation closures | Same pattern; `deepcopy(param)` in `get_cDi_lamwake`, `get_cttr` |
| `src/boundary_layer.jl` | `residual_station`, `residual_transition` | Allocates ~50 intermediate vectors (4-el and 8-el) via `vcat`, array literals, `copy` |
| `src/geometry.jl` | `panel_info` | Returns 2-element vectors `t`, `n` -- called O(N^2) times in `build_gamma!` |
| `src/panels.jl` | Panel influence functions | Allocates 2-element result vectors `[a, b]` on every call |
| `src/inviscid.jl` | `build_gamma!`, `inviscid_velocity`, `build_wake!` | O(N^2) panel calls; `inviscid_velocity` allocates `zeros(2)` per call |
| `src/solver.jl` | `build_glob_sys!`, `solve_glob!`, `init_boundary_layer!` | Allocates/zeros large Jacobian matrices each iteration; warm-start path exists but no sweep driver |
| `src/coupling.jl` | `calc_ue_m!`, `rebuild_ue_m!` | Dense matrix products on structurally sparse matrices |
| `src/types.jl` | All mutable structs | No pre-allocated workspace buffers |

## 2. Current Behaviour

### 2.1 Closure functions (closures_shape, closures_friction, closures_dissipation)

Every closure returns a `(value, linearization)` tuple where the linearization is a freshly allocated `Vector{Float64}` of length 4 (or 8 for combined two-station Jacobians). Typical pattern:

```julia
function get_H(U)
    H = U[2] / U[1]
    H_U = [-H / U[1], 1.0 / U[1], 0.0, 0.0]  # heap allocation
    return H, H_U
end
```

These functions are called hundreds of times per Newton iteration (once per station, per closure, often twice for upwind averaging). Each call allocates 1-3 vectors of length 4.

`get_damp` (`closures_friction.jl:197`) is particularly expensive: it allocates ~10 intermediate 4-vectors (`Hk_U`, `Ret_U`, `Hmi_U`, `aa_U`, `bb_U`, `lrc_U`, `lr_U`, `rn_U`, `rf_U`, `ar_U`, `ex_U`, `da_U`, `af_U`, `damp_U`, `nx_U`, `eex_U`, `ed_U`).

`deepcopy(param)` is used in `get_cDi_lamwake` (`closures_dissipation.jl:102`), `get_cttr` (`closures_dissipation.jl:249`), `build_param` (`boundary_layer.jl:11`), and `wake_init` in `coupling.jl:268` (to flip `turb`/`wake` before calling `get_cttr`). The `Param` struct has 28 fields; deep-copying it is significantly more expensive than flipping a boolean.

### 2.2 Residual computation (boundary_layer.jl)

`residual_station` (`boundary_layer.jl:100-298`) is the single most allocation-heavy function. It:
- Calls `copy(Uin)` to get a mutable 4x2 matrix
- Extracts `th`, `ds`, `sa` as 2-element `Vector{Float64}` from scalar indexing
- Creates `thlog_U`, `uelog_U`, `salog_U` as 8-element vectors from array literals
- Calls ~15 closure functions (each allocating 1-3 vectors)
- Calls `vcat` ~8 times to combine 4-vectors into 8-vectors
- Calls `upwind` ~8 times, each of which calls `vcat` internally
- Builds final `R_U` via `vcat(Rmom_U', Rshape_U', Rlag_U')` (3 transposes + concatenation)

`residual_transition` (`boundary_layer.jl:302-418`) amplifies this: it calls `residual_station` twice (laminar + turbulent sub-intervals), plus performs its own Newton iteration for transition location (calling `get_damp` twice per Newton step), and creates multiple 4x4 matrix copies (`Ut_U1`, `Ut_U2`, `Utl_U1`, etc.).

### 2.3 Panel influence functions (panels.jl, geometry.jl)

`panel_info` (`geometry.jl:92-114`) allocates three 2-element vectors (`t`, `n`, `xz`) on every call. It is called once per panel function, and panel functions are called O(N^2) times in `build_gamma!` and `calc_ue_m!`.

Panel functions (`panel_linvortex_velocity`, `panel_linsource_velocity`, `panel_constsource_velocity`) allocate 2-element result vectors: `a = [ug1 * t[1] + wg1 * n[1], ug1 * t[2] + wg1 * n[2]]`.

`inviscid_velocity` (`inviscid.jl:150`) allocates `V = zeros(2)` and optionally `V_G = zeros(2, N)` on every call. It is called in tight loops within `build_wake!` and `calc_ue_m!`.

### 2.4 Global system assembly and solve (solver.jl)

`build_glob_sys!` (`solver.jl:294-407`) allocates or zeros `R [3*Nsys]`, `R_U [3*Nsys x 4*Nsys]`, and `R_x [3*Nsys x Nsys]` each Newton iteration. For a typical 200-node problem, these are ~600-element, ~480,000-element, and ~120,000-element dense arrays respectively. The code already attempts reuse (`realloc` flag) but still zeros every element each iteration.

`solve_glob!` (`solver.jl:450-505`) allocates `R_V [NN x NN]` where NN ~ 4*Nsys+1. For 200 nodes, this is a ~640,000-element dense matrix. The dense LU solve (`R_V \ R`) is the dominant cost per Newton step, but the matrix has significant structural sparsity (marked "NOTE: sparse candidate" throughout).

### 2.5 Matrix operations in coupling (coupling.jl)

`calc_ue_m!` builds several large matrices (`Cgam [Nw x N]`, `B [N+1 x Npan]`, `Csig [Nw x Npan]`, `sigma_m [Npan x Ntot]`) and performs dense matrix-matrix operations: `Bp = -(M.isol.AIC \ B)` and the final `ue_m = Diagonal(sgue) * ue_sigma * sigma_m`. All are marked as sparse candidates.

### 2.6 No sweep functions (aseq / cseq)

XFoil's `ASEQ` and `CSEQ` commands sweep through sequences of operating points, reusing the previous converged BL solution as the starting guess for the next point. This warm-start is critical: each successive point typically converges in 3-5 Newton iterations instead of 15-30 from cold initialization.

JFoil has the mechanical support for warm-start via `M.oper.initbl`: when `false` and `M.glob.U` has the right dimensions, `init_boundary_layer!` (`solver.jl:19-22`) skips reinitialization and just updates `ue`. However:
- No `aseq` or `cseq` function exists to drive a sweep
- `solve_viscous!` (`solver.jl:679-694`) always calls the full setup sequence (`solve_inviscid!`, `build_wake!`, `calc_ue_m!`, etc.), much of which is redundant between successive points on the same geometry
- There is no polar accumulation (collecting `cl`, `cd`, `cm` vs `alpha` across a sweep)
- The example script (`examples/naca2412.jl`) solves a single operating point

### 2.7 Threading

No multi-threading is used anywhere. Several loops are embarrassingly parallel:
- The O(N^2) panel influence loops in `build_gamma!` (rows are independent)
- The three surface loops (lower, upper, wake) in `build_glob_sys!` and `init_boundary_layer!` (partially independent -- lower/upper are fully independent)
- Station-wise residual computation within each surface (sequential dependency on previous station, so not parallelizable within a surface, but surfaces are independent)

### 2.8 Existing performance annotations

- `@views` is used in 4 places in `boundary_layer.jl` -- good, but inconsistent across the codebase
- No `@inbounds`, `@simd`, or `StaticArrays` usage anywhere
- No pre-allocated workspace buffers in any struct

## 3. Desired Outcome

When this work is complete, JFoil should:

1. **Run the full viscous solve (`solve_viscous!`) as fast as possible**, with zero heap allocations in hot-path functions as the target. Phase 0 will establish a baseline profile to ground expectations.

2. **Eliminate allocations in hot-path functions** through the best-performing approach for each area. Each area will be evaluated via comparison spike (see Section 7) with candidate approaches including `SVector`/`SMatrix`, tuples-of-scalars, pre-allocated workspace buffers, and `@views` into pre-allocated storage. The winner is applied; runners-up are discarded. Target areas:
   - 2-element vectors in `panel_info`, panel functions, and `inviscid_velocity`
   - 4-element Jacobian vectors returned by all closure functions
   - 8-element combined Jacobian vectors in `residual_station` and `upwind`
   - Small fixed-size matrices (4x2, 4x4, 4x8, 3x8, 3x2) in `residual_station`, `residual_transition`, and `stagnation_state`

3. **Eliminate `deepcopy(param)`** in all four call sites (`build_param`, `get_cDi_lamwake`, `get_cttr`, `wake_init` in `coupling.jl`) by passing turbulence/wake/similarity flags as arguments or using save/restore of the toggled fields.

4. **Reuse workspace arrays** for large matrices (`R`, `R_U`, `R_x`, `R_V`) across Newton iterations instead of reallocating. Pre-allocate these in the `Glob` struct and `fill!` to zero each iteration (accumulation via `+=` requires full zeroing; the win is avoiding reallocation, not avoiding zeroing).

5. **Use multi-threading** where independent work exists and benchmarks show net benefit (see `julia-threading` skill for guidelines):
   - Panel influence matrix assembly in `build_gamma!` (row-parallel)
   - Lower/upper surface BL marching in `init_boundary_layer!` (two independent surfaces)
   - Lower/upper surface residual assembly in `build_glob_sys!`

6. **Provide `aseq` and `cseq` sweep functions** that:
   - Accept a range (start, stop, step) of alpha values (`aseq`) or cl values (`cseq`)
   - Reuse the previous converged BL solution as the starting guess for each successive point (warm-start via `initbl = false`)
   - Skip redundant setup work between points on the same geometry (panelling, AIC matrix, coupling matrices stay fixed; only `alpha`/`cltgt` and the inviscid solution change)
   - Accumulate results into a polar table (alpha, cl, cd, cm, cdf, cdp, xtr_upper, xtr_lower)
   - Return the polar as a structured result (NamedTuple of vectors or similar)
   - Handle convergence failure gracefully: skip the failed point (record NaN in polar), revert to the last converged BL state as warm-start for the next point, and continue the sweep

7. **Produce identical numerical results**: all existing tests in `test/runtests.jl` must pass without modification. Floating-point results should be bit-identical (or within existing test tolerances) to the current implementation.

8. **Add `@inbounds` and `@inline`** where appropriate: `@inbounds` on inner loops with verified index safety (panel O(N^2) loops, station marching loops, closure array access); `@inline` on small, frequently-called functions (closures, `panel_info`, `norm2`, `dist`). These become meaningful once allocation overhead is removed and bounds-checking / call overhead become a larger fraction of remaining cost.

9. **Track performance with benchmarks**: use `BenchmarkTools.@btime`/`@benchmark` for function-level micro-benchmarks (closures, `residual_station`, `build_gamma!`); use `@time` for full `solve_viscous!` and sweep benchmarks. A benchmark script reports results and can be used for regression tracking across phases.

10. **Use parametric types for type stability**: function signatures should use parametric type annotations (`f(x::R) where R<:Real` rather than `f(x::Real)`) to allow Julia's specialisation. Verify type stability with `@code_warntype` on key functions after each phase. This applies across the codebase, not just to new code.

## 4. Issues & Challenges

### 4.1 Multiple approaches for eliminating small-vector allocations

- **What:** There are at least four viable approaches for eliminating small temporary vectors in closures and panel functions: (a) `StaticArrays.SVector`/`SMatrix` -- stack-allocated, immutable, zero-cost; (b) tuples-of-scalars (return `t1, t2, n1, n2` instead of `t, n`); (c) pre-allocated mutable workspace buffers passed through the call chain; (d) `@views` into pre-allocated column-major storage. Each has different trade-offs in performance, code complexity, and maintainability.
- **Why it matters:** The best approach may differ by call site and will be determined empirically via comparison spikes (see Section 7). SVectors are natural for 2-element panel geometry; tuples-of-scalars avoid any dependency; workspace buffers avoid type changes but thread through every signature; views are zero-cost but require careful lifetime management. The winner for each area is applied; alternatives are discarded.
- **Evidence:** 2-element vectors in `panel_info` (`geometry.jl:97-99`) could be SVectors or scalar returns. 4-element closure Jacobians are returned ~hundreds of times per Newton step -- SVectors eliminate allocations but require converting ~25 `x = x .* 0.0` patterns. 8-element Jacobians in `residual_station` are built by `vcat` of two 4-vectors and then sliced/transposed.

### 4.2 Closures conditionally zero their Jacobians

- **What:** Many closures contain patterns like `Hk_U = Hk_U .* 0.0` or `Reb_U = Reb_U .* 0.0` to zero out derivatives when a clamp is active. With SVectors, this must become `Hk_U = zero(SVector{4,Float64})` or similar. There are ~25 such patterns across the three closure files.
- **Why it matters:** Each instance must be individually converted. Missing one will cause a method error or type instability.
- **Evidence:** `get_Hs` (`closures_shape.jl:82-86`) has `Hk_U = Hk_U .* 0.0`; `get_damp` (`closures_friction.jl:204`) has `Hk_U = Hk_U .* 0.0`; `get_cteq` (`closures_dissipation.jl:218-229`) has three such patterns.

### 4.3 `vcat` in `upwind` and `residual_station` builds 8-vectors

- **What:** The `upwind` function (`closures_friction.jl:156-160`) uses `vcat((1.0 - upw) .* f1_U1, upw .* f2_U2)` to combine two 4-vectors into an 8-vector. `residual_station` does this ~8 more times. With SVectors, `vcat(SVector{4}, SVector{4})` returns `SVector{8}` -- this works, but the pattern must be consistent.
- **Why it matters:** The 8-vector Jacobians are later sliced, transposed, and assembled into the 3x8 `R_U` matrix. All downstream matrix operations must handle static types.
- **Evidence:** `residual_station` lines 145, 154, 169, 204, 263, 294 all use `vcat` to build 8-vectors.

### 4.4 `deepcopy(param)` is used to toggle 3 boolean flags

- **What:** `build_param` (`boundary_layer.jl:5-16`) does `deepcopy(M.param)` to create a per-side parameter struct, then sets `param.wake`, `param.turb`, `param.simi`. The `Param` struct has 28 scalar fields. `get_cDi_lamwake`, `get_cttr`, and `wake_init` in `coupling.jl` also `deepcopy(param)` to temporarily flip `turb` or `wake`.
- **Why it matters:** `deepcopy` of a 28-field mutable struct is far more expensive than copying 3 booleans. This is called 3 times per Newton iteration in `build_glob_sys!` plus additional times in closure calls and wake initialization.
- **Evidence:** `boundary_layer.jl:11` (`build_param`), `closures_dissipation.jl:102` (`get_cDi_lamwake`), `closures_dissipation.jl:249` (`get_cttr`), `coupling.jl:268` (`wake_init`).

### 4.5 `residual_station` mutates its input copy

- **What:** `residual_station` calls `U = copy(Uin)` then modifies `U[2,1]` and `U[2,2]` (subtracting wake gap). It also extracts column views with `@views U1 = U[:,1]`. The `copy` is necessary to avoid mutating the caller's data, but allocates a 4x2 matrix.
- **Why it matters:** With static arrays, this copy-and-mutate pattern doesn't work. The wake gap subtraction needs to be handled differently (e.g., adjusting `ds` in the closure calls, or using an SMatrix with modified entries).
- **Evidence:** `boundary_layer.jl:113-120`.

### 4.6 Threading in `build_glob_sys!` requires careful index management

- **What:** The three surface loops in `build_glob_sys!` write into different rows/columns of the shared `R`, `R_U`, `R_x` arrays. Lower and upper surfaces write to non-overlapping index ranges, so they can run in parallel. The wake surface also writes to its own range. However, the stagnation equations (lines 346-354) write to indices that overlap with the first stations of both surfaces.
- **Why it matters:** Naive threading will produce data races at the stagnation point indices. The stagnation system must be handled before or after the parallel surface loops.
- **Evidence:** `solver.jl:319-393`. The `Is` arrays for lower/upper/wake surfaces are disjoint (set up by `identify_surfaces!` at `coupling.jl:11-22`), but the stagnation extrapolation at lines 336-354 references indices from both sides.

### 4.7 Panel assembly loop has O(N^2) panel_info calls

- **What:** `build_gamma!` loops `i in 1:N, j in 1:N-1`, calling `panel_linvortex_stream` (which calls `panel_info`) at each iteration. `panel_info` allocates 3 heap vectors (`t`, `n`, `xz`). For N=200, this is ~40,000 calls each allocating 3 vectors = ~120,000 allocations.
- **Why it matters:** This is a one-time cost per solve (not repeated in Newton iterations), but it's still significant and easily reduced with SVectors.
- **Evidence:** `inviscid.jl:18-24` (double loop), `geometry.jl:92-114` (`panel_info`).

### 4.8 Dense linear solve in `solve_glob!` dominates Newton step cost

- **What:** `solve_glob!` performs `dV = -(R_V \ R)` where `R_V` is a dense `[4*Nsys+1 x 4*Nsys+1]` matrix. For 200 nodes, this is ~801x801. The matrix has block-banded structure (3x4 blocks along the diagonal from BL equations, plus dense rows/columns from the ue coupling equation).
- **Why it matters:** The dense LU factorization is O(N^3) and likely dominates the per-iteration cost for larger panels. Exploiting the sparsity structure could yield order-of-magnitude improvements for larger problems. However, the project rules say "no sparse arrays" for CUDA compatibility.
- **Evidence:** `solver.jl:498` (`dV = -(M.glob.R_V \ R)`). The "NOTE: sparse candidate" comments appear at `solver.jl:307,314,468`.

### 4.9 `copy` and `vcat` in `get_cfutstag` and `get_cdutstag`

- **What:** `get_cfutstag` (`closures_friction.jl:73-89`) and `get_cdutstag` (`closures_dissipation.jl:181-199`) both `copy(Uin)` and set `U[4] = 0.0` before calling `get_Hk`. This is a 4-element vector copy just to zero one element.
- **Why it matters:** Minor individually, but these are called during BL initialization which runs a Newton loop per station.
- **Evidence:** `closures_friction.jl:74-75`, `closures_dissipation.jl:182-183`.

### 4.10 Broadcasting on scalars creates unnecessary temporaries

- **What:** Patterns like `num_Hk / Ret .* Hk_U .- num / Ret^2 .* Ret_U` create intermediate vectors at each `.* ` and `.-` operation (2 temporaries per expression). With SVectors, fused broadcasting eliminates these. With heap Vectors, each broadcast allocates.
- **Why it matters:** This pattern appears in virtually every closure function. Eliminating these intermediate allocations requires either SVectors (preferred) or manual loop fusion with `@.`.
- **Evidence:** Present throughout all three `closures_*.jl` files, `get_upw`, `get_uq`, etc.

### 4.11 Sweep functions must separate geometry-fixed from alpha-varying work

- **What:** `solve_viscous!` (`solver.jl:679-694`) calls a sequence of functions, some of which depend only on geometry (panelling, AIC matrix build, wake construction, coupling matrices) and some on operating conditions (`build_gamma!` for the given alpha, `init_boundary_layer!`, `solve_coupled!`). A sweep function must call the geometry-dependent work once, then loop over operating points calling only the alpha/cl-dependent work.
- **Why it matters:** `build_gamma!` does O(N^2) work and a dense LU solve; `calc_ue_m!` does O(N^2) work and multiple dense matrix products. Repeating these for each point in a 20-point polar would waste most of the time on work that doesn't change.
- **Evidence:** `solve_viscous!` calls `solve_inviscid!` (which calls `build_gamma!`), `build_wake!`, `calc_ue_m!`, then `init_boundary_layer!` + `solve_coupled!`. The AIC matrix (`M.isol.AIC`) and coupling matrix (`M.vsol.ue_m`) are geometry-dependent only. `build_gamma!` stores `gamref` at 0 and 90 degrees; the actual `gam` for a given alpha is just `gamref[:,1]*cos(a) + gamref[:,2]*sin(a)` (already done at `inviscid.jl:61`).
- **How `cseq` differs:** For `cseq`, `M.oper.givencl = true` and `M.oper.cltgt` is varied. The coupled solver already handles prescribed-cl via `clalpha_residual` (`solver.jl:417-441`), solving for alpha as part of the Newton system. The sweep function sets `cltgt` and calls the coupled solver.

## 5. Scope

### In scope

- Reducing small-vector allocations in hot-path functions (closures, residuals, panels) -- approach determined per area by comparative benchmarks
- Eliminating `deepcopy(param)` in all four call sites (`build_param`, `get_cDi_lamwake`, `get_cttr`, `wake_init` in `coupling.jl`)
- Pre-allocating workspace arrays in `Glob` struct for Newton iteration reuse
- Adding `@threads` or `@spawn` for independent computations (panel assembly, surface BL loops) where benchmarks show net benefit
- Adding `@inbounds` on verified-safe inner loops and `@inline` on small frequently-called functions
- Using parametric type annotations (`f(x::R) where R<:Real`) and verifying type stability with `@code_warntype`
- `aseq(M, a1, a2, da)` and `cseq(M, cl1, cl2, dcl)` sweep functions with warm-start, polar accumulation, and convergence failure handling
- A benchmark script covering single-point and sweep performance for regression tracking
- Ensuring all existing tests pass unchanged

### Out of scope

- Sparse matrix storage (explicitly excluded per project rules for CUDA compatibility)
- GPU/CUDA porting
- Algorithmic changes (different Newton strategy, different closure formulations)
- Changes to the test suite itself (tests are the correctness oracle)
- Profiling infrastructure beyond a simple benchmark script
- The streamline BL solver (separate functional spec)

## 6. Success Criteria

1. **All existing tests pass**: `julia --project test/runtests.jl` produces zero failures with no test modifications.

2. **Single-point speedup**: `solve_viscous!` on NACA 2412 at Re=1e5, alpha=5 runs as fast as possible. Phase 0 baseline profile will establish where time is spent; subsequent phases each demonstrate measurable improvement against the baseline.

3. **Zero hot-path allocations**: The target is zero heap allocations in all hot-path functions (closures, `residual_station`, `residual_transition`, `panel_info`, panel functions, `inviscid_velocity`). Measured per-function with `BenchmarkTools.@benchmark`.

4. **Sweep functions exist and exploit warm-start**: `aseq(M, 0, 10, 1)` runs a viscous solve at each alpha in 0:1:10, returning a polar table. The total wall-clock time for the 11-point sweep is less than 3x the time for a single cold-start `solve_viscous!` call (demonstrating that warm-start amortises setup cost). Failed points are skipped (NaN in polar), with the last converged state used for the next warm-start.

5. **Benchmark script exists**: A script at `benchmark/bench_viscous.jl` that uses `BenchmarkTools.@benchmark` for function-level micro-benchmarks and `@time` for full `solve_viscous!` / sweep benchmarks, runnable with `julia --project benchmark/bench_viscous.jl`.

6. **Threading correctness**: Running the test suite with `JULIA_NUM_THREADS=4 julia --project test/runtests.jl` produces the same results as single-threaded execution.

7. **Type stability**: Key functions pass `@code_warntype` with no `Any`-typed variables. Function signatures use parametric types (`f(x::R) where R<:Real`) rather than abstract type annotations.

8. **No new dependencies beyond StaticArrays.jl and BenchmarkTools.jl**: `StaticArrays` (if chosen by spikes), `BenchmarkTools` (benchmarks only). Threading uses stdlib `Base.Threads`.

## 7. Phased Schedule

Changes are organised into independent phases so that each can be benchmarked, tested, and assessed in isolation before moving on. Each phase is a separate commit (or branch) with its own benchmark run recorded in the benchmark script output.

Where multiple implementation approaches exist for a phase, a short comparison spike is done first: implement the two or three candidate approaches on a single representative function, benchmark all candidates, then apply the winner to the rest of the phase. The spike code is discarded; only the chosen approach is committed.

### Phase 0 -- Baseline benchmark
- Create `benchmark/bench_viscous.jl` measuring `solve_viscous!` (single point) and a manual 5-point alpha loop (to establish sweep baseline)
- Use `BenchmarkTools.@benchmark` for function-level measurements; `@time` for full solve
- Record: median time, allocation count, allocation bytes for both
- Run `@code_warntype` on key functions (`get_Hk`, `residual_station`, `build_gamma!`, `solve_glob!`) and record any type instabilities
- No code changes -- this is the reference

### Phase 1 -- `deepcopy(param)` elimination
- Replace `deepcopy(M.param)` in `build_param` with a shallow copy or direct field assignment
- Replace `deepcopy(param)` in `get_cDi_lamwake`, `get_cttr`, and `wake_init` (`coupling.jl:268`) with save/restore of the toggled flags
- **Benchmark and test** independently

### Phase 2 -- Panel and geometry small vectors (2-element)
- **Comparison spike:** Convert `panel_info` to return SVectors vs tuples-of-scalars (i.e. return `t1, t2, n1, n2` instead of `t, n`). Benchmark both on `build_gamma!` with `@benchmark`.
- Apply the winner to `panel_info`, all panel functions in `panels.jl`, and `TE_info` in `geometry.jl`
- **Benchmark and test** independently

### Phase 3 -- Closure Jacobians (4-element vectors)
- **Comparison spike:** Convert `get_H` and `get_Hk` to return `SVector{4}` vs `MVector{4}` vs tuples vs pre-allocated workspace. Benchmark with `@benchmark` a standalone loop calling `get_Hk` 10,000 times.
- Apply the winner to all functions in `closures_shape.jl`, `closures_friction.jl`, `closures_dissipation.jl`
- Convert the ~25 `x = x .* 0.0` zero-out patterns to the appropriate idiom
- Convert `copy(Uin)` + mutate in `get_cfutstag` / `get_cdutstag`
- Update function signatures to parametric types where applicable (`f(U::V) where V<:AbstractVector`)
- **Benchmark and test** independently

### Phase 4 -- Residual and upwind 8-vectors
- Convert `upwind`, `get_upw`, `get_uq` to work with the chosen type from Phase 3
- Convert `residual_station` intermediate 8-vectors (`thlog_U`, `uelog_U`, `H_U`, etc.) and the 3x8 / 3x2 output matrices
- Convert `stagnation_state` small matrices (4x8, 4x2)
- Handle `residual_station`'s copy-and-mutate of `Uin` (wake gap subtraction) compatibly with the chosen approach
- **Comparison spike** (approach-dependent): Evaluate candidate representations for `R_U` output (e.g. `SMatrix{3,8}` vs three row-vectors vs mutable workspace). Benchmark `residual_station` in isolation with `@benchmark`.
- Convert `residual_transition` similarly
- **Benchmark and test** independently

### Phase 5 -- Inviscid solver small vectors
- Convert `inviscid_velocity` (`zeros(2)` and `zeros(2, N)`) to use the approach from Phase 2
- Convert `build_wake!` temporary arrays (`xyw`, `tw`, `uewi`)
- Update function signatures to parametric types (`build_gamma!(M, alpha::R) where R<:Real`)
- **Benchmark and test** independently (measure `build_gamma!` + `build_wake!`)

### Phase 6 -- Workspace reuse for large matrices
- Pre-allocate `R`, `R_U`, `R_x`, `R_V` in the `Glob` struct; use `fill!` to zero each iteration (accumulation via `+=` requires full zeroing; the win is avoiding reallocation)
- Avoid re-creating `Aux = zeros(N)` on every `build_glob_sys!` call
- **Benchmark and test** independently (measure full `solve_coupled!` loop)

### Phase 7 -- `@inbounds`, `@inline`, and type stability
- Add `@inbounds` to verified-safe inner loops: panel O(N^2) loops in `build_gamma!`/`calc_ue_m!`, station marching in `build_glob_sys!`, closure array indexing
- Add `@inline` to small frequently-called functions: all closures, `panel_info`, `norm2`, `dist`, `upwind`
- Audit parametric type annotations across all modified functions; verify with `@code_warntype`
- **Benchmark and test** independently

### Phase 8 -- Threading
- **Comparison spike:** Add `@threads` to the row loop in `build_gamma!`. Benchmark single-thread vs 2/4 threads on N=200. Assess whether the overhead is worth it at this problem size. Follow `julia-threading` skill guidelines.
- If beneficial: apply to `build_gamma!` and `calc_ue_m!` panel loops
- Evaluate threading lower/upper surface loops in `init_boundary_layer!` and `build_glob_sys!` (requires stagnation equations handled separately)
- **Benchmark and test** independently with `JULIA_NUM_THREADS=1` and `JULIA_NUM_THREADS=4`

### Phase 9 -- `aseq` / `cseq` sweep functions
- Implement `aseq(M, a1, a2, da)` that:
  - Runs geometry setup once (panelling, AIC, wake, coupling matrices)
  - Loops over alpha values, calling only the alpha-dependent work (update `gamref` combination, `init_boundary_layer!` with `initbl=false`, `solve_coupled!`, `calc_force!`, `get_distributions!`)
  - Accumulates polar results
  - On convergence failure: records NaN for that point, reverts to the last converged BL state, and continues the sweep
- Implement `cseq(M, cl1, cl2, dcl)` similarly with `givencl=true`
- Add tests validating that sweep results match individual `solve_viscous!` calls at each alpha/cl
- **Benchmark:** compare 11-point sweep time vs 11 x single-point cold starts

### Phase 10 -- Final assessment
- Run the full benchmark suite (single-point + sweep) and record final numbers
- Compare against Phase 0 baseline
- Verify `@code_warntype` clean on all key functions
- Identify any regressions introduced by earlier phases
- Clean up any comparison-spike leftovers
- If performance targets are not met, document which phase underperformed and whether further optimization is warranted
