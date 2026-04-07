# Technical Spec: JFoil Performance Optimization

**Date:** 2026-04-08
**Branch:** master
**Functional Spec:** `.claude/projects/performance-optimization-functional-spec.md`
**Status:** Draft

## 1. Approach

Work proceeds bottom-up in 11 phases (0-10). Low-level allocations are eliminated first (closures, panels), then their consumers (residuals, solver), then higher-level features (workspace reuse, annotations, threading, sweeps). Each layer's types must be stable before consumers are converted.

Where multiple strategies exist (SVector vs MVector vs tuples vs workspace), a **comparison spike** benchmarks 2-4 candidates on a single representative function. The winner is applied; spike code is discarded. Use `BenchmarkTools.@benchmark` for function-level measurements; `@time` for full solves.

All `deepcopy(param)` calls are replaced first (Phase 1) because it's the simplest change with the broadest reach. Parametric type annotations (`f(x::R) where R<:Real`) are introduced in each phase as functions are touched, following `/julia-typing` skill conventions. Threading follows `/julia-threading` skill guidelines.

## 2. Key Decisions

| Decision | Rationale |
|----------|-----------|
| Bottom-up phase ordering | Each layer's types must be stable before consumers are converted; avoids cascading rework |
| Comparison spikes before committing | SVector vs MVector vs tuples have different ergonomics; measure before full rewrite |
| `deepcopy` → save/restore of 3 booleans | `Param` has 28 scalar fields that never change during a solve; only `turb`/`wake`/`simi` toggle |
| `BenchmarkTools` for functions, `@time` for solves | BenchmarkTools gives precise stats for micro-benchmarks; `@time` is practical for multi-second solves |
| Parametric types per `/julia-typing` | `f(x::Real)` prevents specialization; `f(x::R) where R<:Real` lets Julia compile efficient code |
| `@inbounds`/`@inline` after allocation elimination | Bounds-checking becomes measurable only after allocation overhead is removed |
| Threading gated by comparison spike | At N=200, thread overhead may exceed benefit; measure first |
| Sweep functions last | They depend on per-point optimizations; cleanest to test independently |
| Convergence failure → NaN + revert | Follows XFoil convention; prevents unconverged state from poisoning subsequent warm-starts |

## 3. Tasks

### Task 0.1: Create baseline benchmark script

**Files:** `benchmark/bench_viscous.jl` (new), `Project.toml` (add BenchmarkTools)
**Depends on:** None

**Description:**
Create `benchmark/bench_viscous.jl` that:
1. Sets up NACA 2412 at Re=1e5, alpha=5 with 199 panels
2. Runs `solve_viscous!` once to compile, then measures with `@time` over 10 runs (median)
3. Runs function-level benchmarks with `@benchmark` on: `get_Hk`, `residual_station`, `build_gamma!`, `panel_info`
4. Runs a 5-point alpha sweep (0, 2, 4, 6, 8) by calling `solve_viscous!` in a loop with `M.oper.initbl = false` after the first point
5. Runs `@code_warntype` on `get_Hk`, `residual_station`, `build_gamma!`, `solve_glob!` and records any type instabilities as comments
6. Prints and records all baseline numbers as comments in the script

No source code changes. This is the reference.

**Task Test:**
```bash
julia --project benchmark/bench_viscous.jl
```
Script runs without error and prints timing/allocation numbers.

**Regression Test:**
```bash
julia --project -e 'using Pkg; Pkg.test()'
```

---

### Task 1.1: Eliminate `deepcopy(param)` in `build_param`

**Files:** `src/boundary_layer.jl`, `src/types.jl`
**Depends on:** Task 0.1

**Description:**
Replace `deepcopy(M.param)` at `boundary_layer.jl:11` with a lightweight copy. Since `Param` contains only 28 scalar fields (Float64, Int, Bool), add a `copy_param` function in `types.jl`:

```julia
function copy_param(p::Param)
    new_p = Param()
    for f in fieldnames(Param)
        setfield!(new_p, f, getfield(p, f))
    end
    return new_p
end
```

Replace `deepcopy(M.param)` with `copy_param(M.param)` in `build_param`.

**Task Test:**
```bash
julia --project -e 'using Pkg; Pkg.test()'
```

**Regression Test:**
```bash
julia --project -e 'using Pkg; Pkg.test()'
```

---

### Task 1.2: Eliminate `deepcopy(param)` in closures and coupling

**Files:** `src/closures_dissipation.jl`, `src/coupling.jl`
**Depends on:** Task 1.1

**Description:**
Replace `deepcopy` with save/restore at all three remaining sites:

1. `get_cDi_lamwake` (`closures_dissipation.jl:102`): Replace `param = deepcopy(paramin)` with saving and restoring `turb`:
   ```julia
   saved_turb = paramin.turb
   paramin.turb = false
   # ... compute ...
   paramin.turb = saved_turb
   ```
   Note: `paramin` is already a `copy_param` from `build_param`, so mutating it is safe within the call.

2. `get_cttr` (`closures_dissipation.jl:249`): Same pattern — save `param.wake`, set to `false`, compute, restore.

3. `wake_init` area in `coupling.jl:268`: Replace `p = deepcopy(param)` with save/restore of `turb` and `wake` on the existing `param`.

After this task, zero `deepcopy` calls remain in the codebase for `Param`.

**Task Test:**
```bash
julia --project -e 'include("test/test_closures.jl")'
julia --project -e 'include("test/test_coupling.jl")'
```

**Regression Test:**
```bash
julia --project -e 'using Pkg; Pkg.test()'
```

---

### Task 2.1: Comparison spike — `panel_info` return type

**Files:** None committed (spike only, results recorded in benchmark script comments)
**Depends on:** Task 0.1

**Description:**
Implement candidate variants of `panel_info` in a scratch file and benchmark each by calling `panel_linvortex_stream` 40,000 times (simulating the N=200 `build_gamma!` loop). Candidates:

1. **SVector{2}**: Return `SVector{2,Float64}` for `t`, `n`; compute `xz` as SVector inline
2. **Scalars**: Return `t1, t2, n1, n2, x, z, d, r1, r2, theta1, theta2` (11 scalars, no vectors)
3. **Tuple of tuples**: Return `((t1,t2), (n1,n2), x, z, d, r1, r2, theta1, theta2)`

Use `@benchmark` for each. Record decision.

**Task Test:**
Decision recorded as comment in `benchmark/bench_viscous.jl`.

**Regression Test:** N/A (no source changes)

---

### Task 2.2: Convert `panel_info`, panel functions, and `TE_info`

**Files:** `src/geometry.jl`, `src/panels.jl`, `src/JFoil.jl` (if adding `using StaticArrays`), `Project.toml` (if adding StaticArrays)
**Depends on:** Task 2.1

**Description:**
Apply the spike winner from Task 2.1. Regardless of winner, eliminate all heap-allocated 2-element vectors:

1. `panel_info` (`geometry.jl:92-114`): Eliminate `t = [xj2 - xj1, ...]`, `n = [-t[2], t[1]]`, `xz = [xi[1] - xj1, ...]` — 3 allocations per call.
2. `TE_info` (`geometry.jl:64-79`): Same pattern for `t1`, `t2`, `t`, `s`, `p`.
3. All 6 panel functions in `panels.jl`: Update to consume new `panel_info` return type. Eliminate `a = [ug1 * t[1] + wg1 * n[1], ...]` at lines 29-30, 88, 145-146.
4. `norm2` in `utils.jl`: Add overload for chosen type if needed.

Update all callers: `build_gamma!`, `inviscid_velocity`, `calc_ue_m!`, `calc_force!`, `build_wake!`.

If SVector wins: add `using StaticArrays` to `JFoil.jl` and add to `Project.toml`.

**Task Test:**
```bash
julia --project -e 'include("test/test_panels.jl")'
julia --project -e 'include("test/test_geometry.jl")'
```
Run `@benchmark` on `build_gamma!` to measure improvement.

**Regression Test:**
```bash
julia --project -e 'using Pkg; Pkg.test()'
```

---

### Task 3.1: Comparison spike — closure Jacobian return type

**Files:** None committed (spike only)
**Depends on:** Task 0.1

**Description:**
Implement candidate variants of `get_Hk` (which calls `get_H` and optionally `get_Mach2`):

1. **SVector{4,Float64}**: Jacobians as immutable SVectors; zero-out becomes `zero(SVector{4,Float64})`; broadcasting works natively
2. **MVector{4,Float64}**: Mutable static vectors; zero-out becomes `fill!(Hk_U, 0.0)`
3. **NTuple{4,Float64}**: Zero-cost like SVector but requires manual arithmetic
4. **Pre-allocated workspace buffer**: Mutable `Vector{Float64}` passed in

Benchmark each calling `get_Hk` 100,000 times with `@benchmark`. Also test the `x = x .* 0.0` zero-out pattern and `vcat(a, b)` for 4+4→8 concatenation with each approach.

**Task Test:**
Decision recorded in benchmark spike results.

**Regression Test:** N/A (no source changes)

---

### Task 3.2: Convert `closures_shape.jl`

**Files:** `src/closures_shape.jl`
**Depends on:** Task 3.1

**Description:**
Convert all 9 functions (`get_H`, `get_Hw`, `get_Mach2`, `get_Hk`, `get_Hs`, `get_Hss`, `get_Ret`, `get_rho`, `get_de`) to return the chosen type from Task 3.1.

Key conversion patterns (example assuming SVector winner — adapt to actual winner):
- `[-H / U[1], 1.0 / U[1], 0.0, 0.0]` → `SVector(-H / U[1], 1.0 / U[1], 0.0, 0.0)` (line 12)
- `Hk_U = Hk_U .* 0.0` → `Hk_U = zero(SVector{4,Float64})` (lines 82, 85)
- `zeros(4)` → `zero(SVector{4,Float64})` (lines 45, 92, 199)
- `copy(Ret_U)` → just reassign (immutable, line 97)
- `Reb_U = Reb_U .* 0.0` → `Reb_U = zero(SVector{4,Float64})` (line 99)
- Broadcasting `.* ` and `.+` → work on SVectors without allocating

Zero-out patterns in this file: lines 82, 85, 99.
`zeros(4)` calls: lines 45, 92, 199.
`copy()` calls: line 97.

Update function signatures to parametric types per `/julia-typing` conventions.

**Task Test:**
```bash
julia --project -e 'include("test/test_closures.jl")'
```

**Regression Test:**
```bash
julia --project -e 'using Pkg; Pkg.test()'
```

---

### Task 3.3: Convert `closures_friction.jl`

**Files:** `src/closures_friction.jl`
**Depends on:** Task 3.2

**Description:**
Convert all 8 functions (`get_cf`, `get_cfxt`, `get_cfutstag`, `get_Us`, `get_upw`, `upwind`, `get_uq`, `get_damp`).

Special attention:
- `get_cfutstag` (line 74): Replace `U = copy(Uin); U[4] = 0.0` with creating a new value of the chosen type with element 4 zeroed. E.g. `U = SVector(Uin[1], Uin[2], Uin[3], 0.0)`.
- `get_upw` (line 132): `Z = zeros(length(Hk1_U1))` → `zero(...)`. `vcat(Z, ...)` at line 136 and `vcat(...)` at line 140 → produce 8-element type from two 4-element types.
- `upwind` (line 158): `vcat((1.0 - upw) .* f1_U1, upw .* f2_U2)` → produces 8-element type from two 4-element inputs. Works naturally with SVectors.
- `get_damp` (lines 197-256): Most complex — ~10 intermediate Jacobian vectors all become chosen type. `zeros(length(U))` at lines 220, 226 → `zero(...)`.
- `get_cfxt` (line 62): `cfxt_U[1] -= cfxt / U[1]` — mutation. With immutable type, reconstruct: `cfxt_U = SVector(cfxt_U[1] - cfxt/U[1], cfxt_U[2], cfxt_U[3], cfxt_U[4])` or use `setindex` equivalent.

Zero-out patterns: lines 28, 103, 106, 114, 117, 143, 174, 177, 182, 203.
`zeros()` calls: lines 11, 132, 220, 226.
`copy()` calls: line 74.
`vcat()` calls: lines 136, 140, 158.
Mutation of elements: lines 62-63.

The `upwind` function signature must accept `upw_U` as either `0` (scalar, symmetric averaging) or 8-element type (from `get_upw`). Handle with dispatch or branch.

**Task Test:**
```bash
julia --project -e 'include("test/test_closures.jl")'
```

**Regression Test:**
```bash
julia --project -e 'using Pkg; Pkg.test()'
```

---

### Task 3.4: Convert `closures_dissipation.jl`

**Files:** `src/closures_dissipation.jl`
**Depends on:** Task 3.2, Task 1.2

**Description:**
Convert all 10 functions (`get_cDi`, `get_cDi_turbwall`, `get_cDi_lam`, `get_cDi_lamwake`, `get_cDi_outer`, `get_cDi_lamstress`, `get_cDixt`, `get_cdutstag`, `get_cteq`, `get_cttr`).

Special attention:
- `get_cDi` (line 12): `cDi_U = zeros(4)` → `zero(...)`. `.+=` accumulation (lines 16, 23, 26) → immutable: `cDi_U = cDi_U + cDi0_U`.
- `get_cDi_turbwall` (line 50): `return 0.0, zeros(4)` → `return 0.0, zero(...)`.
- `get_cDi_outer` (line 126): Same early-return pattern.
- `get_cdutstag` (line 182): Same `copy(Uin)` + `U[4] = 0.0` pattern as `get_cfutstag` — use same approach from Task 3.3.
- `get_cteq` (line 220): `copy(Hk_U)` → unnecessary with immutable.
- `get_cDixt` (line 170): Same element mutation as `get_cfxt`.
- `get_cDi` (line 33): `cDi_U = cDi_U .* 2.0` → `cDi_U = cDi_U * 2.0` (works with SVectors).

Zero-out patterns: lines 218, 223, 228.
`zeros()` calls: lines 12, 50, 126.
`copy()` calls: line 220.
Mutation patterns: lines 16, 23, 26, 170.

**Task Test:**
```bash
julia --project -e 'include("test/test_closures.jl")'
```

**Regression Test:**
```bash
julia --project -e 'using Pkg; Pkg.test()'
```

---

### Task 4.1: Comparison spike — `residual_station` output types

**Files:** None committed (spike only)
**Depends on:** Task 3.3

**Description:**
Now that closures return the chosen static type, evaluate representations for `residual_station` outputs:

1. `R` as `SVector{3}`, `R_U` as `SMatrix{3,8}`, `R_x` as `SMatrix{3,2}`
2. `R` as `SVector{3}`, `R_U` as three `SVector{8}` rows, `R_x` as three `SVector{2}` rows
3. `R` as 3 scalars, `R_U` and `R_x` as the chosen 4/8-element types composed manually

Benchmark `residual_station` in isolation with `@benchmark`. Measure allocation count and time.

**Task Test:**
Decision recorded.

**Regression Test:** N/A (no source changes)

---

### Task 4.2: Convert `residual_station`

**Files:** `src/boundary_layer.jl`
**Depends on:** Task 4.1, Task 3.3, Task 3.4

**Description:**
Convert `residual_station` (`boundary_layer.jl:100-298`) to use static types. This is the largest single conversion.

Key changes:
1. **Input copy** (line 113): Replace `U = copy(Uin)` with constructing modified column vectors. Wake gap subtraction (`U[2,1] -= Aux[1]`) handled by creating new 4-element vectors with adjusted `ds`:
   ```julia
   U1 = SVector(Uin[1,1], Uin[2,1] - Aux[1], Uin[3,1], Uin[4,1])
   U2 = SVector(Uin[1,2], Uin[2,2] - Aux[2], Uin[3,2], Uin[4,2])
   ```
2. **Scalar extraction** (lines 122-124): `th`, `ds`, `sa` become plain scalars.
3. **8-vector Jacobians** (lines 132-134): `thlog_U = SVector(-1.0/th1, 0, 0, 0, 1.0/th2, 0, 0, 0)` — 8-element static type.
4. **vcat calls** (lines 145, 154, 169, etc.): `vcat(H1_U1, H2_U2)` → produces 8-element type from two 4-element types.
5. **upwind calls** (lines 150, 180, etc.): Already return 8-element type from Task 3.3.
6. **Zero-out patterns** (lines 158-161): `thlog_U .* 0.0` → `zero(...)`.
7. **Array literal Jacobians** (lines 204, 215, 239, 248): `[0, 0, -1, 0, 0, 0, 1, 0]` → 8-element static type.
8. **Output assembly** (lines 293-295): `R = [Rmom, Rshape, Rlag]` → `SVector(Rmom, Rshape, Rlag)`. `vcat(Rmom_U', Rshape_U', Rlag_U')` → chosen 3x8 representation from Task 4.1. `vcat(Rmom_x', Rshape_x', Rlag_x')` → chosen 3x2 representation.

`vcat()` calls: lines 145, 154, 169, 263, 294, 295.
Array literals: lines 122-124, 132-135, 175, 204, 215, 239, 248, 264, 275, 279, 293.
Zero-out patterns: lines 158-161.
`copy()`: line 113.

**Task Test:**
```bash
julia --project -e 'include("test/test_boundary_layer.jl")'
```
Run `@benchmark` on `residual_station` to verify zero allocations.

**Regression Test:**
```bash
julia --project -e 'using Pkg; Pkg.test()'
```

---

### Task 4.3: Convert `stagnation_state`

**Files:** `src/boundary_layer.jl`
**Depends on:** Task 4.2

**Description:**
Convert `stagnation_state` (`boundary_layer.jl:31-81`) to static types. Currently returns `(Vector{Float64}(4), Matrix{Float64}(4,8), Matrix{Float64}(4,2), Float64)`.

Key changes:
- `w1_x`, `w2_x`, `dx_x`, `rx_x` are 2-element vectors → 2-element static type (line 44-49)
- `K_U` is an 8-element vector (line 56) → 8-element static type
- `Ust_U = zeros(4, 8)` (line 66) → build as `SMatrix{4,8}` from computed values
- `Ust_x = zeros(4, 2)` (line 74) → build as `SMatrix{4,2}`
- `Ust` is a 4-vector (line 50) → 4-element static type. Note: `Ust[4] = K * xst` at line 61 mutates — with immutable type, reconstruct.

The function currently builds `Ust_U` and `Ust_x` by setting individual elements in loops (lines 67-78). With static matrices, build from column-major data in constructors.

Note: `Ust` is used downstream in `build_glob_sys!` where `Ust[1:3]` is extracted (solver.jl:349). SVector slicing works: `Ust[SVector(1,2,3)]` returns `SVector{3}`.

**Task Test:**
```bash
julia --project -e 'include("test/test_boundary_layer.jl")'
```

**Regression Test:**
```bash
julia --project -e 'using Pkg; Pkg.test()'
```

---

### Task 4.4: Convert `residual_transition`

**Files:** `src/boundary_layer.jl`
**Depends on:** Task 4.2

**Description:**
Convert `residual_transition` (`boundary_layer.jl:302-418`). This function:
1. Runs Newton loop for transition location `xt` (calling `get_damp` twice per iteration)
2. Builds transition-point states `Utl`, `Utt` with 4x4 Jacobians
3. Calls `residual_station` twice (laminar + turbulent sub-intervals)
4. Combines with chain-rule Jacobian multiplication

Key changes:
- `Ut = w1 * U1 + w2 * U2` (line 333) → works on static vectors
- `Z = zeros(4)` (line 320) → `zero(...)`
- `sa = [U[3,1], U[3,2]]` (line 319) → scalars
- `Ut_U1 = w1 * I + (U2 - U1) * xt_U1' / dx` (line 374): `I` is `LinearAlgebra.I` (UniformScaling). With SMatrix: `w1 * SMatrix{4,4}(I) + ...`
- `copy(Ut)`, `copy(Ut_U1)` etc. (lines 378, 383): ~15 copy calls become unnecessary with immutable types — just reassign.
- `Utl_U1[3, :] .= 0.0` (line 380): With SMatrix, reconstruct with zeroed row.
- `hcat(U1, Utl)` (line 400), `hcat(Utt, U2)` (line 404): Build 4x2 matrix for `residual_station` input. With SVectors, `hcat(SVector{4}, SVector{4})` → `SMatrix{4,2}`.
- `[x1, xt]`, `[xt, x2]` (lines 400, 404): 2-element static type.
- `R_U = hcat(R_U1, R_U2)` (line 411): Combine two 3x4 into 3x8.

`vcat()` calls: line 356.
`hcat()` calls: lines 400, 404, 411, 415.
`copy()` calls: lines 378-379, 383-384 (10 total).
`zeros()`: line 320.

**Task Test:**
```bash
julia --project -e 'include("test/test_boundary_layer.jl")'
```

**Regression Test:**
```bash
julia --project -e 'using Pkg; Pkg.test()'
```

---

### Task 5.1: Convert `inviscid_velocity` and `build_wake!`

**Files:** `src/inviscid.jl`
**Depends on:** Task 2.2

**Description:**
Convert remaining inviscid functions to use the static 2-element type from Phase 2:

1. `inviscid_velocity` (`inviscid.jl:150-195`):
   - `V = zeros(2)` → mutable accumulator (`MVector{2}` or two scalars `v1, v2`)
   - `V_G = zeros(2, N)` → keep as heap `Matrix{Float64}` (variable-size, O(N))
   - `Vinf .* [cosd(alpha), sind(alpha)]` → static 2-element type
   - Inner loop accumulation into `V` needs mutable target if using scalars

2. `build_wake!` (`inviscid.jl:203-251`):
   - `xyw = zeros(2, Nw)`, `tw = zeros(2, Nw)` → keep as heap Matrix (variable-size)
   - `xyte`, `n_vec`, `t_vec` → static 2-element type
   - Small temporary 2-vectors eliminated

3. `stagpoint_find!` (`inviscid.jl:110-141`):
   - `M.isol.xstag` computation: use SVector or view for 2-element stag location

4. Update `build_gamma!(M, alpha::Real)` signature to parametric: `build_gamma!(M::Mfoil, alpha::R) where R<:Real`.

**Task Test:**
```bash
julia --project -e 'include("test/test_inviscid.jl")'
```
Run `@benchmark` on `build_gamma!` + `build_wake!`.

**Regression Test:**
```bash
julia --project -e 'using Pkg; Pkg.test()'
```

---

### Task 6.1: Workspace reuse for large matrices in `build_glob_sys!`

**Files:** `src/solver.jl`
**Depends on:** Task 4.2

**Description:**
The `Glob` struct already stores `R`, `R_U`, `R_x`, `R_V` and has a `realloc` flag. The current code at `solver.jl:294-317` already conditionally allocates vs zeros. Optimize:

1. Replace `M.glob.R .= 0.0` with `fill!(M.glob.R, 0.0)` (line 302) — minor, avoids broadcast machinery.
2. Same for `R_U` (line 309) and `R_x` (line 316): `fill!(M.glob.R_U, 0.0)`, `fill!(M.glob.R_x, 0.0)`.
3. Same for `R_V` in `solve_glob!` (line 470).
4. Remove `Aux = zeros(N)` allocation inside the surface loop (`solver.jl:324`). For surfaces 1 and 2, `Aux` values are 0. For surface 3 (wake), use `M.vsol.wgap` directly. Avoid creating a new vector per surface.
5. Remove repeated `[i-1, i]` array allocation in the march loop (`solver.jl:375`). Pass indices directly or use a tuple/SVector.
6. Set `M.glob.realloc = false` after first allocation so subsequent Newton iterations reuse.

Note: `+=` accumulation (lines 367-389) requires full zeroing. The win is avoiding reallocation, not avoiding zeroing.

**Task Test:**
```bash
julia --project -e 'include("test/test_solver.jl")'
```
Run `@time` on `solve_coupled!` to measure allocation reduction.

**Regression Test:**
```bash
julia --project -e 'using Pkg; Pkg.test()'
```

---

### Task 6.2: Update `build_glob_sys!` consumers for static types

**Files:** `src/solver.jl`
**Depends on:** Task 4.2, Task 4.3, Task 4.4

**Description:**
Now that `residual_station` and `residual_transition` return static types, update the assembly code in `build_glob_sys!` (lines 319-407):

1. **Station residual insertion** (lines 385-389): `Ri_U[:, 1:4]` slicing on static 3x8 type → extract first 4 columns. Assignment into dense `M.glob.R_U[Ig, jcols]` works via `setindex!`.
2. **Stagnation system** (lines 336-354): `R1_Ust * Ust_U` is static 3x4 * 4x8 → 3x8. Assignment into dense R_U works.
3. **`hcat(Ust, Ust)`** (line 338): Build 4x2 input for `residual_station` from two static 4-vectors.
4. **Array literal `[i0, i0+1]`** (line 335), `[Is[i0], Is[i0+1]]` (line 343): → `SVector(i0, i0+1)` or tuple.
5. **`zeros(0, 0)`** (line 358): Empty matrix for wake `R1_x` — keep as-is (branch-specific).

Also update `init_boundary_layer!` (lines 10-208):
- `hcat(Ust, Ust)` (line 69) → static 4x2.
- `zeros(2)` for Aux (line 69) → static 2-element type.
- `[th, ds, 0.0, K * xst]` (line 65) → static 4-vector.
- `vcat(A \ b, 0.0)` (lines 74, 126) → construct 4-element type from 3-element solve + scalar.
- `[i-1, i]` (line 96) → tuple or SVector.

**Task Test:**
```bash
julia --project -e 'include("test/test_solver.jl")'
```

**Regression Test:**
```bash
julia --project -e 'using Pkg; Pkg.test()'
```

---

### Task 7.1: Add `@inline` to hot functions

**Files:** `src/closures_shape.jl`, `src/closures_friction.jl`, `src/closures_dissipation.jl`, `src/geometry.jl`, `src/utils.jl`, `src/panels.jl`
**Depends on:** Task 3.4, Task 2.2

**Description:**
Add `@inline` to small, frequently-called functions:

- **All 9 closures_shape functions**: `get_H`, `get_Hw`, `get_Mach2`, `get_Hk`, `get_Hs`, `get_Hss`, `get_Ret`, `get_rho`, `get_de`
- **All 8 closures_friction functions**: `get_cf`, `get_cfxt`, `get_cfutstag`, `get_Us`, `get_upw`, `upwind`, `get_uq`, `get_damp`
- **All 10 closures_dissipation functions**: `get_cDi` through `get_cttr`
- **Geometry**: `panel_info`, `TE_info`
- **Utils**: `norm2`, `dist`
- **Panel functions**: all 6 functions in `panels.jl`

Note: `@inline` is a hint, not a guarantee. Julia usually inlines small functions anyway, but being explicit helps when the optimizer is borderline.

**Task Test:**
```bash
julia --project -e 'using Pkg; Pkg.test()'
```
Run `@benchmark` on `residual_station` and `build_gamma!` to verify no regression and potential improvement.

**Regression Test:**
```bash
julia --project -e 'using Pkg; Pkg.test()'
```

---

### Task 7.2: Add `@inbounds` to verified-safe inner loops

**Files:** `src/inviscid.jl`, `src/solver.jl`, `src/coupling.jl`
**Depends on:** Task 5.1, Task 6.2

**Description:**
Add `@inbounds` to inner loops where index bounds are verified by construction:

1. `build_gamma!` (`inviscid.jl`): The `for i in 1:N, j in 1:N-1` double loop — all indices are within matrix bounds by construction.
2. `build_glob_sys!` (`solver.jl`): The station march loop `for i in (i0+1):N` — indices derived from `Is` which is validated by `identify_surfaces!`.
3. `calc_ue_m!` (`coupling.jl`): The `for i in 1:Nw, j in ...` loops — indices within matrix bounds.
4. `inviscid_velocity` (`inviscid.jl`): Inner panel loop.
5. `init_boundary_layer!` (`solver.jl`): Station march Newton loop.

Do NOT add `@inbounds` to user-facing functions or functions where indices come from external input.

**Task Test:**
```bash
julia --project -e 'using Pkg; Pkg.test()'
```

**Regression Test:**
```bash
julia --project -e 'using Pkg; Pkg.test()'
```

---

### Task 7.3: Audit parametric types and verify with `@code_warntype`

**Files:** All `src/*.jl` files touched in prior phases
**Depends on:** Task 7.1, Task 7.2

**Description:**
Systematic audit of type annotations following `/julia-typing` skill conventions:

1. Convert remaining `::Real` annotations to parametric form: `f(x::R) where R<:Real`
2. Convert `::AbstractVector` to parametric: `f(U::V) where V<:AbstractVector`
3. Known targets:
   - `build_gamma!(M::Mfoil, alpha::Real)` → `build_gamma!(M::Mfoil, alpha::R) where R<:Real`
   - `inviscid_velocity(X, G, Vinf::Real, alpha::Real, ...)` → parametric
   - `thwaites_init(K::Float64, nu::Float64)` → `K::R, nu::R) where R<:Real`
   - All closure functions with `::AbstractVector` and `::Real` args
   - `upwind(upw::Real, ...)`, `get_uq(ds::Real, ...)`, etc.
4. Run `@code_warntype` on: `get_Hk`, `residual_station`, `build_gamma!`, `solve_glob!`, `panel_info`, `upwind`
5. Fix any `Any`-typed variables revealed by `@code_warntype`

**Task Test:**
Add a section to `benchmark/bench_viscous.jl` that runs `@code_warntype` and asserts no `Any` types (or document remaining ones with justification).

**Regression Test:**
```bash
julia --project -e 'using Pkg; Pkg.test()'
```

---

### Task 8.1: Comparison spike — threading `build_gamma!`

**Files:** None committed (spike only)
**Depends on:** Task 5.1

**Description:**
Add `Threads.@threads` to the outer `for i in 1:N` loop in `build_gamma!`. Each row `i` writes to `A[i, :]` and `rhs[i, :]` independently (no data races between rows).

Follow `/julia-threading` skill guidelines for thread safety.

Benchmark with `JULIA_NUM_THREADS=1` vs `JULIA_NUM_THREADS=4` on N=200 using `@benchmark`. If the threaded version is not at least 1.5x faster, threading is not worth the complexity at this problem size.

Also assess `calc_ue_m!` panel loops (`coupling.jl:65-71`, `78-118`).

**Task Test:**
Decision recorded: "threading beneficial at N=200: yes/no" with measured speedup ratios.

**Regression Test:** N/A (no source changes)

---

### Task 8.2: Apply threading (if beneficial)

**Files:** `src/inviscid.jl`, `src/coupling.jl`, optionally `src/solver.jl`
**Depends on:** Task 8.1

**Description:**
If Task 8.1 shows benefit, apply following `/julia-threading` skill guidelines:

1. `build_gamma!` (`inviscid.jl`): `@threads` on outer `for i in 1:N`. TE source/vortex contributions (lines 31-38) write to `A[i, 1]` and `A[i, N]` — per-row, safe.

2. `calc_ue_m!` Cgam loop (`coupling.jl:65-71`): `@threads` on `for i in 1:Nw`. Each row independent.

3. `calc_ue_m!` B matrix loop (`coupling.jl:78-118`): `@threads` on `for i in 1:N`. Each row independent.

4. `init_boundary_layer!` surface loop (`solver.jl:30-206`): Surfaces si=1 (lower) and si=2 (upper) are independent, write to disjoint `Is` indices. Could run as two `@spawn` tasks. Surface si=3 (wake) depends on both → sequential after parallel. Requires splitting the loop. Only attempt if surface march is long enough to justify spawn overhead.

If Task 8.1 shows no benefit: skip this task entirely.

**Task Test:**
```bash
JULIA_NUM_THREADS=4 julia --project -e 'using Pkg; Pkg.test()'
```
Results must match single-threaded. Run benchmark with 1 and 4 threads.

**Regression Test:**
```bash
julia --project -e 'using Pkg; Pkg.test()'
```

---

### Task 9.1: Implement `aseq` sweep function

**Files:** `src/solver.jl`, `src/JFoil.jl` (export)
**Depends on:** Task 6.1

**Description:**
Add `aseq(M::Mfoil, a1::R, a2::R, da::R) where R<:Real`:

```julia
function aseq(M::Mfoil, a1::R, a2::R, da::R) where R<:Real
    # 1. One-time geometry setup (if not already done)
    #    assert M.foil.N > 0
    #    solve_inviscid!(M) to get AIC, gamref
    #    init_thermo!, build_wake!, stagpoint_find!, identify_surfaces!,
    #    set_wake_gap!, calc_ue_m!

    # 2. Build alpha sequence
    alphas = a1:da:a2

    # 3. Save last converged state for revert on failure
    U_last = copy(M.glob.U)
    turb_last = copy(M.vsol.turb)

    # 4. Loop over operating points
    #    For each alpha:
    #    a) Set M.oper.alpha, M.oper.givencl = false
    #    b) Recompute gam: M.isol.gam = gamref[:,1]*cosd(a) + gamref[:,2]*sind(a)
    #    c) Recompute ueinv, rebuild stagnation sensitivities
    #    d) First point: M.oper.initbl = true; subsequent: false
    #    e) try: init_boundary_layer!, stagpoint_move!, solve_coupled!, calc_force!, get_distributions!
    #       On convergence failure: record NaN, revert to U_last/turb_last
    #       On success: save U_last = copy(M.glob.U), turb_last = copy(M.vsol.turb)
    #    f) Accumulate results

    # 5. Return polar as NamedTuple
    return (alpha=..., cl=..., cd=..., cm=..., cdf=..., cdp=...,
            xtr_upper=..., xtr_lower=...)
end
```

Key insight: `build_gamma!` stores `gamref` at 0 and 90 degrees. The AIC factorization is reused — only the linear combination changes per alpha. Extract an `update_alpha!(M, alpha)` helper from `solve_inviscid!`/`rebuild_isol!` that updates `gam`, `sgnue`, `sstag_ue` without re-solving the panel system.

Export `aseq` from `JFoil.jl`.

**Task Test:**
Add test in `test/test_solver.jl`:
```julia
@testset "aseq" begin
    M = Mfoil(); naca_points!(M, "2412"); make_panels!(M, 199)
    M.oper.Re = 1e5; M.param.verb = 0
    polar = aseq(M, 0.0, 4.0, 2.0)  # alpha = 0, 2, 4

    # Verify against individual solve_viscous! calls
    for (i, a) in enumerate([0.0, 2.0, 4.0])
        M2 = Mfoil(); naca_points!(M2, "2412"); make_panels!(M2, 199)
        M2.oper.Re = 1e5; M2.oper.alpha = a; M2.param.verb = 0
        solve_viscous!(M2)
        @test polar.cl[i] ≈ M2.post.cl rtol=1e-3
        @test polar.cd[i] ≈ M2.post.cd rtol=1e-3
    end
end
```
Note: warm-start may converge to slightly different solutions (different Newton paths), so use `rtol=1e-3`.

**Regression Test:**
```bash
julia --project -e 'using Pkg; Pkg.test()'
```

---

### Task 9.2: Implement `cseq` sweep function

**Files:** `src/solver.jl`, `src/JFoil.jl` (export)
**Depends on:** Task 9.1

**Description:**
Add `cseq(M::Mfoil, cl1::R, cl2::R, dcl::R) where R<:Real`:

Same structure as `aseq` but:
- Sets `M.oper.givencl = true`, `M.oper.cltgt = cl_target` for each point
- `solve_coupled!` solves for alpha via `clalpha_residual` (already handles prescribed-cl at `solver.jl:417-441`)
- Alpha is an output read from `M.oper.alpha` after convergence
- Same convergence failure handling: NaN + revert to last converged state

Export `cseq` from `JFoil.jl`.

**Task Test:**
Add test in `test/test_solver.jl`:
```julia
@testset "cseq" begin
    M = Mfoil(); naca_points!(M, "2412"); make_panels!(M, 199)
    M.oper.Re = 1e5; M.param.verb = 0
    polar = cseq(M, 0.2, 0.6, 0.2)  # cl = 0.2, 0.4, 0.6

    for (i, clt) in enumerate([0.2, 0.4, 0.6])
        @test polar.cl[i] ≈ clt rtol=1e-2
    end
    @test all(polar.cd .> 0)
    @test all(polar.cd .< 0.1)
end
```

**Regression Test:**
```bash
julia --project -e 'using Pkg; Pkg.test()'
```

---

### Task 9.3: Update benchmark script with sweep benchmarks

**Files:** `benchmark/bench_viscous.jl`
**Depends on:** Task 9.1, Task 9.2

**Description:**
Update the benchmark script to:
1. Benchmark `aseq(M, 0, 10, 1)` (11-point alpha sweep) with `@time`
2. Compare against 11 individual cold-start `solve_viscous!` calls
3. Report the speedup ratio
4. Benchmark `cseq(M, 0.2, 1.0, 0.1)` similarly
5. Record final numbers alongside Phase 0 baseline

**Task Test:**
```bash
julia --project benchmark/bench_viscous.jl
```

**Regression Test:**
```bash
julia --project -e 'using Pkg; Pkg.test()'
```

---

### Task 10.1: Final benchmark and assessment

**Files:** `benchmark/bench_viscous.jl` (update with final numbers)
**Depends on:** All previous tasks

**Description:**
1. Run the full benchmark suite (single-point + sweeps) and record final numbers
2. Compare against Phase 0 baseline:
   - Single-point: time, allocations, bytes
   - 5-point sweep: time, allocations, bytes
   - 11-point aseq: time vs 11 x cold-start
3. Run `@benchmark` on hot functions to verify zero allocations: `get_Hk`, `residual_station`, `panel_info`, `upwind`
4. Run `@code_warntype` on key functions — confirm no `Any` types
5. Run tests with `JULIA_NUM_THREADS=1` and `JULIA_NUM_THREADS=4`
6. Update baseline comments in benchmark script with before/after numbers
7. Identify any phases that underperformed; document whether further optimization is warranted
8. Clean up comparison-spike leftovers

**Task Test:**
All success criteria from functional spec verified:
1. All tests pass
2. Single-point measurably faster than baseline
3. Zero allocations in hot-path functions
4. 11-point sweep < 3x single cold-start time
5. Thread-safe with `JULIA_NUM_THREADS=4`
6. `@code_warntype` clean on key functions

**Regression Test:**
```bash
julia --project -e 'using Pkg; Pkg.test()'
JULIA_NUM_THREADS=4 julia --project -e 'using Pkg; Pkg.test()'
```

---

## 4. Regression Test Suite

After every task:

```bash
# Full test suite (must always pass)
julia --project -e 'using Pkg; Pkg.test()'

# Quick checks for specific phases:
# Phase 1 (deepcopy): closures + boundary layer + solver + coupling
julia --project -e 'using Test; include("test/test_closures.jl"); include("test/test_boundary_layer.jl"); include("test/test_solver.jl"); include("test/test_coupling.jl")'

# Phase 2 (panels): panels + geometry + inviscid
julia --project -e 'using Test; include("test/test_panels.jl"); include("test/test_geometry.jl"); include("test/test_inviscid.jl")'

# Phase 3 (closures): closures
julia --project -e 'using Test; include("test/test_closures.jl")'

# Phase 4 (residuals): boundary layer + solver
julia --project -e 'using Test; include("test/test_boundary_layer.jl"); include("test/test_solver.jl")'

# Phase 8 (threading): full suite with threads
JULIA_NUM_THREADS=4 julia --project -e 'using Pkg; Pkg.test()'
```

## 5. Final Verification

Once all tasks are complete, verify the functional spec's success criteria:

1. **All tests pass:**
   ```bash
   julia --project -e 'using Pkg; Pkg.test()'
   JULIA_NUM_THREADS=4 julia --project -e 'using Pkg; Pkg.test()'
   ```

2. **Single-point speedup:**
   ```bash
   julia --project benchmark/bench_viscous.jl
   ```
   Compare against Phase 0 baseline. Each phase should show measurable improvement.

3. **Zero hot-path allocations:**
   `@benchmark` on `get_Hk`, `residual_station`, `panel_info`, `upwind` shows 0 allocations.

4. **11-point sweep < 3x single cold-start:**
   `aseq(M, 0, 10, 1)` timing from benchmark script.

5. **Type stability:**
   `@code_warntype` on `get_Hk`, `residual_station`, `build_gamma!`, `solve_glob!`, `panel_info`, `upwind` — no `Any`-typed variables.

6. **Convergence failure handling:**
   `aseq` with an alpha range that includes near-stall points records NaN and continues.
