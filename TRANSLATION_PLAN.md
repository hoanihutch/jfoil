# JFoil Translation Plan

Translate mfoil (Python) to Julia function-by-function, ordered by dependency.

Primary reference: `original/mfoil.py`. Cross-reference: `original/mfoil.m`.

## Workflow (every step)

```
1. Build test   — run function in Python, capture reference output
2. Julia test   — write @test assertions from reference data
3. Implement    — translate the function to Julia
4. Pass test    — run tests, iterate until green
5. Commit       — git commit with the function name
```

---

## Phase 0: Data Structures & Utilities

Foundation that everything else depends on.

| Step | Function | Description | Status |
|------|----------|-------------|--------|
| 0.1 | `Geom, Panel, Oper, Isol, Vsol, Glob, Post, Param, Mfoil` | All mutable structs + default constructors | [ ] |
| 0.2 | `vprint` | Conditional printing | [ ] |
| 0.3 | `norm2` | 2-norm of 2-vector | [ ] |
| 0.4 | `dist` | Distance between two points | [ ] |

**Notes:**
- Julia has built-in `sind`, `cosd`, `atan(y,x)` — no translation needed
- Structs go in `src/types.jl`, utils in `src/utils.jl`

---

## Phase 1: Geometry & Paneling

| Step | Function | Deps | Description | Status |
|------|----------|------|-------------|--------|
| 1.1 | `naca_points` | — | Generate NACA 4/5-digit coordinates | [ ] |
| 1.2 | `set_coords` | — | Set geometry from coordinate array | [ ] |
| 1.3 | `quadseg` | — | Quadrature points/weights (constants) | [ ] |
| 1.4 | `spline2d` | `norm2`, `quadseg` | 2D parametric spline with arclength | [ ] |
| 1.5 | `splineval` | `spline2d` | Evaluate spline at parameter values | [ ] |
| 1.6 | `splinetan` | `spline2d` | Spline tangent at parameter values | [ ] |
| 1.7 | `spline_curvature` | `spline2d`, `splineval`, `splinetan`, `quadseg` | Curvature-adaptive point distribution | [ ] |
| 1.8 | `space_geom` | — | Geometric spacing solver | [ ] |
| 1.9 | `TE_info` | `norm2` | Trailing edge geometry | [ ] |
| 1.10 | `panel_info` | `norm2`, `dist` | Panel geometry (angles, distances) | [ ] |
| 1.11 | `clear_solution` | structs | Reset all solution fields | [ ] |
| 1.12 | `make_panels` | `spline_curvature`, `clear_solution` | Create panel mesh from geometry | [ ] |

**Notes:**
- Python uses `scipy.interpolate.CubicSpline` — Julia equivalent: `Interpolations.jl` or manual cubic spline
- `spline2d` returns a dict with spline objects; in Julia use a struct
- Files: `src/naca.jl`, `src/spline.jl`, `src/geometry.jl`

---

## Phase 2: Panel Influence Functions

All pure functions, depend only on `panel_info`.

| Step | Function | Description | Status |
|------|----------|-------------|--------|
| 2.1 | `panel_linvortex_stream` | Linear vortex panel streamfunction | [ ] |
| 2.2 | `panel_linvortex_velocity` | Linear vortex panel velocity | [ ] |
| 2.3 | `panel_constsource_stream` | Constant source panel streamfunction | [ ] |
| 2.4 | `panel_constsource_velocity` | Constant source panel velocity | [ ] |
| 2.5 | `panel_linsource_stream` | Linear source panel streamfunction | [ ] |
| 2.6 | `panel_linsource_velocity` | Linear source panel velocity | [ ] |

**Notes:**
- These implement the analytical formulas from the PDF Appendix A (Eqs. 37-44)
- Each returns influence coefficients (scalars or small arrays)
- File: `src/panels.jl`

---

## Phase 3: Inviscid Solver

| Step | Function | Deps | Description | Status |
|------|----------|------|-------------|--------|
| 3.1 | `init_thermo` | `Param` | Set thermodynamic params from Ma | [ ] |
| 3.2 | `get_cp` | — | Pressure coefficient + linearization | [ ] |
| 3.3 | `get_uk` | — | Karman-Tsien speed + linearization | [ ] |
| 3.4 | `build_gamma` | `TE_info`, `panel_linvortex_stream`, `panel_constsource_stream` | Build AIC, solve for vortex strengths | [ ] |
| 3.5 | `get_ueinv` | `Isol` | Inviscid edge velocity at current alpha | [ ] |
| 3.6 | `get_ueinvref` | `Isol` | Reference ue at 0 and 90 deg | [ ] |
| 3.7 | `stagpoint_find` | `Isol`, `Panel` | Locate stagnation, compute xi | [ ] |
| 3.8 | `inviscid_velocity` | `TE_info`, `panel_linvortex_velocity`, `panel_constsource_velocity` | Velocity at arbitrary point | [ ] |
| 3.9 | `build_wake` | `space_geom`, `inviscid_velocity`, `norm2` | Trace wake streamline | [ ] |
| 3.10 | `calc_force` | `get_ueinv`, `get_cp`, `get_uk`, + viscous closures | Compute cl, cm, cd | [ ] |
| 3.11 | `solve_inviscid` | `init_thermo`, `build_gamma`, `calc_force` | Top-level inviscid solve | [ ] |

**Notes:**
- `build_gamma` creates AIC matrix (N+1)x(N+1) — `# NOTE: sparse candidate`
- `calc_force` has both inviscid and viscous paths; inviscid-only path works without closures
- File: `src/inviscid.jl`, `src/thermo.jl`

**MILESTONE: Inviscid analysis works.** Verify NACA 2412 at alpha=2 gives correct cl, cm, cdpi.

---

## Phase 4: Closure Relations

~30 get_* functions returning `(value, linearization)`. Ordered by dependency.

| Step | Function | Deps | Description | Status |
|------|----------|------|-------------|--------|
| 4.1 | `get_H` | — | Shape parameter delta*/theta | [ ] |
| 4.2 | `get_Hw` | — | Wake gap shape parameter | [ ] |
| 4.3 | `get_uq` | — | Equilibrium velocity gradient | [ ] |
| 4.4 | `upwind` | — | Upwind averaging utility | [ ] |
| 4.5 | `get_Mach2` | `get_uk` | Local Mach squared | [ ] |
| 4.6 | `get_Hk` | `get_H`, `get_Mach2` | Kinematic shape parameter | [ ] |
| 4.7 | `get_Ret` | `get_Mach2`, `get_uk` | Re_theta | [ ] |
| 4.8 | `get_upw` | `get_Hk` | Upwind factor | [ ] |
| 4.9 | `get_Hs` | `get_Hk`, `get_Ret` | H* kinetic energy shape param | [ ] |
| 4.10 | `get_Hss` | `get_Mach2`, `get_Hk` | H** density shape param | [ ] |
| 4.11 | `get_cf` | `get_Hk`, `get_Ret`, `get_Mach2` | Skin friction | [ ] |
| 4.12 | `get_de` | `get_Hk` | BL thickness | [ ] |
| 4.13 | `get_rho` | `get_Mach2` | Density | [ ] |
| 4.14 | `get_Us` | `get_Hs`, `get_Hk`, `get_H` | Wall slip velocity | [ ] |
| 4.15 | `get_damp` | `get_Hk`, `get_Ret` | Amplification rate | [ ] |
| 4.16 | `get_cfxt` | `get_cf` | cf * x/theta | [ ] |
| 4.17 | `get_cfutstag` | `get_Hk` | cf*ue*theta at stagnation | [ ] |
| 4.18 | `get_cdutstag` | `get_Hk` | cDi*ue*theta at stagnation | [ ] |
| 4.19 | `get_cDi_lam` | `get_Hk`, `get_Ret` | Laminar dissipation | [ ] |
| 4.20 | `get_cDi_lamwake` | `get_Hk`, `get_Hs`, `get_Ret` | Laminar wake dissipation | [ ] |
| 4.21 | `get_cDi_outer` | `get_Hs`, `get_Us` | Outer dissipation | [ ] |
| 4.22 | `get_cDi_lamstress` | `get_Hs`, `get_Us`, `get_Ret` | Laminar stress dissipation | [ ] |
| 4.23 | `get_cDi_turbwall` | `get_cf`, `get_Hk`, `get_Hs`, `get_Us`, `get_Ret` | Turbulent wall dissipation | [ ] |
| 4.24 | `get_cDi` | all cDi sub-functions | Total dissipation | [ ] |
| 4.25 | `get_cDixt` | `get_cDi` | cDi * x/theta | [ ] |
| 4.26 | `get_cteq` | `get_Hk`, `get_Hs`, `get_H`, `get_Ret`, `get_Us` | Equilibrium ctau | [ ] |
| 4.27 | `get_cttr` | `get_cteq`, `get_Hk` | Transition ctau | [ ] |

**Notes:**
- Every function returns (value, value_U) where value_U is the Jacobian w.r.t. state U
- These are the closure relations from PDF Appendix B (Eqs. 45-79)
- Files: `src/closures_shape.jl` (4a-4c), `src/closures_friction.jl` (4d), `src/closures_dissipation.jl` (4e-4f)

---

## Phase 5: BL Residuals & Initialization

| Step | Function | Deps | Description | Status |
|------|----------|------|-------------|--------|
| 5.1 | `build_param` | structs | Create parameter struct for a surface | [ ] |
| 5.2 | `station_param` | structs | Update param for specific station | [ ] |
| 5.3 | `stagnation_state` | — | Extrapolate state to stagnation | [ ] |
| 5.4 | `thwaites_init` | — | Thwaites laminar init correlation | [ ] |
| 5.5 | `residual_station` | many get_* | BL residual at a station | [ ] |
| 5.6 | `residual_transition` | `residual_station`, `get_damp`, `get_cttr` | Transition station residual | [ ] |
| 5.7 | `store_transition` | — | Store transition location | [ ] |
| 5.8 | `march_amplification` | `get_damp`, `upwind` | March amplification equation | [ ] |
| 5.9 | `update_transition` | `march_amplification`, `get_cttr` | Find new transition | [ ] |

**Notes:**
- `residual_station` is the most complex — calls ~15 closure functions
- File: `src/boundary_layer.jl`

---

## Phase 6: Viscous Matrices & Coupling

| Step | Function | Deps | Description | Status |
|------|----------|------|-------------|--------|
| 6.1 | `identify_surfaces` | — | Split airfoil into lower/upper/wake | [ ] |
| 6.2 | `set_wake_gap` | `TE_info` | TE dead-air thickness | [ ] |
| 6.3 | `calc_ue_m` | panel influence fns | Build ue-mass sensitivity D matrix | [ ] |
| 6.4 | `rebuild_ue_m` | — | Rebuild after stagnation moves | [ ] |
| 6.5 | `wake_sys` | `TE_info`, `get_cttr` | Wake first-node residual | [ ] |
| 6.6 | `wake_init` | `wake_sys` | Initialize first wake point | [ ] |

**Notes:**
- `calc_ue_m` builds D matrix (N+Nw)x(N+Nw) — `# NOTE: sparse candidate`
- Also builds sigma_m, ue_sigma internally
- File: `src/coupling.jl`

---

## Phase 7: Coupled Solver

| Step | Function | Deps | Description | Status |
|------|----------|------|-------------|--------|
| 7.1 | `init_boundary_layer` | `stagnation_state`, `thwaites_init`, `residual_station` | March BL from stagnation | [ ] |
| 7.2 | `stagpoint_move` | `identify_surfaces`, `rebuild_ue_m` | Move stagnation point | [ ] |
| 7.3 | `build_glob_sys` | `build_param`, `stagnation_state`, `residual_station`, `wake_sys`, `residual_transition` | Assemble residuals + Jacobian | [ ] |
| 7.4 | `clalpha_residual` | `get_ueinvref` | cl-constraint residual | [ ] |
| 7.5 | `solve_glob` | `get_ueinv`, `clalpha_residual` | Solve linear system | [ ] |
| 7.6 | `update_state` | `get_Hk` | Under-relaxed Newton update | [ ] |
| 7.7 | `solve_coupled` | `build_glob_sys`, `calc_force`, `solve_glob`, `update_state`, `stagpoint_move`, `update_transition` | Newton iteration loop | [ ] |
| 7.8 | `solve_viscous` | everything | Top-level viscous solve | [ ] |

**Notes:**
- `build_glob_sys`: R_U Jacobian 3N^tot x 4N^tot — `# NOTE: sparse candidate`
- `solve_glob`: global solve — `# NOTE: sparse candidate` (Python uses spsolve)
- File: `src/solver.jl`

**MILESTONE: Full viscous analysis.** Verify NACA 2412, Re=1e6, alpha=2, Ma=0.4:
cl=0.4889, cd=0.00617, cm=-0.0501

---

## Phase 8: Post-Processing & Polish

| Step | Function | Description | Status |
|------|----------|-------------|--------|
| 8.1 | `get_distributions` | Extract BL distributions | [ ] |
| 8.2 | `mgeom_flap` | Flap deployment | [ ] |
| 8.3 | `mgeom_addcamber` | Camber addition | [ ] |
| 8.4 | `mgeom_derotate` | Derotate chord line | [ ] |
| 8.5 | `ping_test` / `check_ping` | FD derivative verification | [ ] |
| 8.6 | Plotting (defer) | Plots.jl / Makie.jl | [ ] |

---

## File Layout

```
src/
  JFoil.jl              # module, includes, exports
  types.jl              # mutable structs (Phase 0)
  utils.jl              # utilities (Phase 0)
  spline.jl             # spline utilities (Phase 1b)
  naca.jl               # NACA coordinate generation (Phase 1a)
  geometry.jl           # set_coords, make_panels (Phase 1a,1c)
  panels.jl             # panel influence functions (Phase 2)
  thermo.jl             # thermodynamics (Phase 3a)
  inviscid.jl           # inviscid solver (Phase 3b-d)
  closures_shape.jl     # shape parameter closures (Phase 4a-4c)
  closures_friction.jl  # friction & amplification closures (Phase 4d)
  closures_dissipation.jl # dissipation & shear stress closures (Phase 4e-4f)
  boundary_layer.jl     # BL residuals, transition (Phase 5)
  coupling.jl           # ue_m matrices, surface ID (Phase 6)
  solver.jl             # coupled Newton solver (Phase 7)
  postprocess.jl        # distributions, geometry mods (Phase 8)

test/
  runtests.jl
  test_utils.jl
  test_spline.jl
  test_naca.jl
  test_geometry.jl
  test_panels.jl
  test_thermo.jl
  test_inviscid.jl
  test_closures.jl      # covers all three closures_*.jl
  test_boundary_layer.jl
  test_coupling.jl
  test_solver.jl
  test_postprocess.jl
```

## Implementation Notes

- **Sparse candidates** (use dense + comment for CUDA.jl):
  1. `build_gamma` — AIC matrix (N+1)x(N+1)
  2. `calc_ue_m` — D, sigma_m, ue_sigma matrices
  3. `build_glob_sys` — R_U, R_x Jacobians
  4. `solve_glob` — global linear solve

- **Array indexing**: Python 0-based -> Julia 1-based (MATLAB uses 1-based too, cross-reference helps)

- **Return tuples**: Python `(value, linearization)` -> Julia `(value, linearization)` or named tuples

- **In-place mutation**: Functions modifying M should use Julia convention `function_name!(M, ...)` with `!`

- **Python splines**: `scipy.interpolate.CubicSpline` -> Julia `Interpolations.jl` or manual implementation
