# Functional Spec: Streamline Boundary Layer Solver for 3D Coupling

**Date:** 2026-04-07
**Branch:** master
**Status:** Draft

## 1. Background

JFoil currently implements a **tightly coupled 2D inviscid-viscous solver**. The inviscid panel method and the integral boundary layer equations are solved simultaneously via Newton iteration, with displacement thickness feeding back into the inviscid solution through source-distribution coupling matrices (`ue_m`, `ue_sigma`, `sigma_m`).

The goal is to enable JFoil's boundary layer solver to work with an **external 3D fluid solver** (e.g. a panel method, Euler, or RANS solver). In this mode, the 3D solver provides the pressure distribution along surface streamlines, and JFoil marches the integral boundary layer equations along each streamline independently — from a stagnation point downstream to the trailing edge. The BL solver returns displacement thickness, momentum thickness, skin friction, and transition location, which the 3D solver can use to update its solution. The two solvers iterate until convergence.

This is the **streamline analogy** (or quasi-3D) approach: the 2D integral boundary layer equations remain valid along a streamline in 3D, provided the streamline curvature is moderate and crossflow effects are secondary. It is the standard approach used in production 3D panel-method codes (e.g. VSAERO, QUADPAN, PAN AIR + BL strips).

### Relevant code paths

| File | Role | Reusable? |
|------|------|-----------|
| `src/boundary_layer.jl` | BL residuals, transition, amplification marching | Core logic reusable |
| `src/closures_shape.jl` | Shape parameter closures (H, Hk, Hs, Hss, etc.) | Fully reusable |
| `src/closures_friction.jl` | Skin friction, amplification rate | Fully reusable |
| `src/closures_dissipation.jl` | Dissipation closures | Fully reusable |
| `src/thermo.jl` | Compressibility / thermodynamics | Fully reusable |
| `src/solver.jl` | `init_boundary_layer!`, `solve_coupled!`, `build_glob_sys!`, `solve_glob!`, `update_state!` | **Not reusable** — tightly coupled to 2D panel method |
| `src/coupling.jl` | `calc_ue_m!`, `ue_sigma`, `sigma_m` coupling matrices | **Not reusable** — specific to 2D panel method |
| `src/inviscid.jl` | Panel method, wake, stagnation finding | **Not needed** — replaced by external 3D solver |
| `src/types.jl` | `Mfoil` struct aggregating all state | **Not reusable as-is** — new lightweight container needed |

## 2. Current Behaviour

### 2.1 Boundary layer marching (what works correctly)

The BL solver marches three integral equations station-by-station from a stagnation point:

- **Momentum equation**: `d(theta)/dx` in terms of `cf`, `H`, `ue`
- **Shape parameter equation**: `d(H*)/dx` in terms of dissipation `cDi`, `cf`, `Hs`
- **Third equation**: amplification factor `n` (laminar) or shear-lag `sqrt(ctau)` (turbulent)

The state vector at each station is `U = [theta, delta_star, n_or_ctau, ue]` (4 components). Stations are connected by two-point backward Euler with upwinding. The residual function `residual_station` (`boundary_layer.jl:100`) computes the 3-equation residual and its 3x8 Jacobian w.r.t. the two bounding states.

### 2.2 Initialization

`init_boundary_layer!` (`solver.jl:10`) marches from stagnation on each surface:

1. Thwaites correlation sets initial `theta`, `delta_star` from stagnation velocity gradient `K = ue/x`
2. Newton iteration at the similarity station to satisfy the BL equations exactly
3. Forward march with prescribed `ue`, solving for `(theta, delta_star, n_or_ctau)` at each station
4. If near separation (`Hk > Hmax`), switches to inverse mode (prescribes `Hk`, solves for `ue`)
5. Detects transition when amplification `n > n_crit` and restarts the station as turbulent

### 2.3 Coupled Newton iteration

`solve_coupled!` (`solver.jl:631`) assembles a global residual of size `3*Nsys` (BL equations) plus `Nsys` coupling equations (ue = ue_inv + D*(ue*delta_star)), forming a `4*Nsys x 4*Nsys` system. This is specific to the 2D panel method and **does not apply** to the streamline solver.

### 2.4 Known limitations

- The `Mfoil` struct bundles all 2D state (foil geometry, panel method, inviscid solution, coupling matrices, BL state, post-processing). There is no way to use the BL solver without constructing the full 2D problem.
- Stagnation point location is found from the panel-method circulation (`gam`), not from an external pressure distribution.
- The march direction, surface splitting (upper/lower), and wake handling are all tied to the 2D airfoil topology.
- `residual_station` and the closures are pure functions of `(param, x, U, Aux)` with no dependence on `Mfoil` — they are already decoupled from the 2D solver.

## 3. Desired Outcome

When this work is complete, JFoil should provide a **streamline boundary layer solver** that:

1. **Accepts a streamline definition**: a sequence of stations with arclength coordinates `s[1:N]` and prescribed edge velocities `ue[1:N]` (or equivalently, pressure coefficients `Cp[1:N]`). The stagnation point is the first station (or indicated by `ue ≈ 0`).

2. **Marches the BL equations** from stagnation to the last station, solving the same three integral equations as the current solver (momentum, shape parameter, amplification/shear-lag), using the same closure relations.

3. **Handles transition** using the same `e^n` amplification method and `n_crit` criterion.

4. **Returns BL results** at every station:
   - Momentum thickness `theta`
   - Displacement thickness `delta_star`
   - Skin friction coefficient `cf`
   - Shape parameter `Hk`
   - Amplification factor `n` (laminar) or shear stress coefficient `ctau` (turbulent)
   - Transition location `x_tr` (if transition occurred)
   - Separation flag (if `Hk` exceeded separation threshold)

5. **Is callable per-streamline**: a single function call processes one streamline. The external 3D solver calls it once per streamline, collects results, updates its own solution, and calls again. No global Newton coupling inside JFoil — the outer loop lives in the 3D solver.

6. **Supports iterative convergence**: since the 3D solver updates `ue` between calls, the BL solver must accept a previous BL solution as an initial guess to accelerate convergence (warm start). On the first call it initializes from Thwaites; on subsequent calls it reuses the prior state and re-marches.

7. **Does not depend on `Mfoil`**: the streamline solver should work with a lightweight input structure that contains only what it needs — no panel geometry, no AIC matrix, no coupling matrices.

8. **Coexists with the existing 2D solver**: the 2D coupled solver (`solve_viscous!`) continues to work unchanged. The streamline solver is an additional entry point, not a replacement.

## 4. Issues & Challenges

### 4.1 Stagnation point initialization from external data

- **What:** The current solver finds the stagnation point from the panel-method circulation zero-crossing. In the streamline solver, the stagnation point is implicit in the edge velocity distribution — it is where `ue → 0`. The Thwaites initialization requires the velocity gradient `K = due/ds` at stagnation, which must be estimated from the first few stations of the provided `ue` distribution.
- **Why it matters:** Poor estimation of `K` leads to incorrect initial `theta` and `delta_star`, which can cause the march to diverge or produce wrong transition locations.
- **Evidence:** `thwaites_init` at `boundary_layer.jl:85` takes `K` and `nu`. The current solver computes `K = ue[2]/xi[2]` at `solver.jl:59`, assuming the first station is at or very near stagnation.

### 4.2 Inverse mode and separation

- **What:** During initialization, if `Hk` exceeds a threshold the current solver switches to inverse mode — prescribing `Hk` and solving for `ue` instead of vice versa. In the streamline solver, `ue` is prescribed by the 3D solver and cannot be adjusted locally.
- **Why it matters:** If the prescribed `ue` distribution causes separation (`Hk > ~3.8` laminar or `~2.5` turbulent), the BL march may fail or produce physically meaningless results. The 3D solver may need to be informed that separation occurred so it can adjust.
- **Evidence:** Inverse-mode switch at `solver.jl:152-168`. The coupled solver handles this gracefully because it can modify `ue` through the Newton iteration; the streamline solver cannot.

### 4.3 Wake handling

- **What:** The current solver marches a wake downstream of the trailing edge as a third "surface" with special initialization (`wake_init`, `wake_sys`) and a wake gap parameter (`wgap`). It is unclear whether the streamline solver should handle wakes.
- **Why it matters:** For 3D applications, the wake is typically handled by the 3D solver (vortex sheets, wake panels). The streamline BL solver may only need to march to the trailing edge. However, some applications may want to continue the BL into the near-wake for drag integration.
- **Evidence:** Wake initialization at `solver.jl:86-90`, wake system in `coupling.jl:wake_sys`, wake gap in `coupling.jl:set_wake_gap!`.

### 4.4 Crossflow effects are not captured

- **What:** The streamline analogy assumes the boundary layer develops independently along each streamline, ignoring crossflow (the velocity component perpendicular to the streamline within the BL). This is a fundamental limitation of the approach.
- **Why it matters:** On swept wings, crossflow instabilities can trigger transition earlier than the streamwise Tollmien-Schlichting mechanism captured by the `e^n` method. The streamline solver will underpredict drag and miss transition on swept geometries.
- **Evidence:** The current closures (`get_damp`, `get_cf`, `get_cDi`) are all 2D — they have no crossflow terms. This is inherent to the integral BL formulation and not a bug.

### 4.5 Metric terms from streamline divergence

- **What:** In 3D, streamlines converge and diverge. The integral BL equations along a streamline should include a metric (or divergence) term that accounts for the changing width of the streamtube. Without this correction, momentum and displacement thickness will be wrong on 3D bodies with significant streamline convergence/divergence (e.g. near the nose of a fuselage or wingtip).
- **Why it matters:** Ignoring the metric term is acceptable for mild 3D geometries (wings at moderate sweep) but becomes a significant source of error on highly 3D bodies.
- **Evidence:** The current momentum equation in `residual_station` (`boundary_layer.jl:267`) has no metric/divergence term — it is purely 2D.

### 4.6 Warm-start convergence

- **What:** When the 3D solver updates `ue` between iterations, the BL solver should reuse the previous solution as an initial guess rather than re-initializing from Thwaites. The current `init_boundary_layer!` supports this partially (`if !M.oper.initbl && size(M.glob.U,2) == M.glob.Nsys` at `solver.jl:19`), but the logic is coupled to the `Mfoil` state.
- **Why it matters:** Without warm start, each BL solve is a cold start from Thwaites, wasting computation and potentially giving inconsistent transition locations between iterations, slowing convergence of the outer 3D loop.
- **Evidence:** `solver.jl:19-23` shows the warm-start path exists but is tied to `Mfoil` state management.

### 4.7 Edge velocity sign convention

- **What:** The current solver uses signed edge velocities (`sgnue` at `inviscid.jl`) to distinguish upper and lower surfaces, with the stagnation point splitting them. The streamline solver receives a single streamline where `ue` is always positive (magnitude of the velocity at the BL edge).
- **Why it matters:** The closures and residual functions use `ue` magnitude internally, but the coupling and stagnation extrapolation assume specific sign conventions. The streamline solver must ensure `ue > 0` everywhere except possibly at stagnation.
- **Evidence:** `stagnation_state` at `boundary_layer.jl:31-80` uses `ue` linearly near stagnation (`ue = K*x`), which is always positive.

## 5. Scope

### In scope

- A new public function (e.g. `solve_streamline_bl`) that takes a streamline definition and operating conditions and returns BL distributions
- A lightweight input/output data structure for the streamline solver (independent of `Mfoil`)
- Reuse of existing closures, `residual_station`, `residual_transition`, and `thwaites_init` without modification
- Transition detection via the existing `e^n` method
- Warm-start capability for iterative 3D coupling
- Separation detection and reporting (flag, not inverse-mode fix)
- Tests validating the streamline solver against the existing 2D solver on known cases (e.g. NACA 2412 — the streamline solver with 2D panel-method `ue` should reproduce the 2D coupled result)

### Out of scope

- Crossflow transition or crossflow BL equations
- Streamline metric/divergence terms (future enhancement; see Issue 4.5)
- The 3D fluid solver itself
- Wake marching beyond the trailing edge
- Modifications to the existing 2D coupled solver
- 3D geometry handling, streamline tracing, or surface meshing

## 6. Success Criteria

1. **Reproduces 2D results**: Running the streamline solver with edge velocities extracted from the existing 2D coupled solution for NACA 0012 at Re=1e5, alpha=0 and NACA 2412 at Re=1e5, alpha=5 produces `delta_star`, `theta`, `cf`, and transition location within 1% of the 2D solver output at every station.

2. **Iterative convergence**: A simple test loop that (a) runs the 2D inviscid panel method for `ue`, (b) calls `solve_streamline_bl` to get `delta_star`, (c) adds mass sources `ue*delta_star` to perturb `ue`, and (d) repeats — converges (BL state change < 1e-6) within 30 iterations. This validates the solver works in an iterative coupling context.

3. **Warm start acceleration**: The iterative test from criterion 2 converges in fewer iterations with warm start enabled than with cold Thwaites re-initialization each time.

4. **Separation reporting**: When given a `ue` distribution that causes separation (e.g. strong adverse pressure gradient), the solver returns a separation flag and the station index where separation occurred, rather than crashing or producing NaN.

5. **Independence from `Mfoil`**: The streamline solver function signature does not accept or require an `Mfoil` instance. It operates on its own input structure.

6. **Existing tests pass**: All tests in `test/runtests.jl` continue to pass without modification.
