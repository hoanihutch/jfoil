# mfoil -> JFoil Translation Reference

## Project Overview

**mfoil** is a coupled inviscid-viscous airfoil analysis solver by Krzysztof Fidkowski (University of Michigan). It reimplements the core ideas of XFOIL in a more accessible, documented codebase.

### What it does

Computes aerodynamic properties (lift, drag, moment, pressure distribution, boundary layer quantities) for 2D airfoils at subsonic speeds. It couples:

1. **Inviscid solver** -- vortex-panel potential flow method with Kutta condition
2. **Viscous solver** -- integral boundary layer equations (momentum, shape parameter, amplification/shear-lag)
3. **Transition prediction** -- e^N method for laminar-to-turbulent transition
4. **Compressibility correction** -- Karman-Tsien for subcritical Mach numbers

## Source Files

| File | Role |
|------|------|
| `mfoil.m` (~75k tokens) | MATLAB class -- the canonical, most complete version |
| `mfoilui.m` | MATLAB GUI (not needed for translation) |
| `mfoil.py` (~63k tokens) | Python translation -- functional style, uses numpy/scipy |
| `mfoil_notex.py` | Variant of mfoil.py with TeX disabled + some bugs (not reliable) |
| `mfoil.pdf` (20 pages) | Paper describing the physics, equations, and algorithms |

## Architecture (9 data structures + solver)

- **Geom** -- airfoil geometry (coordinates, chord, name)
- **Panel (foil/wake)** -- panel mesh (nodes, arclength, tangents)
- **Oper** -- operating conditions (alpha, Re, Ma, viscous flag)
- **Isol** -- inviscid solution (AIC matrix, vortex strengths, stagnation)
- **Vsol** -- viscous solution (BL thicknesses, transition, surface indices)
- **Glob** -- global Newton system (state U=[theta, delta*, n_tilde/c_tau^1/2, u_e], Jacobian, residuals)
- **Post** -- outputs (cl, cd, cm, cp, cf, Hk distributions)
- **Param** -- solver parameters (tolerances, physical constants, closure coefficients)

## Data Structure Details

### Geom
- `chord`: chord length
- `wakelen`: wake extent (in chords)
- `npoint`: number of geometry points
- `name`: airfoil name (e.g., "NACA 0012")
- `xpoint`: point coordinates [2 x npoint]
- `xref`: moment reference point

### Panel (foil and wake)
- `N`: number of nodes
- `x`: node coordinates [2 x N]
- `s`: arclength at nodes
- `t`: tangent vectors (dx/ds, dz/ds)

### Oper
- `Vinf`: freestream velocity magnitude
- `alpha`: angle of attack (degrees)
- `rho`: density
- `cltgt`: target lift coefficient
- `givencl`: boolean (cl given instead of alpha)
- `initbl`: initialize boundary layer
- `viscous`: boolean (viscous analysis)
- `redowake`: rebuild wake after alpha changes
- `Re`: Reynolds number
- `Ma`: Mach number
- `xft`: forced transition x/c values [lower, upper]

### Isol (Inviscid Solution)
- `AIC`: aerodynamic influence coefficient matrix
- `gamref`: 0,90-deg vorticity distributions at nodes (N x 2)
- `gam`: gamma at current alpha
- `sstag`: s-location of stagnation point
- `sstag_g`, `sstag_ue`: linearizations of stagnation point
- `Istag`: node indices before/after stagnation
- `sgnue`: +/-1 signs on upper/lower surfaces
- `xi`: distance from stagnation point
- `uewi`: inviscid edge velocity in wake
- `uewiref`: 0,90-deg reference ue solutions on wake
- `xstag`: physical (x,z) location of stagnation point

### Vsol (Viscous Solution)
- `th`, `ds`: momentum and displacement thickness arrays
- `Is`: 3 arrays of surface indices (lower, upper, wake)
- `wgap`: wake gap at wake points
- `ue_m`: linearization of ue w.r.t. mass
- `sigma_m`: d(source)/d(mass) matrix
- `ue_sigma`: d(ue)/d(source) matrix
- `turb`: turbulent flag array (1=turbulent, 0=laminar)
- `xt`: transition location (xi) on current surface
- `Xt`: transition xi/x values for lower/upper surfaces

### Glob (Global System)
- `Nsys`: number of equations and states
- `U`: primary states [4 x Nsys] containing [th; ds; sa; ue] at each node
- `dU`: state update from Newton solver
- `dalpha`: angle of attack update
- `conv`: convergence flag
- `R`: residual vector [3*Nsys x 1]
- `R_U`: residual Jacobian w.r.t. states [3*Nsys x 4*Nsys]
- `R_x`: residual Jacobian w.r.t. xi values [3*Nsys x Nsys]

### Post (Post-Processing)
- `cp`, `cpi`: pressure coefficient (viscous and inviscid)
- `cl`, `cm`: lift and moment coefficients
- `cdpi`, `cd`, `cdf`, `cdp`: drag components
- `cl_ue`, `cl_alpha`: linearizations of cl
- Distribution arrays: `th`, `ds`, `sa`, `ue`, `uei`, `cf`, `Ret`, `Hk`

### Param (Solver Parameters)
- `verb`: verbosity level
- `rtol`: residual tolerance for Newton (default 1e-10)
- `niglob`: max global iterations (default 50)
- `doplot`: boolean for plotting results
- `ncrit`: critical amplification factor (default 9)
- G-beta locus parameters: `GA` (6.7), `GB` (0.75), `GC` (18.0)
- Shear stress: `CtauC` (1.8), `CtauE` (3.3)
- Thermodynamic: `gam` (1.4), `Tsrat` (0.35), `KTb`, `KTl`, `cps`
- Wake: `wakelen` (1 chord), `ewpow` (1e-5 offset)

## Solver Flow

### Inviscid Pipeline
1. `solve_inviscid()` -- main entry
2. `init_thermo()` -- set thermodynamic parameters from Ma
3. `build_gamma()` -- build AIC matrix, solve at 0 and 90 deg reference
4. Superpose reference solutions for desired alpha
5. `stagpoint_find()` -- locate stagnation point, compute xi distances
6. `calc_force()` -- compute cl, cm, cdpi from pressure integration
7. Optional: `cltrim_inviscid()` -- iterate alpha to match target cl

### Viscous Pipeline
1. `solve_viscous()` -- main entry
2. Solve inviscid first
3. `build_wake()` -- trace wake streamline from TE, create wake panels
4. `identify_surfaces()` -- split airfoil into lower/upper + wake
5. `set_wake_gap()` -- trailing edge dead-air thickness
6. `calc_ue_m()` -- build sensitivity of ue to mass (transpiration BC)
7. `init_boundary_layer()` -- march BL from stagnation on each surface
8. `stagpoint_move()` -- update stagnation location from viscous solution
9. `solve_coupled()` -- Newton iteration loop:
   - `build_glob_sys()` -- assemble residuals and Jacobian
   - `solve_glob()` -- solve linear system (with cl constraint if needed)
   - `update_state()` -- apply update with under-relaxation
   - `stagpoint_move()` -- move stagnation point
   - `update_transition()` -- re-march amplification, find new transition

### Boundary Layer Equations (per node)
Three residuals + one coupling equation:

1. **Momentum** (Eq. 9): d(theta)/d(xi) from cf, H, Hw, Me
2. **Shape parameter** (Eq. 10): d(H*)/d(xi) from cDi, cf
3. **Amplification** (Eq. 11, laminar): d(n_tilde)/d(xi) from empirical damp
   -- OR **Shear lag** (Eq. 12, turbulent): d(c_tau)/d(xi) from equilibrium ctau
4. **Edge velocity coupling** (Eq. 27): ue = ue_inv + D * (ue * delta*)

### State Vector
At each node: u_i = [theta, delta*, n_tilde or c_tau^{1/2}, u_e]

## Functions List (~120 total)

### Geometry & Paneling
- `set_coords()` / `naca_points()` -- set/generate airfoil coordinates
- `make_panels()` -- curvature-adaptive paneling
- `spline_curvature()` / `spline2d()` / `splineval()` / `splinetan()` -- spline utilities
- `space_geom()` -- geometric point spacing
- `TE_info()` / `panel_info()` -- panel geometry helpers
- `mgeom_flap()` / `mgeom_addcamber()` / `mgeom_derotate()` -- geometry modification

### Panel Influence Functions
- `panel_linvortex_stream()` / `panel_linvortex_velocity()`
- `panel_constsource_stream()` / `panel_constsource_velocity()`
- `panel_linsource_stream()` / `panel_linsource_velocity()`

### Inviscid Solver
- `build_gamma()` -- build and solve AIC system
- `get_ueinv()` / `get_ueinvref()` -- inviscid edge velocity
- `inviscid_velocity()` -- velocity at arbitrary point
- `rebuild_isol()` -- rebuild after alpha change

### Viscous Solver
- `solve_coupled()` -- main Newton loop
- `build_glob_sys()` -- assemble residuals + Jacobian
- `solve_glob()` -- solve linear system
- `update_state()` -- under-relaxed state update
- `calc_ue_m()` / `rebuild_ue_m()` -- ue-mass sensitivity matrices

### Boundary Layer
- `init_boundary_layer()` -- initialization by marching
- `stagnation_state()` -- extrapolate state to stagnation
- `thwaites_init()` -- Thwaites laminar correlation
- `residual_station()` -- BL residual at a station
- `residual_transition()` -- transition station residual
- `march_amplification()` -- march amplification equation
- `update_transition()` / `store_transition()` -- transition bookkeeping
- `wake_init()` / `wake_sys()` -- wake initialization

### Closure Relations (~30 "get" functions)
All return (value, linearization) for Newton Jacobian:
- `get_Hk()` -- kinematic shape parameter
- `get_Hs()` -- kinetic energy shape parameter (H*)
- `get_Hss()` -- density shape parameter (H**)
- `get_H()` -- shape parameter (delta*/theta)
- `get_Hw()` -- wake gap shape parameter
- `get_cf()` -- skin friction coefficient
- `get_cDi()` -- dissipation coefficient (with sub-components: wall, outer, lam, stress)
- `get_Us()` -- normalized wall slip velocity
- `get_damp()` -- amplification rate (dn/dxi)
- `get_Ret()` -- momentum-thickness Reynolds number
- `get_rho()` -- density from isentropic relations
- `get_Mach2()` -- local Mach number squared
- `get_cp()` / `get_uk()` -- pressure coefficient / Karman-Tsien speed
- `get_upw()` / `upwind()` -- upwinding factors
- `get_uq()` -- equilibrium velocity gradient
- `get_cttr()` / `get_cteq()` -- shear stress coefficients
- `get_de()` -- BL thickness measure
- `get_cfxt()` / `get_cfutstag()` / `get_cdutstag()` / `get_cDixt()` -- combined quantities

### Post-Processing
- `calc_force()` -- cl, cm, cd from pressure/friction integration
- `get_distributions()` -- extract BL distributions for plotting
- `plot_cpplus()` / `plot_airfoil()` / `plot_boundary_layer()` / `plot_results()`

### Testing
- `ping_test()` / `check_ping()` -- finite-difference derivative verification

## Key Notes for Julia Translation

- **Primary reference**: `mfoil.py` (Python -> Julia is most natural)
- **Cross-reference**: `mfoil.m` (MATLAB) for any ambiguity -- it is the canonical version
- **Do NOT use**: `mfoil_notex.py` has computational bugs (wrong array indexing, chord calc error)
- **Sparse matrices**: numpy/scipy sparse -> Julia `SparseArrays` + `LinearAlgebra`
- **Closure functions**: The (value, linearization) return pattern maps well to Julia tuples
- **Array indexing**: Python 0-based -> Julia 1-based (matches MATLAB, which helps)
- **Julia struct design**: The 8 data containers map naturally to Julia mutable structs
- **Plotting**: matplotlib -> Plots.jl or Makie.jl (defer to later)

## Validation Reference (from PDF Table 1)

NACA 2412, Re=1e6, alpha=2 deg, Ma=0.4, N=200:

| Output | mfoil | XFOIL |
|--------|-------|-------|
| lift coefficient | 0.4889 | 0.4910 |
| quarter-chord moment | -0.0501 | -0.0506 |
| drag coefficient | 0.00617 | 0.00618 |
| skin friction drag | 0.00421 | 0.00421 |
| upper surface transition | 0.4908c | 0.4901c |
| lower surface transition | 0.9487c | 0.9486c |

## Physical Constants (from PDF Table 2)

| Parameter | Description | Value |
|-----------|-------------|-------|
| gamma_air | Ratio of specific heats | 1.4 |
| n_crit | Critical amplification factor | 9 |
| G_A | G-beta locus A constant | 6.7 |
| G_B | G-beta locus B constant | 0.75 |
| G_C | G-beta locus C constant | 18.0 |
| eta_D | wall/wake dissipation length ratio | 0.9 in wake |
| K_lag | shear lag constant | 5.6 |
| C_tau | shear stress initialization constant | 1.8 |
| E_tau | shear stress initialization exponent | 3.3 |
| r_Su | Sutherland temperature ratio | 0.35 |
| f^w | wake gap continuation factor | 2.5 |
| d^w | wake length, in airfoil chords | 1 |
| epsilon^w | first wake point offset, in airfoil chords | 1e-5 |
