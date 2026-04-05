# JFoil - Translation Rules

## Translation Workflow (per function)

1. **Run in Python** -- execute the function in `mfoil.py` with a known input to get reference output
2. **Generate Julia tests** -- write tests in Julia that assert the expected output from step 1
3. **Translate to Julia** -- write the Julia implementation
4. **Test** -- run the Julia tests, iterate until they pass
5. **Move on** -- only proceed to the next function once tests pass

## Code Rules

- **No sparse arrays** -- use dense `Matrix{Float64}` / `Vector{Float64}` throughout. Add a `# NOTE: sparse candidate` comment wherever a sparse array would be natural (for future CUDA.jl compatibility). Sparse arrays are problematic on GPU.
- **Primary reference**: `original/mfoil.py` (Python -> Julia is most natural translation path)
- **Cross-reference**: `original/mfoil.m` (MATLAB) when `mfoil.py` is ambiguous -- it is the canonical version
- **Do NOT reference**: `original/mfoil_notex.py` -- it contains computational bugs
- **Reference doc**: `original/REFERENCE.md` has the full architecture, data structures, function list, and validation data
- **Translation plan**: `TRANSLATION_PLAN.md` has the step-by-step schedule with status checkboxes

## Project Structure

```
src/
  JFoil.jl              # module definition, includes, exports
  types.jl              # all mutable structs (Geom, Panel, Oper, Isol, Vsol, Glob, Post, Param, Mfoil)
  utils.jl              # vprint, norm2, dist
  spline.jl             # quadseg, spline2d, splineval, splinetan, spline_curvature
  naca.jl               # naca_points (NACA coordinate generation)
  geometry.jl           # set_coords, make_panels, clear_solution
  panels.jl             # panel_info, panel_linvortex_*, panel_constsource_*, panel_linsource_*
  thermo.jl             # init_thermo, get_cp, get_uk
  inviscid.jl           # build_gamma, get_ueinv, get_ueinvref, stagpoint_find, inviscid_velocity, build_wake, calc_force, solve_inviscid
  closures_shape.jl     # get_H, get_Hw, get_Mach2, get_Hk, get_Hs, get_Hss, get_Hk, get_Ret, get_rho, get_de
  closures_friction.jl  # get_cf, get_cfxt, get_cfutstag, get_Us, get_upw, upwind, get_uq, get_damp
  closures_dissipation.jl # get_cDi, get_cDi_lam, get_cDi_lamwake, get_cDi_outer, get_cDi_lamstress, get_cDi_turbwall, get_cDixt, get_cdutstag, get_cteq, get_cttr
  boundary_layer.jl     # build_param, station_param, stagnation_state, thwaites_init, residual_station, residual_transition, store_transition, march_amplification, update_transition
  coupling.jl           # identify_surfaces, set_wake_gap, calc_ue_m, rebuild_ue_m, wake_sys, wake_init
  solver.jl             # init_boundary_layer, stagpoint_move, build_glob_sys, clalpha_residual, solve_glob, update_state, solve_coupled, solve_viscous
  postprocess.jl        # get_distributions, mgeom_flap, mgeom_addcamber, mgeom_derotate

test/
  runtests.jl           # master test runner
  test_utils.jl
  test_spline.jl
  test_naca.jl
  test_geometry.jl
  test_panels.jl
  test_thermo.jl
  test_inviscid.jl
  test_closures.jl      # tests all three closures_*.jl files
  test_boundary_layer.jl
  test_coupling.jl
  test_solver.jl
  test_postprocess.jl
```

## Julia Conventions

- Mutating functions use `!` suffix: `solve_inviscid!(M)`, `build_gamma!(M, alpha)`
- Return tuples for (value, linearization) pairs from closure functions
- Use `@views` for array slices to avoid allocations
- No sparse arrays -- dense only, mark sparse candidates with `# NOTE: sparse candidate`
