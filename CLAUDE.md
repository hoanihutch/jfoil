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
