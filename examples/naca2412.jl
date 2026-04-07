using JFoil
using Plots

# Create airfoil
M = Mfoil()
naca_points!(M, "2412")
make_panels!(M, 199)

# Set operating conditions
M.oper.alpha = 2.0    # angle of attack (degrees)
M.oper.Re = 1e6       # Reynolds number
M.oper.Ma = 0.0       # Mach number

# Solve viscous
solve_viscous!(M)

# Print results
println("NACA 2412 — Viscous Solution")
println("  alpha = $(M.oper.alpha) deg")
println("  cl    = $(round(M.post.cl, digits=4))")
println("  cd    = $(round(M.post.cd, digits=6))")
println("  cm    = $(round(M.post.cm, digits=4))")
println("  cdf   = $(round(M.post.cdf, digits=5))")
println("  cdp   = $(round(M.post.cdp, digits=5))")

# Plot
fig = plot_results(M)
savefig(fig, joinpath(@__DIR__, "naca2412.png"))
println("\nPlot saved to examples/naca2412.png")
