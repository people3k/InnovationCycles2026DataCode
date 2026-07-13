# Innovators Dilemma Project

## Overview
Dynamical systems model of innovation, capital (K), resources (R), and population (N).
Originally written for xppaut, now being translated to Julia for more advanced analysis.

## Key Files
- `innovation_with_pop.ode` — 3D ODE system (K, R, N) in xppaut format (the main model)
- `innovation.ode` — 2D version (K, R) without population
- `equilibrium_analysis.tex` — Analytical equilibrium derivations (illustrative only, uses simplified cases)
- `eigenvalue_sweep.jl` — 2-parameter sweep over (ϕ, A), computes Jacobian eigenvalues at equilibria, classifies real/complex/unstable, plots overshoot severity. Saves to `eigenvalue_sweep.pdf`
- `basins_direct.jl` — Basins of attraction via direct integration. Classifies ICs by whether they converge to the interior equilibrium. Has `compute_compare_N()` and `plot_compare_N()` for overlaying basins at different N values
- `basins_of_attraction.jl` — Alternative approach using DynamicalSystems.jl `AttractorsViaRecurrences` (less reliable for this system, kept for reference)

## Model Details
- **State variables**: K (capital), R (resource), N (population)
- **Output**: Y = N·(1 + A·(K/N)^α)·R
- **ODEs**: dK/dt = s·Y − δ·K; dR/dt = tanh(20R)·R·(1−R) − β·Y; dN/dt = −dₙ·N + ϕ·Y
- **Default params**: A=3, s=1, δ=0.5, α=0.5, β=0.05, dₙ=0.2, ϕ=0.1
- **Key equilibria**: saddle at (0,0,0), node at (0,1,0), interior equilibrium (the one we track)
- System undergoes Hopf bifurcation as A and ϕ increase
- Two equilibria can collide (saddle-node) due to tanh nonlinearity

## Julia Notes
- Uses ForwardDiff.jl for Jacobian computation
- NLsolve.jl for equilibrium finding — must reject saddles (real positive eigenvalue) and trivial equilibria
- Heatmap transpose: Makie `heatmap!(ax, x, y, z)` expects `z[i,j]` at `(x[i], y[j])` — watch array orientation
- `cpow(x, a) = x > 0 ? x^a : zero(x)` — clipped power matching xppaut, uses `zero(x)` for ForwardDiff compatibility
- xppaut syntax `20R` must be `20*R` in Julia
- Makie colormaps: use `Reverse(:RdBu)` not `:RdBu_r` (matplotlib syntax)
- Categorical colormaps via `cgrad` can be finicky; manual `MarkerElement` legends work better

## Workflow Tips
- `compute_compare_N()` / `plot_compare_N(data)` split lets you tweak plots without recomputing
- After editing a function in VS Code, must re-evaluate the definition before calling it
- ImageMagick available at `/opt/local/bin/convert` for cropping PNGs
