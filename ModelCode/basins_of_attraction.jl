using DynamicalSystems
using OrdinaryDiffEq
using CairoMakie

# ── Parameters ──────────────────────────────────────────────────────
struct Params
    A::Float64
    s::Float64
    δ::Float64
    α::Float64
    β::Float64
    dn::Float64
    ϕ::Float64
end

default_params() = Params(3.0, 1.0, 0.5, 0.5, 0.05, 0.2, 0.1)

# ── Clipped power (matches xppaut pow) ──────────────────────────────
cpow(x, a) = x > 0 ? x^a : zero(x)

# ── RHS of the 3D system ───────────────────────────────────────────
# u = [K, R, N], p = Params
function innovation_rule(u, p, t)
    K, R, N = u
    Y = N * (1 + p.A * cpow(K / N, p.α)) * R
    dK = p.s * Y - p.δ * K
    dR = tanh(20 * R) * R * (1 - R) - p.β * Y
    dN = -p.dn * N + p.ϕ * Y
    return SVector(dK, dR, dN)
end

# ── Compute basins of attraction on a (K, R) slice ─────────────────
function compute_basins(;
    p = default_params(),
    N0 = 1.0,            # fixed N for the slice
    K_range = range(0.01, 15.0, length=200),
    R_range = range(0.01, 0.99, length=200),
    T_att = 500.0,        # integration time for attractor finding
    Δt = 0.1,
)
    # Initial condition (will be overwritten by the grid)
    u0 = SVector(1.0, 0.5, N0) 

    # Define the dynamical system
    ds = CoupledODEs(innovation_rule, u0, p; diffeq=(alg=Vern9(), reltol=1e-9, abstol=1e-9))

    # 3D grid: the recurrence algorithm needs resolution in all dimensions
    N_range = range(0.1, 5.0, length=50)
    grid = (K_range, R_range, N_range)
    mapper = AttractorsViaRecurrences(ds, grid;
        consecutive_recurrences=1000, consecutive_lost_steps=100, sparse=false)

    basins, attractors = basins_of_attraction(mapper)

    # Extract the 2D (K, R) slice at the N index closest to N0
    N_idx = argmin(abs.(collect(N_range) .- N0))
    basins_2d = basins[:, :, N_idx]

    return K_range, R_range, basins_2d, attractors
end

# ── Plot ───────────────────────────────────────────────────────────
function main(;
    A = 5.0, ϕ = 0.25, N0 = 1.0,
    K_range = range(0.01, 8.0, length=200),
    R_range = range(0.05, 0.4, length=200),
)
    p = Params(A, 1.0, 0.5, 0.5, 0.05, 0.2, ϕ)

    println("Computing basins of attraction...")
    println("  A=$A, ϕ=$ϕ, N0=$N0")
    K_range, R_range, basins, attractors = compute_basins(;
        p=p, N0=N0, K_range=K_range, R_range=R_range)
    println("Done. Found $(length(attractors)) attractor(s).")

    for (k, v) in attractors
        println("  Attractor $k: $(round.(v[end]; digits=4))")
    end

    fig = Figure(size=(650, 550))
    ax = Axis(fig[1, 1];
        xlabel="K", ylabel="R",
        title="Basins of attraction",
        aspect=1,
        xlabelsize=18, ylabelsize=18,
        xticklabelsize=15, yticklabelsize=15,
        titlesize=18)

    hm = heatmap!(ax, collect(K_range), collect(R_range), basins';
        colormap=:Set1_4)

    # Parameter annotation
    param_str = "A=$A, ϕ=$ϕ, N0=$N0\ns=$(p.s), δ=$(p.δ), α=$(p.α), β=$(p.β), dₙ=$(p.dn)"
    text!(ax, 0.02, 0.98; text=param_str, align=(:left, :top),
        space=:relative, fontsize=12,
        color=:black, font=:regular)

    # Attractor legend
    n_att = length(attractors)
    if n_att > 0
        att_str = join(["Attr $k: $(length(v)) pts" for (k, v) in attractors], "\n")
        text!(ax, 0.98, 0.98; text=att_str, align=(:right, :top),
            space=:relative, fontsize=12,
            color=:black, font=:regular)
    end

    save("basins_of_attraction.png", fig, px_per_unit=2)
    println("Saved basins_of_attraction.png")
    display(fig)
    return fig
end

main()
