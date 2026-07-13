using OrdinaryDiffEq
using LinearAlgebra
using Statistics
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

# ── RHS of the 3D system (in-place for DifferentialEquations) ──────
function innovation!(du, u, p::Params, t)
    K, R, N = u
    Y = N * (1 + p.A * cpow(K / N, p.α)) * R
    du[1] = p.s * Y - p.δ * K
    du[2] = tanh(20 * R) * R * (1 - R) - p.β * Y
    du[3] = -p.dn * N + p.ϕ * Y
end

# ── Callback to stop integration if state blows up or goes unphysical
function make_callback()
    condition(u, t, integrator) = (u[2] < 1e-6 || u[2] > 1.5 ||
                                   u[3] < 1e-6 || u[1] < -1.0 ||
                                   norm(u) > 1e6)
    DiscreteCallback(condition, terminate!)
end

# ── Classify a single IC: does it converge to the equilibrium? ─────
function classify_ic(u0, p::Params, u_eq;
    T=100.0, tol=0.1)
    prob = ODEProblem(innovation!, u0, (0.0, T), p)
    sol = solve(prob, Vern9(); reltol=1e-9, abstol=1e-9,
                callback=make_callback(), maxiters=1e7)

    u_final = sol.u[end]
    dist = norm(u_final .- u_eq)

    # Converged if final state is within tol of equilibrium
    return dist < tol ? 1 : -1
end

# ── Find equilibrium numerically ───────────────────────────────────
using NLsolve

function find_equilibrium(p::Params)
    guesses = [
        [1.0, 0.5, 1.0],
        [5.0, 0.3, 2.0],
        [0.5, 0.8, 0.5],
        [10.0, 0.1, 5.0],
        [0.1, 0.9, 0.1],
    ]
    function f!(F, u)
        du = similar(u)
        innovation!(du, u, p, 0.0)
        F .= du
    end
    for u0 in guesses
        result = nlsolve(f!, u0, autodiff=:forward, method=:trust_region,
                         ftol=1e-12, iterations=1000)
        if converged(result)
            K, R, N = result.zero
            # Reject trivial/unphysical equilibria
            if K > 1e-4 && R > 1e-4 && R < 1 - 1e-4 && N > 1e-4
                return result.zero
            end
        end
    end
    error("Could not find nontrivial equilibrium")
end

# ── Sweep over (K, R) grid with N fixed ────────────────────────────
function compute_basins(;
    p = default_params(),
    N0 = 1.0,
    K_range = range(0.01, 8.0, length=200),
    R_range = range(0.05, 0.4, length=200),
    T = 100.0,
    tol = 0.5,
)
    # Find the equilibrium for these parameters
    u_eq = find_equilibrium(p)
    println("  Equilibrium: K=$(round(u_eq[1]; digits=3)), R=$(round(u_eq[2]; digits=3)), N=$(round(u_eq[3]; digits=3))")

    nK = length(K_range)
    nR = length(R_range)
    basins = zeros(Int, nK, nR)
    count = 0

    for j in 1:nR
        for i in 1:nK
            u0 = [K_range[i], R_range[j], N0]

            prob = ODEProblem(innovation!, u0, (0.0, T), p)
            sol = solve(prob, Vern9(); reltol=1e-9, abstol=1e-9,
                        callback=make_callback(), maxiters=1e7,
                        saveat=T/100)
            t_final = sol.t[end]

            # A trajectory cut short by the callback (resource collapse R→0, blow-up,
            # N→0, …) left the physical domain and is NOT in the basin. Scoring it by a
            # tail-average over a truncated trajectory produces spurious "converged"
            # speckle in e.g. the high-K/high-R corner — see diagnose_ic / termination_map.
            if sol.retcode != ReturnCode.Success || t_final < T * (1 - 1e-6)
                basins[i, j] = -1
                count += 1
                continue
            end

            # Average distance over last 20% of trajectory (robust for spirals)
            n_pts = length(sol.u)
            tail_start = max(1, n_pts - div(n_pts, 5))
            avg_dist = mean(norm(u .- u_eq) for u in sol.u[tail_start:end])
            basins[i, j] = avg_dist < tol ? 1 : -1

            count += 1
        end
    end

    return K_range, R_range, basins, u_eq
end

# ── Diagnostics ────────────────────────────────────────────────────
# Integrate a single IC and return everything compute_basins discards, so you can
# see WHY a cell was classified as it was: retcode, whether the terminate! callback
# cut the run short (`early`), the final state, and the tail-average distance.
# Returns the full `sol` too, so you can lines!(sol.t, ...) to view the trajectory.
function diagnose_ic(p::Params, u0; u_eq = find_equilibrium(p), T = 100.0, tol = 0.5)
    sol = solve(ODEProblem(innovation!, u0, (0.0, T), p), Vern9();
                reltol=1e-9, abstol=1e-9, callback=make_callback(),
                maxiters=1e7, saveat=T/100)
    n_pts = length(sol.u)
    tail_start = max(1, n_pts - div(n_pts, 5))
    tailavg = mean(norm(u .- u_eq) for u in sol.u[tail_start:end])
    early = sol.retcode != ReturnCode.Success || sol.t[end] < T * (1 - 1e-6)
    (; retcode = sol.retcode, t_end = sol.t[end], early, n_pts,
       u_final = sol.u[end], tailavg, class = (early ? -1 : (tailavg < tol ? 1 : -1)),
       max_norm = maximum(norm(u) for u in sol.u),
       min_R = minimum(u[2] for u in sol.u), sol)
end

# Grid version: alongside the basin classification, record which cells had their
# integration cut short by the callback, and save a 2-panel PNG. The "odd" basin
# regions show up as early-terminations (they leave the physical domain).
function termination_map(; p = default_params(), N0 = 1.0,
        K_range = range(0.1, 40.0, length=120),
        R_range = range(0.01, 1.0, length=120),
        T = 100.0, tol = 0.5, filename = "termination_map.png")
    u_eq = find_equilibrium(p)
    nK, nR = length(K_range), length(R_range)
    cls  = zeros(nK, nR)
    term = zeros(nK, nR)
    for j in 1:nR, i in 1:nK
        d = diagnose_ic(p, [K_range[i], R_range[j], N0]; u_eq=u_eq, T=T, tol=tol)
        cls[i, j]  = d.class
        term[i, j] = d.early ? 1.0 : 0.0
    end
    println("converged=$(count(==(1), cls))  diverged=$(count(==(-1), cls))  " *
            "terminated_early=$(Int(sum(term)))/$(length(cls))")

    Rv, Kv = collect(R_range), collect(K_range)
    fig = Figure(size=(1000, 480))
    ax1 = Axis(fig[1, 1]; xlabel="R", ylabel="K", aspect=1,
        title="Classification (blue=converged)")
    heatmap!(ax1, Rv, Kv, cls'; colormap=cgrad([:gray70, :dodgerblue]), colorrange=(-1, 1))
    scatter!(ax1, [u_eq[2]], [u_eq[1]]; color=:red, marker=:star5, markersize=12)
    ax2 = Axis(fig[1, 2]; xlabel="R", ylabel="K", aspect=1,
        title="Terminated early by callback (yellow=yes)")
    heatmap!(ax2, Rv, Kv, term'; colormap=cgrad([:navy, :yellow]), colorrange=(0, 1))
    save(filename, fig; px_per_unit=2)
    println("Saved $filename")
    return (; K_range, R_range, cls, term, u_eq)
end

# ── Plot ───────────────────────────────────────────────────────────
function main(;
    A = 5.0, ϕ = 0.25, N0 = 3.5,
    K_range = range(0.1, 20.0, length=200),
    R_range = range(0.01, 0.25, length=200),
    T = 100.0, tol = 0.5,
)
    p = Params(A, 1.0, 0.5, 0.5, 0.05, 0.2, ϕ)

    println("Computing basins of attraction (direct integration)...")
    println("  A=$A, ϕ=$ϕ, N0=$N0")
    K_range, R_range, basins, u_eq = compute_basins(;
        p=p, N0=N0, K_range=K_range, R_range=R_range, T=T, tol=tol)

    n_conv = count(==(1), basins)
    n_div  = count(==(-1), basins)
    println("Done. Converged: $n_conv, Diverged: $n_div")

    fig = Figure(size=(650, 550))
    ax = Axis(fig[1, 1];
        xlabel="R", ylabel="K",
        title="Basins of attraction (direct integration)",
        aspect=1,
        xlabelsize=18, ylabelsize=18,
        xticklabelsize=15, yticklabelsize=15,
        titlesize=18)

    hm = heatmap!(ax, collect(R_range), collect(K_range), Float64.(basins');
        colormap=cgrad([:gray70, :dodgerblue]),
        colorrange=(-1.0, 1.0))

    # Mark the equilibrium
    scatter!(ax, [u_eq[2]], [u_eq[1]]; color=:red, markersize=10, marker=:star5)

    # Parameter and equilibrium annotation
    eq_str = "K̄=$(round(u_eq[1]; digits=2)), R̄=$(round(u_eq[2]; digits=2)), N̄=$(round(u_eq[3]; digits=2))"
    param_str = "A=$A, ϕ=$ϕ, N0=$N0\n$eq_str"
    text!(ax, 0.02, 0.98; text=param_str, align=(:left, :top),
        space=:relative, fontsize=12, color=:black)

    # Legend
    #elem_conv = [MarkerElement(color=:dodgerblue, marker=:circle, markersize=15)]
    #elem_div  = [MarkerElement(color=:gray70, marker=:circle, markersize=15)]
    #elem_eq   = [MarkerElement(color=:red, marker=:star5, markersize=15)]
    #Legend(fig[1, 2],
    #    [elem_conv, elem_div, elem_eq],
    #    ["Converges to eqbm", "Diverges", "Equilibrium"])

    save("basins_direct.png", fig, px_per_unit=2)
    println("Saved basins_direct.png")
    display(fig)
    return fig
end

# ── Compute basins for multiple N values ───────────────────────────
function compute_compare_N(;
    A = 5.0, ϕ = 0.25, β = 0.05, s = 1.0,
    N_values = [1.0, 2.3, 3.5],
    Kmin = 0.1, Kmax = 20.0,
    Rmin = 0.01, Rmax = 0.25,
    npts = 200,
    T = 100.0, tol = 0.5,
)
    K_range = range(Kmin, Kmax, length=npts)
    R_range = range(Rmin, Rmax, length=npts)
    p = Params(A, s, 0.5, 0.5, β, 0.2, ϕ)
    u_eq = find_equilibrium(p)
    println("Equilibrium: K=$(round(u_eq[1]; digits=2)), R=$(round(u_eq[2]; digits=2)), N=$(round(u_eq[3]; digits=2))")

    all_basins = []
    for N0 in N_values
        println("Computing basins for N0=$N0...")
        _, _, basins, _ = compute_basins(;
            p=p, N0=N0, K_range=K_range, R_range=R_range, T=T, tol=tol)
        push!(all_basins, basins)
    end

    return (;  # named tuple
        A, ϕ, β, s, N_values, K_range, R_range, u_eq, all_basins
    )
end

# ── Plot precomputed basins ────────────────────────────────────────
function plot_compare_N(data, filename="basins_compare_N.pdf";
    colors = [(:red, 0.5), (:green, 0.5), (:blue, 0.5)],
)
    (; A, ϕ, β, s, N_values, K_range, R_range, u_eq, all_basins) = data
    Rv = collect(R_range)
    Kv = collect(K_range)

    fig = Figure(size=(650, 550), figure_padding=(4, 14, 4, 4))
    ax = Axis(fig[1, 1];
        xlabel="Natural infrastructure (R)", ylabel="Built infrastructure (K)",
        title="A=$A, ϕ=$ϕ, β=$β, s=$s",
        width=500, height=500,   # fixed square size lets resize_to_layout! trim the slack
        xlabelsize=27, ylabelsize=27,
        xticklabelsize=22.5, yticklabelsize=22.5,
        titlesize=27)

    legend_elems = []
    legend_labels = String[]

    for (k, N0) in enumerate(N_values)
        basin_mask = Float64.(all_basins[k]' .== 1)
        basin_alpha = replace(basin_mask, 0.0 => NaN)

        # rasterize: in PDF output each heatmap cell is otherwise a separate quad,
        # and anti-aliased seams between quads show as hairline white lines.
        heatmap!(ax, Rv, Kv, basin_alpha;
            colormap=cgrad([colors[k][1], colors[k][1]]),
            colorrange=(0.0, 1.0), alpha=colors[k][2], rasterize=10)

        # Swatch color = fill color pre-blended with white at the fill's alpha,
        # i.e. what a single basin looks like over the white axis background.
        # An opaque swatch avoids double-blending with the legend background.
        c = Makie.to_color(colors[k][1])
        α = colors[k][2]
        swatch = Makie.RGBf(α * c.r + (1 - α), α * c.g + (1 - α), α * c.b + (1 - α))
        push!(legend_elems, [MarkerElement(color=swatch, marker=:circle, markersize=15)])
        push!(legend_labels, "N₀ = $N0")
    end

    scatter!(ax, [u_eq[2]], [u_eq[1]]; color=:black, markersize=20, marker=:star5)
    push!(legend_elems, [MarkerElement(color=:black, marker=:star5, markersize=15)])
    push!(legend_labels, "Equilibrium:\n K=$(round(u_eq[1]; digits=3)), R=$(round(u_eq[2]; digits=3)), N=$(round(u_eq[3]; digits=3))")

    axislegend(ax, legend_elems, legend_labels, "Initial population (N₀)";
        position=:rt, margin=(10, 10, 10, 10), backgroundcolor=(:white, 0.8), framecolor=:gray70, labelsize=18,
        titlesize=19.5)

    xlims!(ax, (first(R_range), 1.0))
    ylims!(ax, (first(K_range), last(K_range)))

    resize_to_layout!(fig)          # shrink the figure to fit content → trims whitespace
    save(filename, fig, px_per_unit=2)
    println("Saved $filename")
    #display(fig)
    return fig
end

data1 = compute_compare_N(A=2, ϕ=0.25, β=0.1, s=1, N_values=[1, 2, 3], Rmax=1, Kmax=40)
data2 = compute_compare_N(A=4, ϕ=0.25, β=0.1, s=1, N_values=[1, 2, 3], Rmax=1, Kmax=40)
data3 = compute_compare_N(A=4, ϕ=0.25, β=0.05, s=1, N_values=[1, 2, 3], Rmax=1, Kmax=40)
data4 = compute_compare_N(A=4, ϕ=0.15, β=0.05, s=1, N_values=[1, 2, 3], Rmax=1, Kmax=40)


plot_compare_N(data1, "figure_3a.pdf")
plot_compare_N(data2, "figure_3b.pdf")
plot_compare_N(data3, "figure_3c.pdf")
plot_compare_N(data4, "figure_3d.pdf")