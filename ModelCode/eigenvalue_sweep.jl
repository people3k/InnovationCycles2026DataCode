using ForwardDiff
using LinearAlgebra
using NLsolve
using CairoMakie
using Statistics
using OrdinaryDiffEq

# ── Parameters ──────────────────────────────────────────────────────
Base.@kwdef struct Params
    A::Float64  = 3.0
    s::Float64  = 1.0
    δ::Float64  = 0.5    # d in the ode file
    α::Float64  = 0.5    # al
    β::Float64  = 0.05   # be
    dn::Float64 = 0.2
    ϕ::Float64  = 0.2    # phi
end

# ── Clipped power (matches xppaut pow) ──────────────────────────────
cpow(x, a) = x > 0 ? x^a : zero(x)

# ── RHS of the 3D system ───────────────────────────────────────────
# u = [K, R, N]
function rhs(u, p::Params)
    K, R, N = u
    Y = N * (1 + p.A * cpow(K / N, p.α)) * R
    dK = p.s * Y - p.δ * K
    dR = tanh(20*R) * R * (1 - R) - p.β * Y
    dN = -p.dn * N + p.ϕ * Y
    return [dK, dR, dN]
end

# ── Find equilibrium numerically ───────────────────────────────────
function find_equilibrium(p::Params; u0=[1.0, 0.5, 1.0])
    f!(F, u) = (F .= rhs(u, p))
    result = nlsolve(f!, u0, autodiff=:forward, method=:trust_region,
                     ftol=1e-12, iterations=1000)
    if converged(result)
        return result.zero
    else
        return nothing
    end
end

# Try multiple initial conditions to find a nontrivial equilibrium
function find_nontrivial_equilibrium(p::Params)
    guesses = [
        [1.0, 0.5, 1.0],
        [0.1, 0.9, 0.1],
        [5.0, 0.3, 2.0],
        [10.0, 0.1, 5.0],
        [0.5, 0.8, 0.5],
    ]
    for u0 in guesses
        u★ = find_equilibrium(p; u0=u0)
        if u★ !== nothing
            K, R, N = u★
            # Accept only physically meaningful: K≥0, 0<R<1, N>0
            if K > 1e-8 && R > 1e-8 && R < 1 - 1e-8 && N > 1e-8
                return u★
            end
        end
    end
    return nothing
end

# ── Jacobian via ForwardDiff ───────────────────────────────────────
function jacobian(u, p::Params)
    ForwardDiff.jacobian(v -> rhs(v, p), u)
end

# ── Sweep over (ϕ, A) and classify eigenvalues ────────────────────
function sweep(;
    ϕ_range  = range(0.001, 0.3, length=300),
    A_range  = range(0.01, 6.0, length=300),
    base_params = Params()
)
    nϕ = length(ϕ_range)
    nA = length(A_range)

    # 0 = no valid equilibrium, 1 = all real, 2 = complex pair
    classification = zeros(Int, nA, nϕ)
    # |Im(λ)| / |Re(λ)| for the complex eigenvalue pair (overshoot severity)
    overshoot_ratio = fill(NaN, nA, nϕ)
    # Equilibrium R value for diagnostics
    eq_R = fill(NaN, nA, nϕ)

    # Store equilibria for continuation seeding
    last_eq = Vector{Union{Nothing, Vector{Float64}}}(nothing, nA)

    for (j, ϕ) in enumerate(ϕ_range)
        for (i, A) in enumerate(A_range)
            p = Params(; A=A, ϕ=ϕ,
                        s=base_params.s, δ=base_params.δ,
                        α=base_params.α, β=base_params.β,
                        dn=base_params.dn)

            # Try continuation from neighbors first, then fall back to standard guesses
            candidates = Vector{Float64}[]
            # Try previous A neighbor
            if i > 1 && last_eq[i-1] !== nothing
                u = find_equilibrium(p; u0=last_eq[i-1])
                if u !== nothing && u[1] > 1e-8 && u[2] > 1e-8 && u[2] < 1 - 1e-8 && u[3] > 1e-8
                    push!(candidates, u)
                end
            end
            # Try previous ϕ neighbor
            if j > 1 && last_eq[i] !== nothing
                u = find_equilibrium(p; u0=last_eq[i])
                if u !== nothing && u[1] > 1e-8 && u[2] > 1e-8 && u[2] < 1 - 1e-8 && u[3] > 1e-8
                    push!(candidates, u)
                end
            end
            # Standard guesses
            for u0 in [[1.0, 0.5, 1.0], [0.1, 0.9, 0.1], [5.0, 0.3, 2.0],
                        [10.0, 0.1, 5.0], [0.5, 0.8, 0.5]]
                u = find_equilibrium(p; u0=u0)
                if u !== nothing && u[1] > 1e-8 && u[2] > 1e-8 && u[2] < 1 - 1e-8 && u[3] > 1e-8
                    push!(candidates, u)
                end
            end

            # Among all valid equilibria, reject saddles (real positive eigenvalue)
            # and pick the one with highest R (the "upper" equilibrium)
            u★ = nothing
            best_R = -Inf
            for u in candidates
                res = norm(rhs(u, p))
                res > 1e-8 && continue
                J_test = jacobian(u, p)
                λ_test = eigvals(J_test)
                # Reject saddles: any real eigenvalue that is positive
                is_saddle = any(abs(imag(λ)) < 1e-10 && real(λ) > 1e-10 for λ in λ_test)
                if !is_saddle && u[2] > best_R
                    best_R = u[2]
                    u★ = u
                end
            end

            if u★ === nothing
                classification[i, j] = 0
                last_eq[i] = nothing
                continue
            end
            last_eq[i] = u★
            eq_R[i, j] = u★[2]

            # Already verified as fixed point above
            J = jacobian(u★, p)
            λs = eigvals(J)

            has_complex = any(abs(imag(λ)) > 1e-10 for λ in λs)

            # Find the complex pair's real part specifically
            re_complex = NaN
            im_complex = NaN
            for λ in λs
                if abs(imag(λ)) > 1e-10
                    re_complex = real(λ)
                    im_complex = abs(imag(λ))
                    break
                end
            end

            if has_complex && re_complex > 0
                classification[i, j] = 3  # complex unstable (past Hopf)
                overshoot_ratio[i, j] = log10(im_complex / (abs(re_complex) + 0.001) + 0.001)
            elseif has_complex
                classification[i, j] = 2  # complex stable (oscillatory attractor)
                overshoot_ratio[i, j] = log10(im_complex / (abs(re_complex) + 0.001) + 0.001)
            else
                classification[i, j] = 1  # all real (stable node)
                overshoot_ratio[i, j] = NaN
            end
        end
    end

    return ϕ_range, A_range, classification, overshoot_ratio, eq_R
end

# ── Histogram equalization: map each value to its empirical quantile in [0,1] ──
# Spreads color resolution by density so the bunched mid-range gets full contrast.
function equalize(field)
    out = fill(NaN, size(field))
    sv = sort(filter(!isnan, vec(field)))
    n = length(sv)
    n == 0 && return out
    for idx in eachindex(field)
        x = field[idx]
        isnan(x) && continue
        out[idx] = (searchsortedlast(sv, x) - 0.5) / n
    end
    return out
end

# ── Example trajectories: standard params & init, three A values ────
# Init (K=0, R=1, N=0.1) from innovation_with_pop.ode; total time trimmed to 125
# (the .ode uses 200) so the transient/overshoot fills the panels.
function example_trajectories(; A_values=[0.25, 1.25, 3.0],
                                u0=[0.0, 1.0, 0.1], T=125.0, base=Params())
    map(A_values) do A
        p = Params(; A=A, s=base.s, δ=base.δ, α=base.α, β=base.β, dn=base.dn, ϕ=base.ϕ)
        prob = ODEProblem((u, pp, t) -> rhs(u, pp), u0, (0.0, T), p)
        sol = solve(prob, Vern9(); reltol=1e-9, abstol=1e-9, saveat=T/400)
        (; A, sol)
    end
end

# ── Run and plot ───────────────────────────────────────────────────
function main()
    println("Running eigenvalue sweep...")
    ϕ_range, A_range, class, overshoot, eq_R = sweep()
    n0 = count(==(0), class)
    n1 = count(==(1), class)
    n2 = count(==(2), class)
    n3 = count(==(3), class)
    nn = count(isnan, overshoot)
    println("Sweep complete.")
    println("  class=0 (no eqbm): $n0")
    println("  class=1 (real):     $n1")
    println("  class=2 (cplx stable): $n2")
    println("  class=3 (cplx unstable): $n3")
    println("  NaN in overshoot: $nn")
    println("Plotting...")

    fig = Figure(size=(1560, 450), figure_padding=8)
    colgap!(fig.layout, 10)

    # Panel 1: Real vs Complex eigenvalues (scatter-based for legend)
    ax1 = Axis(fig[1, 1];
        xlabel="Reproductive Efficiency (ϕ)", ylabel="Technological Efficiency (A)",
        title="Eigenvalue type at equilibrium",
        aspect=1, xlabelsize=18, ylabelsize=18, xticklabelsize=15, yticklabelsize=15, titlesize=18)

    ϕv = collect(ϕ_range)
    Av = collect(A_range)

    # Eigenvalue-type field as a categorical heatmap (gap-free, replacing the old
    # markersize=2 scatter). Same colormap drives the legend colorbar below.
    eig_colors = cgrad([:gray70, :dodgerblue, :red, :magenta], categorical=true)
    heatmap!(ax1, ϕv, Av, class'; colormap=eig_colors, colorrange=(-0.5, 3.5))

    # Unstable (class 3) points — also overlaid on Panel 2 below
    ϕ_unst, A_unst = Float64[], Float64[]
    for (j, ϕ) in enumerate(ϕv), (i, A) in enumerate(Av)
        class[i, j] == 3 && (push!(ϕ_unst, ϕ); push!(A_unst, A))
    end

    # Categorical "colorbar" legend for eigenvalue type (mirrors panels 2 & 3)
    Colorbar(fig[1, 2];
        colormap = eig_colors,
        limits = (0, 4),
        ticks = (0.5:1.0:3.5,
                 ["none", "real\nstable", "complex\nstable", "complex\nunstable"]),
        ticklabelsize = 11,
        width = 12)

    # Panel 2: Overshoot severity with unstable overlay
    ax2 = Axis(fig[1, 3];
        xlabel="Reproductive Efficiency (ϕ)", ylabel="Technological Efficiency (A)",
        title="Overshoot severity",
        aspect=1, xlabelsize=18, ylabelsize=18, xticklabelsize=15, yticklabelsize=15, titlesize=18)

    valid = sort(filter(!isnan, overshoot[:]))
    oeq = equalize(overshoot)              # values → empirical quantile in [0,1]
    hm2 = heatmap!(ax2, ϕv, Av, oeq';
        colormap=:turbo, colorrange=(0, 1))
    scatter!(ax2, ϕ_unst, A_unst; color=:magenta, markersize=2)
    # Color axis is quantile-spaced; relabel ticks with the underlying log₁₀ values
    cb_q = [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]
    cb_vals = isempty(valid) ? cb_q : quantile(valid, cb_q)
    Colorbar(fig[1, 4], hm2;
        label="log₁₀(|Im(λ)| / (|Re(λ)| + 0.001))  (equalized)",
        ticks=(cb_q, [string(round(v, digits=2)) for v in cb_vals]),
        width=12)

    # Panel 3 (DISABLED — kept for reference): Equilibrium R̄ diagnostic heatmap
    # ax3 = Axis(fig[1, 5];
    #     xlabel="ϕ", ylabel="A",
    #     title="Equilibrium R̄",
    #     aspect=1, xlabelsize=18, ylabelsize=18, xticklabelsize=15, yticklabelsize=15, titlesize=18)
    # hm3 = heatmap!(ax3, ϕv, Av, eq_R';
    #     colormap=:viridis)
    # Colorbar(fig[1, 6], hm3; label="R̄", width=12)

    # Panel 3 region: two narrow example-trajectory panels (K vs t, N vs t) for the
    # standard parameter set at three A values (init & total time from the .ode file).
    trajs = example_trajectories()
    traj_colors = [:steelblue, :darkorange, :firebrick]
    traj_gl = fig[1, 5] = GridLayout()
    axK = Axis(traj_gl[1, 1]; xlabel="Scaled Time", ylabel="Human-made infrastructure (K)", title="K(t)",
        xlabelsize=16, ylabelsize=16, xticklabelsize=11, yticklabelsize=11, titlesize=16)
    axN = Axis(traj_gl[1, 2]; xlabel="Scaled Time", ylabel="Human infrastructure (N)", title="N(t)",
        xlabelsize=16, ylabelsize=16, xticklabelsize=11, yticklabelsize=11, titlesize=16)
    for (k, tr) in enumerate(trajs)
        lines!(axK, tr.sol.t, [u[1] for u in tr.sol.u]; color=traj_colors[k], linewidth=2)
        lines!(axN, tr.sol.t, [u[3] for u in tr.sol.u]; color=traj_colors[k], linewidth=2)
    end
    colgap!(traj_gl, 8)
    colsize!(fig.layout, 5, Aspect(1, 1.0))

    # Shared A-value legend for both K(t) and N(t): a categorical colorbar to the
    # right, using the same convention as Panel 1's eigenvalue-type colorbar.
    nA = length(trajs)
    Colorbar(fig[1, 6];
        colormap = cgrad(traj_colors[1:nA], categorical=true),
        limits = (0, nA),
        ticks = (collect(0.5:1.0:nA), ["A = $(tr.A)" for tr in trajs]),
        ticklabelsize = 11,
        width = 12)

    # Match axis ranges across all panels
    ϕ_lims = (first(ϕv), last(ϕv))
    A_lims = (first(Av), last(Av))
    for ax in [ax1, ax2]
        xlims!(ax, ϕ_lims)
        ylims!(ax, A_lims)
    end

    save("eigenvalue_sweep.pdf", fig) 
    println("Saved eigenvalue_sweep.pdf")

    #display(fig)
    return fig
end

main()
