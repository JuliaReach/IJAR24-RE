using ReachabilityAnalysis, Plots, ProbabilityBoundsAnalysis, ProbabilisticReachability
import PyPlot
using Plots: plot
using Plots.PlotMeasures
using LaTeXStrings

plot_dir = "$(@__DIR__)/plots"
if !isdir(plot_dir)
    mkdir(plot_dir)
end

global BENCHMARK_MODE
if !@isdefined(BENCHMARK_MODE)
    BENCHMARK_MODE = true  # run full benchmarks
end

println("==================")
println("Duffing oscillator")
println("==================")

# ==========================================
# Version with parametric uncertainty
# ==========================================

@taylorize function duffing!(du, u, p, t)
    local α = -1.0
    local δ = 0.3
    local γ = 0.37

    x, v, β = u
    f = γ * cos(ω * t)

    du[1] = u[2]
    du[2] = -α * x - δ * v - β * x^3 + f
    return du[3] = zero(u[3])
end

ω = 1.2
T = 2 * pi / ω

failure_dom = Ball2([-0.3, 0.35], 0.18)
U_approx = overapproximate(failure_dom, 0.01)

function initial_distribution(; NSTEPS)
    ProbabilityBoundsAnalysis.setSteps(NSTEPS)
    x1 = beta(2 .. 3, 3 .. 4) .* 0.2 .+ 0.9
    x2 = beta(7 .. 8, 2 .. 3) .* 0.2 .- 0.1
    x3 = U(1 - 0.01, 1 + 0.01)

    c = GaussianCopula([1 -0.8 0; -0.8 1 0; 0 0 1])

    ## multi-dimensional p-box
    return U0 = c([x1, x2, x3])
end

# Version without probabilities
# ----------------------------------------

X0 = Singleton([1.0, 0.0, 0.0]) ⊕ Hyperrectangle([0.0, 0, 1], [0.1, 0.1, 0.01])
prob = @ivp(x' = duffing!(x), x(0) ∈ X0, dim = 3);
sol = solve(prob; tspan=(0.0, 10 * T), alg=TMJets21a());

cols = colormap("RdBu", length(sol); mid=0.5)

fig = Plots.plot(; xlab=L"x(t)", ylab=L"y(t)", xtickfontsize=14, ytickfontsize=14,
                 xguidefontsize=20, yguidefontsize=20, right_margin=1.2cm, bottom_margin=0.2cm)
for (i, s) in enumerate(sol)
    plot!(fig, s; vars=(1, 2), lw=0.8, color=cols[i], alpha=0.8)
end

plot!(fig, failure_dom; color=:grey, lw=1.2)
clims = (0, tspan(sol).hi)
p_all = scatter!([0], [0]; zcolor=[NaN], clims=clims, label="", c=:RdBu, colorbar_title="",
                 background_color_subplot=:transparent, lw=0.0, markerstrokecolor=:transparent,
                 markersize=0, subplot=1)
annotate!(2.4, 0.05, Plots.text(L"time [s]", 20; rotation=90))
Plots.savefig(fig, "$plot_dir/duffing_no_prob.pdf")
Plots.savefig(fig, "$plot_dir/duffing_no_prob")

# Version with beta distributions: confidence plot
# ---------------------------------------------------

U0 = initial_distribution(; NSTEPS=BENCHMARK_MODE ? 200 : 10)
prob = @ivp(x' = duffing!(x), x(0) ∈ U0, dim = 3);
sol = solve(prob; tspan=(0.0, 10 * T), alg=TMJets21a());
conf_thresholds = BENCHMARK_MODE ? [1.0, 0.95, 0.8, 0.6, 0.4, 0.2, 0.001] : [1, 0.6]
confs = confidence(sol, conf_thresholds);

cols = colormap("RdBu", length(confs); mid=0.5)

fig = Plots.plot(; xlab=L"x(t)", ylab=L"y(t)", xtickfontsize=14, ytickfontsize=14,
                 xguidefontsize=20, yguidefontsize=20, right_margin=1.2cm, bottom_margin=0.2cm)

for (i, s) in enumerate(confs)
    plot!(fig, s; vars=(1, 2), color=cols[i], alpha=0.3, lw=1.2, lc=:match)
end
clims = (0, 1)
plot!(fig, failure_dom; color=:grey, lw=1.2)
p_all = scatter!([0], [0]; zcolor=[NaN], clims=clims, label="", c=:RdBu,
                 background_color_subplot=:transparent, lw=0.0, markerstrokecolor=:transparent,
                 markersize=0, subplot=1)
annotate!(2.4, 0.05, Plots.text(L"\alpha", 22))

Plots.savefig("$plot_dir/duffing_confidence.pdf")
Plots.savefig("$plot_dir/duffing_confidence")

# Version with beta distributions: pbox computations
# ----------------------------------------------------

Unew = initial_distribution(; NSTEPS=BENCHMARK_MODE ? 50 : 10)

pbs = get_pbox(sol, 20, Unew)

PyPlot.plot(pbs[1])
PyPlot.savefig("$plot_dir/duffing_pba1.pdf")
PyPlot.savefig("$plot_dir/duffing_pba1")

PyPlot.plot(pbs[2])
PyPlot.savefig("$plot_dir/duffing_pba2.pdf")
PyPlot.savefig("$plot_dir/duffing_pba2")

PyPlot.plot(pbs[3])
PyPlot.savefig("$plot_dir/duffing_pba3.pdf")
PyPlot.savefig("$plot_dir/duffing_pba3")

PyPlot.close("all")
