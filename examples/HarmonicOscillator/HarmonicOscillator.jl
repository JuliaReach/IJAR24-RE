# # Harmonic oscillator

# ## Model

# We consider a two-dimensional harmonic oscillator.
# The differential equation is:
#
# ```math
# \left\{ \begin{aligned}
# x'(t) &= y(t) \\
# y'(t) &= -x(t)
# \end{aligned} \right. \qquad \textrm{subject to } (x(0), y(0)) \in X_0 = [0.8, 1.2] \times [0.8, 1.2] \subset \mathbb{R}^2
# ```
# and the analytic solution is:
#
# ```math
# \left\{ \begin{aligned}
# x(t) &= x_0 \cos~ t + y_0 \sin~ t \\
# y(t) &= -x_0 \sin t + y_0 \cos t
# \end{aligned} \right.
# ```

using ReachabilityAnalysis,
      ProbabilisticReachability,
      ProbabilityBoundsAnalysis, Plots, Plots.Measures, FileIO
import PyPlot

using LaTeXStrings

const PBA = ProbabilityBoundsAnalysis

global BENCHMARK_MODE
if !@isdefined(BENCHMARK_MODE)
    BENCHMARK_MODE = true  # run full benchmarks
end

println("===================")
println("Harmonic oscillator")
println("===================")

# The corresponding Julia function describing the ODE is:

@taylorize function f!(du, u, p, t)
    du[1] = u[2]
    return du[2] = -u[1]
end

# ## Probabilistic initial-value problem

# We assume that the initial condition belongs to the p-box
# given by a uniform distribution with independent marginals,
# $(x(0), y(0)) \in U(0.8 .. 1.2) × U(0.8 .. 1.2)$.

## fix the discretization steps

ns = BENCHMARK_MODE ? 200 : 10
ProbabilityBoundsAnalysis.setSteps(ns)

x1 = beta(3, 3) .* 0.4 .+ 0.8
x2 = beta(8, 2) .* 0.4 .+ 0.8

c = GaussianCopula([1 -0.8; -0.8 1])

## multi-dimensional p-box
U0 = c([x1, x2])

## initial-value problem
prob = @ivp(u' = f!(u), u(0) ∈ U0, dim = 2);

# We can plot the cdf of the initial distribution.

##ProbabilityBoundsAnalysis.plot(U0, save=true, name="univariate_oscillator_U0")

# ## p-box propagation

# We obtain a flowpipe corresponding to the probabilistic initial-value problem
# defined above from time $t = 0$ to $T = 15$.

alg = TMJets21b(; orderT=10, orderQ=2, abstol=1e-15)

sol = solve(prob; alg=alg, T=10);

## number of computed reach-sets

length(sol)

# ## Compute the two p-boxes at t = 5
#pb = get_pbox(sol, 3.0)

u_failure = Ball2([0.0, 0.0], 1.2)
U_approx = overapproximate(u_failure, 0.01)

Plots.plot(sol; vars=(1, 2), ratio=1, xlab=L"x(t)", ylab=L"y(t)", xtickfontsize=18,
           ytickfontsize=18, xguidefontsize=22, yguidefontsize=22, right_margin=1.2cm,
           bottom_margin=0.5cm)
Plots.plot!(U_approx; xlims=(-2, 2), ylims=(-2, 2), color=:yellow)

Plots.savefig("$plot_dir/Flowpipe_with_failure.pdf")
Plots.savefig("$plot_dir/Flowpipe_with_failure")

#-

#PBA.plot(pb[1], col= "red", plotting=false)
#PBA.gcf() #hide

#-

#-

#PBA.plot(pb[1], col= "blue", plotting=false)
#PBA.gcf() #hide

#-
println("********************")
println("Running confidence plots")
println("********************")

confs = confidence(sol,
                   BENCHMARK_MODE ? [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.001] :
                   [1, 0.5])

cols = colormap("RdBu", length(confs); mid=0.5)

fig = Plots.plot(; xlab=L"x(t)", ylab=L"y(t)", ratio=1, xtickfontsize=18, ytickfontsize=18,
                 xguidefontsize=22, yguidefontsize=22, right_margin=1.2cm, bottom_margin=0.5cm)

#fig = Plots.plot(xlab="x(t)", ylab="y(t)",ratio =1)
vv = (1, 2)

for (i, s) in enumerate(confs)
    plot!(fig, s; vars=(1, 2), lw=0, color=cols[i], lc=cols[i], alpha=1)
end

Plots.plot!(U_approx; color=:yellow, xlims=(-2, 2), ylims=(-2, 2))
clims = (0, 1)
p_all = scatter!([0], [0]; zcolor=[NaN], clims=clims, label="", c=:RdBu, colorbar_title="",
                 label_fontsize=22, background_color_subplot=:transparent, lw=0.0,
                 markerstrokecolor=:transparent, markersize=0, subplot=1)

annotate!(3.3, 0.05, text(L"\alpha", 22))

Plots.savefig("$plot_dir/H2_confidence.pdf")
Plots.savefig("$plot_dir/H2_confidence")

### Failure prob

# u_failure = Ball2([0.0, 0.0], 1.2)
# U_approx = overapproximate(u_failure, 0.01)
#
#
# fig = Plots.plot(xlab=L"x(t)", ylab=L"y(t)", ratio=1)
# vv = (1, 2)
#
#
# Plots.plot!(U_approx)
#
# Plots.savefig("$plot_dir/H2_confidence_with_failure.pdf")
# Plots.savefig("$plot_dir/H2_confidence_with_failure")

### Show convergence in step size

# println("********************")
# println("Running example with varying step size")
# println("********************")
#
# using PyPlot, LaTeXStrings
# using ProbabilityBoundsAnalysis: plot
#
# ns = 100
# ProbabilityBoundsAnalysis.setSteps(ns)
#
#
# c = GaussianCopula([1 -0.8; -0.8 1])
#
# steps = [5, 10, 50, 100, 200]
#
# for s in steps
#
#     ProbabilityBoundsAnalysis.setSteps(s)
#
#     UNew = c([x1, x2])
#
#     pbNew = get_pbox(sol, 10, UNew)
#
#     fontsize = 24
#
#     PBA.plot(pbNew[1], name = "x1", fontsize = fontsize);
#     PyPlot.xlabel(L"x", fontsize = fontsize)
#     PyPlot.savefig("steps1.pdf")
#     PyPlot.savefig("steps1")
#
#     PBA.plot(pbNew[2], name = "x2", fontsize = fontsize)
#     PyPlot.xlabel(L"y", fontsize = fontsize)
#     PyPlot.savefig("steps2.pdf")
#     PyPlot.savefig("steps2")
#
# end

# Change the descritization and compute again (without solving the Reachability problem again)

println("********************")
println("Running example with varying correlation")
println("********************")

ns = BENCHMARK_MODE ? 100 : 10
ProbabilityBoundsAnalysis.setSteps(ns)

x1 = beta(3, 3) .* 0.4 .+ 0.8
x2 = beta(8, 2) .* 0.4 .+ 0.8

c = GaussianCopula([1 -0.8; -0.8 1])

UNew = c([x1, x2])

pbNew = get_pbox(sol, 10, UNew)

fontsize = 44
ticksize = 30

PBA.plot(pbNew[1]; name=L"x1", col="red", fontsize=ticksize);
PyPlot.xlabel(L"x"; fontsize=fontsize)

PBA.plot(pbNew[2]; name=L"x2", col="red", fontsize=ticksize)
PyPlot.xlabel(L"y"; fontsize=fontsize)

cors = range(-0.999, 0.999; length=5)
colors = ["red", "blue", "green", "black", "purple"]

for (i, cor) in enumerate(cors)
    local c = GaussianCopula([1 cor; cor 1])

    UNew2 = c([x1, x2])
    pbOut = get_pbox(sol, 10, UNew2)

    PBA.plot(pbOut[1]; col=colors[i], name="x1_diff", fontsize=ticksize)
    PyPlot.tight_layout()
    PyPlot.savefig("$plot_dir/X_1_pboxes_diff.pdf")
    PyPlot.savefig("$plot_dir/X_1_pboxes_diff")

    PBA.plot(pbOut[2]; col=colors[i], name="x2_diff", fontsize=ticksize)
    PyPlot.tight_layout()
    PyPlot.savefig("$plot_dir/X_2_pboxes_diff.pdf")
    PyPlot.savefig("$plot_dir/X_2_pboxes_diff")
end

PyPlot.close("all")

### How does the correlation effect the failure prob
println("********************")
println("Running failure probability calculation")
println("********************")

ns = BENCHMARK_MODE ? 50 : 10
ProbabilityBoundsAnalysis.setSteps(ns)

U_time = 8 .. 10

prob_f_baseline = failure_probability(sol, U_time, U_approx)
println("Baseline pf: $prob_f_baseline")
save("$plot_dir/harmonic_oscillator_baseline.jld2", "baseline", prob_f_baseline)

x1 = beta(3, 3) .* 0.4 .+ 0.8
x2 = beta(8, 2) .* 0.4 .+ 0.8

cors = range(-0.999, 0.999; length=BENCHMARK_MODE ? 50 : 5)

probs_f = Vector{IntervalArithmetic.Interval{Float64}}(undef, length(cors))

for (i, cor) in enumerate(cors)
    local c = GaussianCopula([1 cor; cor 1])

    UNew2 = c([x1, x2])
    probs_f[i] = failure_probability(sol, U_time, U_approx, UNew2)
    println("***********************")
    println("Correlation: $cor")
    println("failureprob: $(probs_f[i])")
end

save("$plot_dir/harmonic_oscillator_corrs.jld2", "corrs", cors)
save("$plot_dir/harmonic_oscillator_probs_f.jld2", "probs", probs_f)

PyPlot.plot(cors, sup.(probs_f); color="red")
PyPlot.plot(cors, inf.(probs_f); color="red")
PyPlot.fill_between(cors, sup.(probs_f), inf.(probs_f); color="grey", alpha=0.3)

PyPlot.xticks([-1, -0.5, 0, 0.5, 1]; fontsize=18)
PyPlot.yticks(; fontsize=18)

PyPlot.xlabel(L"\rho"; fontsize=24)
PyPlot.ylabel("failure probability"; fontsize=24)

PyPlot.tight_layout()

PyPlot.savefig("$plot_dir/varying_correlation_prob_f.pdf")
PyPlot.savefig("$plot_dir/varying_correlation_prob_f")
PyPlot.close("all")

#-
#PBA.plot(pbNew[1], col= "red", plotting=false)
#PBA.gcf() #hide
#-

#-
#PBA.plot(pbNew[1], col= "blue", plotting=false)
#PBA.gcf() #hide
#-

#=

U_time = interval(7.5, 10)

U_failure = interval(-1.5, -1) × interval(0, 0.5)

U_failure = convert(Hyperrectangle, U_failure)

#ProbabilityBoundsAnalysis.setSteps(10)
ProbabilityBoundsAnalysis.setSteps(200)
failure_prob = failure_probability(sol, U_time, U_failure, UNew)
failure_prob = failure_probability(sol, U_time, U_failure)
failure_prob_max = failure_probability(sol, U_time, U_failure, UNew, maxitive = true)

failure_prob = maximum([failure_probability(sol, ts, U_failure, UNew) for ts in mince(U_time, 10)])
=#

#-

##cdf(pb, -0.1)

#-

##PBA.plot(pb, save=true, name="univariate_oscillator_t15", plotting=false) #hide
##PBA.plot(pb, plotting=false)
##PBA.gcf() #hide

# ## Confidence flowpipes

# @time out = confidence(sol, [1.0, 0.9, 0.5, 0.1]);
#fig = Plots.plot(xlab="x(t)", ylab="y(t)", title="U0 ~ U([0.8, 1.2]^2), p ∈ [.1, .5, .9, 1.]")
#vv = (1, 2)
#Plots.plot!(fig, out[1], vars=vv, c=:green, alpha=1., lw=0.0, lab="p=1.0")
#Plots.plot!(fig, out[2], vars=vv, c=:red, alpha=1., lw=0.0, lab="p=0.9")
#Plots.plot!(fig, out[3], vars=vv, c=:blue, alpha=1., lw=0.0, lab="p=0.5")
#Plots.plot!(fig, out[4], vars=vv, c=:yellow, alpha=1., lw=0.0, lab="p=0.1", ratio=1.)
#Plots.plot!(fig, ensemble(sol), vars=(0, 1), c=:black, lab="xx")
#fig
