using ProbabilityBoundsAnalysis
const PBA = ProbabilityBoundsAnalysis
using ProbabilisticReachability, LazySets, IntervalArithmetic
using Plots, LaTeXStrings
import PyPlot
using Plots: plot
using Plots.Measures
import ReachabilityAnalysis
using ReachabilityAnalysis: tstart, tend

MODELs = "UnivariateOscillator"

include("$(MODELs).jl")
using .UnivariateOscillator
MODEL = UnivariateOscillator

plot_dir = "$(@__DIR__)/plots"
if !isdir(plot_dir)
    mkdir(plot_dir)
end

name = "$plot_dir/$(MODELs).txt"

# Initial distributions:
# -> precise distribution uniform(-1, 1)
# -> imprecise distribution uniform(interval(-1, 0), interval(0, 1))

# ========== p-boxes ==========

PBA.plot(MODEL.pb; plotting=false, fontsize=30)
PBA.PyPlot.xticks([-1, -0.5, 0, 0.5, 1])
PBA.savefig("$plot_dir/$(MODELs)_precise_prob.pdf")
PBA.savefig("$plot_dir/$(MODELs)_precise_prob")
PBA.gcf()

PBA.plot(MODEL.pb_imp; plotting=false, fontsize=30)
PBA.PyPlot.xticks([-1, -0.5, 0, 0.5, 1])
PBA.savefig("$plot_dir/$(MODELs)_imprecise_prob.pdf")
PBA.savefig("$plot_dir/$(MODELs)_imprecise_prob")
PBA.gcf()

# ========== Flowpipe with failure domain ==========

fig = plot(; xlab=L"t", ylab=L"x(t)", xtickfontsize=14, ytickfontsize=14,
           xguidefontsize=20, yguidefontsize=20, right_margin=1.2cm,
           bottom_margin=0.2cm)
plot!(fig, MODEL.sol; vars=(0, 1), lw=0, c=:yellow, lc=:match, linealpha=0.0)
[plot!(fig, MODEL.dt, MODEL.analytic_sol(x0).(MODEL.dt); lab="", c=:blue,
       lw=2.0) for x0 in MODEL.x0vals];
failure_set1 = cartesian_product(BallInf([MODEL.failure_time_domain], 0.07),
                                 LazySets.Interval(MODEL.failure_domain))
plot!(fig, failure_set1; c=:red, alpha=0.9)
failure_set2 = cartesian_product(BallInf([MODEL.failure_time_domain], 0.07),
                                 LazySets.Interval(MODEL.failure_domain2))
plot!(fig, failure_set2; c=:red, alpha=0.9)
fig
savefig(fig, "$plot_dir/$(MODELs)_flowpipe.pdf")
savefig(fig, "$plot_dir/$(MODELs)_flowpipe")

# ========== Flowpipe initial sets ==========
fig = plot(; xlab=L"t", ylab=L"x(t)", xtickfontsize=14, ytickfontsize=14,
           xguidefontsize=20, yguidefontsize=20, right_margin=1.2cm,
           bottom_margin=0.2cm)
NR = 20
t1 = 10.0
t2 = 11.5
plot!(fig, MODEL.sol(t1 .. t2); vars=(0, 1), lw=1.0, c=:yellow, lc=:black, linealpha=1.0, ls=:dash)
# Highlight color.
plot!(fig, MODEL.sol(t1); vars=(0, 1), lw=1.0, c=:orange, lc=:black, linealpha=1.0, ls=:dash)
@show ts = tstart(MODEL.sol(t1))
@show te = tend(MODEL.sol(t2))
dta = ts:0.0001:te
[plot!(fig, dta, MODEL.analytic_sol(x0).(dta); lab="", c=:blue, lw=2.0)
 for x0 in MODEL.x0vals[1:2:end]];
fig
savefig(fig, "$plot_dir/$(MODELs)_flowpipe_zoom.pdf")
savefig(fig, "$plot_dir/$(MODELs)_flowpipe_zoom")
#= Chosen reach-set for the example:
julia> MODEL.sol(t1)
ReachabilityAnalysis.TaylorModelReachSet{Float64, Float64}
TaylorModels.TaylorModel1{TaylorSeries.TaylorN{Float64}, Float64}
[  0.1561947947609466 x₁ + ( 0.08058023256116119 x₁) t + ( 0.08768774238565437 x₁) t² + ( 0.024658928220114718 x₁) t³ + ( 0.011188253317172145 x₁) t⁴ + ( - 0.0007740283315820886 x₁) t⁵ + ( - 0.0011413464531049775 x₁) t⁶ + ( - 0.0007586181252353169 x₁) t⁷ + ( - 0.00021188682563597408 x₁) t⁸ + ( - 1.361963773664984e-6 x₁) t⁹ + ( 2.0848295480730633e-5 x₁) t¹⁰ + [-7.59378e-15, 7.59378e-15]], [9.96683, 10.0598])

=#

# ========== Precision as a function of the discretization step ==========

plot(; xlab="discretization granularity", ylab="failure probability",
     xtickfontsize=14, ytickfontsize=14, xguidefontsize=14, yguidefontsize=14,
     right_margin=1.2cm, bottom_margin=0.2cm)
plot!(MODEL.steps_list, inf.(MODEL.p_lower_precise_uniform_tune); c=:blue, lab="interval estimate")
plot!(MODEL.steps_list, sup.(MODEL.p_lower_precise_uniform_tune); c=:blue, lab="")
plot!(MODEL.steps_list,
      repeat([MODEL.p_lower_analytic_uniform_larger_domain], length(MODEL.steps_list)); c=:red,
      lab="analytic")
savefig("$plot_dir/$(MODELs)_discretization.pdf")
savefig("$plot_dir/$(MODELs)_discretization")

# Write output file
open(name, "w") do io
    case = "$(MODELs)_p_lower_analytic_uniform"
    print(io, case, ", ")
    println(io, MODEL.p_lower_analytic_uniform)

    case = "$(MODELs)_p_lower_precise"
    print(io, case, ", ")
    println(io, MODEL.p_lower_precise)

    for (i, steps) in enumerate(MODEL.steps_list)
        case = "$(MODELs)_p_lower_precise_$steps"
        print(io, case, ", ")
        println(io, MODEL.p_lower_precise_uniform_tune[i])
    end

    case = "$(MODELs)_p_lower_imprecise"
    print(io, case, ", ")
    println(io, MODEL.p_lower_imprecise)

    case = "$(MODELs)_p_lower_analytic_uniform_larger_domain"
    print(io, case, ", ")
    println(io, MODEL.p_lower_analytic_uniform_larger_domain)

    case = "$(MODELs)_p_lower_precise_larger_domain"
    print(io, case, ", ")
    return println(io, MODEL.p_lower_precise_larger_domain)
end

# ========== Confidence plots ==========

c = colormap("RdBu", length(MODEL.conf); mid=0.5)

# Confidence plots using precise distribution
fig = Plots.plot(; xlab=L"t", ylab=L"x(t)", xtickfontsize=14, ytickfontsize=14, xguidefontsize=20,
                 yguidefontsize=20, right_margin=1.2cm, bottom_margin=0.2cm)
[Plots.plot!(fig, MODEL.conf[i]; vars=(0, 1), c=c[i], alpha=0.9, lw=0, lc=:match)
 for i in 1:length(c)]
clims = (0, 1)
p_all = scatter!([0], [0]; zcolor=[NaN], clims=clims, label="", c=:RdBu, colorbar_title="",
                 label_fontsize=22, background_color_subplot=:transparent, lw=0.0,
                 markerstrokecolor=:transparent, markersize=0, subplot=1)
annotate!(19.2, 0.03, text(L"\alpha", 22))
savefig("$plot_dir/$(MODELs)_confidence.pdf")
savefig("$plot_dir/$(MODELs)_confidence")

# Confidence plots using imprecise distribution
fig = Plots.plot(; xlab=L"t", ylab=L"x(t)", xtickfontsize=14, ytickfontsize=14, xguidefontsize=20,
                 yguidefontsize=20, right_margin=1.2cm, bottom_margin=0.2cm)
[Plots.plot!(fig, MODEL.conf_imp[i]; vars=(0, 1), c=c[i], alpha=0.9, lw=0, lc=:match)
 for i in 1:length(c)]
clims = (0, 1)
p_all = scatter!([0], [0]; zcolor=[NaN], clims=clims, label="", c=:RdBu, colorbar_title="",
                 label_fontsize=22, background_color_subplot=:transparent, lw=0.0,
                 markerstrokecolor=:transparent, markersize=0, subplot=1)
annotate!(19.2, 0.03, text(L"\alpha", 22))
savefig("$plot_dir/$(MODELs)_confidence_imprecise.pdf")
savefig("$plot_dir/$(MODELs)_confidence_imprecise")

# Confidence plots using precise distribution and consonance method
fig = Plots.plot(; xlab=L"t", ylab=L"x(t)", xtickfontsize=14, ytickfontsize=14, xguidefontsize=20,
                 yguidefontsize=20, right_margin=1.2cm, bottom_margin=0.2cm)
[Plots.plot!(fig, MODEL.conf_cons[i]; vars=(0, 1), c=c[i], alpha=0.9, lw=0, lc=:match)
 for i in 1:length(c)]
clims = (0, 1)
p_all = scatter!([0], [0]; zcolor=[NaN], clims=clims, label="", c=:RdBu, colorbar_title="",
                 label_fontsize=22, background_color_subplot=:transparent, lw=0.0,
                 markerstrokecolor=:transparent, markersize=0, subplot=1)
annotate!(19.2, 0.03, text(L"\alpha", 22))
savefig("$plot_dir/$(MODELs)_confidence_imprecise.pdf")
savefig("$plot_dir/$(MODELs)_confidence_imprecise")

# Confidence plots using imprecise distribution and consonance method
fig = Plots.plot(; xlab=L"t", ylab=L"x(t)", xtickfontsize=14, ytickfontsize=14, xguidefontsize=20,
                 yguidefontsize=20, right_margin=1.2cm, bottom_margin=0.2cm)
[Plots.plot!(fig, MODEL.conf_imp_cons[i]; vars=(0, 1), c=c[i], alpha=0.9, lw=0, lc=:match)
 for i in 1:length(c)]
clims = (0, 1)
p_all = scatter!([0], [0]; zcolor=[NaN], clims=clims, label="", c=:RdBu, colorbar_title="",
                 label_fontsize=22, background_color_subplot=:transparent, lw=0.0,
                 markerstrokecolor=:transparent, markersize=0, subplot=1)
annotate!(19.2, 0.03, text(L"\alpha", 22))
savefig("$plot_dir/$(MODELs)_confidence_imprecise.pdf")
savefig("$plot_dir/$(MODELs)_confidence_imprecise")

# explanation of algorithm
if BENCHMARK_MODE
    using ProbabilisticReachability: Hyperrectangle, mid, radius, normalize, Interval, evaluate
    function interval2box(x, T)
        return Hyperrectangle([T, mid(x)], [0.01, radius(x)])
    end
    X0 = Hyperrectangle([0.0, 0], [0.0, 1])
    T = 2 * π
    Tmax = T + 0.3
    Udom = MODEL.failure_domain
    ε = 0.05  # plotting offset
    UT = interval2box(Udom, T + ε)
    R = MODEL.sol(T)
    pb0 = get_initial_pbox(MODEL.sol)
    focal_elems, masses = split(normalize(pb0))
    X1 = focal_elems[24]
    X2 = focal_elems[26]
    X3 = focal_elems[28]
    fX1 = evaluate(R, T, X1)
    fX2 = evaluate(R, T, X2)
    fX3 = evaluate(R, T, X3)

    Plots.plot(; xlab=L"t", ylab=L"x(t)", xtickfontsize=12, ytickfontsize=12, xguidefontsize=16,
               yguidefontsize=16, right_margin=1.2cm, bottom_margin=0.2cm)
    c = :yellow
    Plots.plot!(MODEL.sol(Interval(0, Tmax)); vars=[0, 1], c=c, lc=c)
    dt = range(0, Tmax; length=1000)
    x0vals = range(-1, 1; length=25)
    [Plots.plot!(dt, MODEL.analytic_sol(x0).(dt); lab="", c=:blue, lw=2) for x0 in x0vals]
    c = :red
    Plots.plot!(UT; c=c, lc=c, lw=4, alpha=1)
    c = :black
    Plots.plot!(interval2box(X1, 0); vars=[0, 1], c=c, lc=c, lw=4, alpha=1)
    Plots.plot!(interval2box(fX1, T); vars=[0, 1], c=c, lc=c, lw=4, alpha=1)
    c = :gray
    Plots.plot!(interval2box(X1, 0); vars=[0, 1], c=c, lc=c, lw=2, alpha=1)
    Plots.plot!(interval2box(fX1, T); vars=[0, 1], c=c, lc=c, lw=2, alpha=1)
    c = :black
    Plots.plot!(interval2box(X2, 0); vars=[0, 1], c=c, lc=c, lw=4, alpha=1)
    Plots.plot!(interval2box(fX2, T); vars=[0, 1], c=c, lc=c, lw=4, alpha=1)
    c = :orange
    Plots.plot!(interval2box(X2, 0); vars=[0, 1], c=c, lc=c, lw=2, alpha=1)
    Plots.plot!(interval2box(fX2, T); vars=[0, 1], c=c, lc=c, lw=2, alpha=1)
    c = :black
    Plots.plot!(interval2box(X3, 0); vars=[0, 1], c=c, lc=c, lw=4, alpha=1)
    Plots.plot!(interval2box(fX3, T); vars=[0, 1], c=c, lc=c, lw=4, alpha=1)
    c = :green
    Plots.plot!(interval2box(X3, 0); vars=[0, 1], c=c, lc=c, lw=2, alpha=1)
    Plots.plot!(interval2box(fX3, T); vars=[0, 1], c=c, lc=c, lw=2, alpha=1)
    xlims!(-0.3, Tmax)
    ylims!(-0.6, -0.4)
    savefig("$plot_dir/$(MODELs)_algorithm_explanation.pdf")
    savefig("$plot_dir/$(MODELs)_algorithm_explanation")
end
