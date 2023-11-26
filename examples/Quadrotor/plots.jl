using ProbabilisticReachability, ProbabilityBoundsAnalysis, LazySets, IntervalArithmetic
using Plots, LaTeXStrings, Plots.PlotMeasures
const PBA = ProbabilityBoundsAnalysis

MODELs = "Quadrotor"
include("$(MODELs).jl")
using .Quadrotor
MODEL = Quadrotor

plot_dir = "$(@__DIR__)/plots"
if !isdir(plot_dir)
    mkdir(plot_dir)
end

name = "$plot_dir/$(MODELs).txt"

# -----------------------------------------------------------
# Plot reachable states for the variable x3 vs time (height)
#
fig = Plots.plot(MODEL.solz_case1; vars=(0, 3), xlab="t", ylab="x3")
Plots.savefig(fig, "$plot_dir/$(MODELs)_case1_t_x3_uniform.pdf")
Plots.savefig(fig, "$plot_dir/$(MODELs)_case1_t_x3_uniform")

fig = Plots.plot(MODEL.solz_case2; vars=(0, 3), xlab="t", ylab="x3")
Plots.savefig(fig, "$plot_dir/$(MODELs)_case2_t_x3_uniform.pdf")
Plots.savefig(fig, "$plot_dir/$(MODELs)_case2_t_x3_uniform")

# -----------------------------------------------------------
# Print failure probability results
#
pcases = ((:pfail_cond1_case1, MODEL.pfail_cond1_case1),
          (:pfail_cond1_case2, MODEL.pfail_cond1_case2),
          (:pfail_cond2_case1, MODEL.pfail_cond2_case1),
          (:pfail_cond2_case2, MODEL.pfail_cond2_case2),
          (:pfail_cond3_case1, MODEL.pfail_cond3_case1),
          (:pfail_cond3_case2, MODEL.pfail_cond3_case2),
          (:pfail_cond13_case1, MODEL.pfail_cond13_case1),
          (:pfail_cond13_case2, MODEL.pfail_cond13_case2))

# Write output file
open(name, "w") do io
    for (psym, pint) in pcases
        case = "$(MODELs)_$(psym)"
        print(io, case, ", ")
        println(io, pint)
    end
end

# ------------------------------------------------------------------------------
# Plots
#
using ProbabilisticReachability: subscript

colors = colormap("RdBu", length(MODEL.pconf); mid=0.5)

U1 = convert(Hyperrectangle, interval(0, 5) × interval(1.4, 5))
U2 = convert(Hyperrectangle, interval(1, 5) × interval(-2, 0.9))
U3 = convert(Hyperrectangle, interval(5, 5) × interval(0.98, 1.02))

vv = (0, 3) # variable x3 vs time
for (i, out) in enumerate((MODEL.conf_case1, MODEL.conf_case2))
    local fig = Plots.plot(; ylab=L"x_3(t)", xlab=L"t", xtickfontsize=14, ytickfontsize=14,
                           xguidefontsize=20, yguidefontsize=20, right_margin=1.2cm,
                           bottom_margin=0.2cm, xlims=(0.0, 5.0), ylims=(-1.0, 1.5))
    #=
    xlab=L"t",
    ylab=L"x_3",
    tickfont=font(30, "Times"), guidefontsize=45)
    xtick=[0.0, 1.0, 2.0, 3.0, 4.0, 5.0], ytick=[-1.0, -0.5, 0.0, 0.5, 1.0, 1.5],
    ,
    bottom_margin=6mm, left_margin=2mm, right_margin=4mm, top_margin=3mm,
    size=(1000, 1000))
    # \raisebox{-0.5mm}{\textcolor{white}{.}}"
    # \raisebox{2mm}{\textcolor{white}{.}}",
    =#

    [Plots.plot!(fig, out[i]; vars=vv, c=colors[i], alpha=0.9, lw=0, lc=:match)
     for i in 1:length(colors)]

    local clims = (0, 1)
    local p_all = scatter!([0], [0]; zcolor=[NaN], clims=clims, label="", c=:RdBu,
                           colorbar_title="", label_fontsize=22,
                           background_color_subplot=:transparent, lw=0.0,
                           markerstrokecolor=:transparent, markersize=0, subplot=1)
    annotate!(6.2, 0.55, text(L"\alpha", 22))

    #=
    Plots.plot!(x -> x, x -> 0.98, 0.0, 5.0, line=1.5, color=:black, linestyle=:dot, legend=nothing)
    Plots.plot!(x -> x, x -> 1.02, 0.0, 5.0, line=1.5, color=:black, linestyle=:dot, legend=nothing)
    Plots.plot!(x -> x, x -> 0.9, 0.0, 5.0, line=1.5, color=:black, linestyle=:dot, legend=nothing)
    Plots.plot!(x -> x, x -> 1.4, 0.0, 5.0, line=1.5, color=:black, linestyle=:dot, legend=nothing)
    =#

    Plots.plot!(U1; vars=vv, lw=1.5, c=:orange)
    Plots.plot!(U2; vars=vv, lw=1.5, c=:green)
    Plots.plot!([5, 5], [0.98, 1.02]; lw=10, legend=nothing, color=:black)

    annotate!(2.5, 1.7, text(L"\mathcal{U}_1", 20))
    annotate!(3, 0, text(L"\mathcal{U}_2", 20))

    start_x = 4.3
    start_y = 1.2
    annotate!(start_x - 0.2, start_y + 0.05, text(L"\mathcal{U}_3", 20))
    quiver!([start_x], [start_y]; quiver=([4.9 - start_x], [1.05 - start_y]), color=:black, lw=1.5)

    Plots.ylims!((-1, 2))
    Plots.xlims!((0, 5))

    Plots.savefig(fig, "$plot_dir/quadrotor_confidence_case$i.pdf")
    Plots.savefig(fig, "$plot_dir/quadrotor_confidence_case$i")
end

# ------------------------------------------------------------------------------
# p-box plots
#
import PyPlot

PBA.PyCall.pyimport_conda("mpl_toolkits.mplot3d", "mpl_toolkits")
art3d = PBA.PyCall.PyObject(PyPlot.art3D)
poly3d = PBA.PyCall.PyObject(PyPlot.matplotlib.collections.PolyCollection)

# case 1
fig = make_step_plots(MODEL.sol_case1, Quadrotor.BENCHMARK_MODE ? 80 : 3; dims=3:3,
                      figure=PyPlot.figure, poly3d=poly3d,
                      xlabel=PyPlot.xlabel, ylabel=PyPlot.ylabel, zlabel=PyPlot.zlabel,
                      xlim=PyPlot.xlim,
                      ylim=PyPlot.ylim);
PyPlot.zlabel("CDF"; rotation=90)

PyPlot.savefig("$plot_dir/pbox_3d_quadrotor_case1.pdf")
PyPlot.savefig("$plot_dir/pbox_3d_quadrotor_case1")

# case 2
fig = make_step_plots(MODEL.sol_case2, Quadrotor.BENCHMARK_MODE ? 80 : 3; dims=3:3,
                      figure=PyPlot.figure, poly3d=poly3d,
                      xlabel=PyPlot.xlabel, ylabel=PyPlot.ylabel, zlabel=PyPlot.zlabel,
                      xlim=PyPlot.xlim,
                      ylim=PyPlot.ylim);
PyPlot.zlabel("CDF"; rotation=90)

PyPlot.savefig("$plot_dir/pbox_3d_quadrotor_case2.pdf")
PyPlot.savefig("$plot_dir/pbox_3d_quadrotor_case2")
PyPlot.close("all")
