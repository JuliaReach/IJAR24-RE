using ProbabilityBoundsAnalysis
const PBA = ProbabilityBoundsAnalysis
using ProbabilisticReachability, LazySets
using Plots, LaTeXStrings
using IntervalArithmetic: (..)
import PyPlot
using Plots: plot
using Plots.Measures

MODELs = "Lorenz"

include("$(MODELs).jl")
using .Lorenz # runs the model
MODEL = Lorenz

plot_dir = "$(@__DIR__)/plots"
if !isdir(plot_dir)
    mkdir(plot_dir)
end

name = "$plot_dir/$(MODELs).txt"

# -----------------------------------------------------------
# Plot reachable states for the variables x and y
fig = Plots.plot(MODEL.sol; vars=(1, 2), xlab=L"x", ylab=L"y")
savefig(fig, "$plot_dir/$(MODELs)_x_y")

# Plot reachable states for the variable z vs time
fig = Plots.plot(MODEL.sol; vars=(0, 3), xlab=L"t", ylab=L"z")
savefig(fig, "$plot_dir/$(MODELs)_t_z")

# FIXME:
#  - Add failure domain in the same plot.
#  - Overlap two different cases in the same plot?
#  - Plot random trajectories inside the flowpipes.

# -----------------------------------------------------------
# confidence plots

colors = colormap("RdBu", length(MODEL.c); mid=0.5)

fig = Plots.plot(; xlab=L"t", ylab=L"x")
Plots.plot!(fig, MODEL.sol(19 .. 19.3); vars=(0, 1), lw=0.0, alpha=0.9, c=:black, lab=1.0)
for (ci, vi, coli) in zip(MODEL.c, MODEL.confidences, colors)
    Plots.plot!(fig, ci(19 .. 19.3); vars=(0, 1), lw=0.0, alpha=0.9, c=coli)
end
savefig(fig, "$plot_dir/$(MODELs)_confidence")

fig = Plots.plot(; xlab=L"x", ylab=L"y")
Plots.plot!(fig, MODEL.sol(19 .. 19.3); vars=(1, 2), lw=0.0, alpha=0.9, c=:black, lab=1.0)
for (ci, vi, coli) in zip(MODEL.c, MODEL.confidences, colors)
    Plots.plot!(fig, ci(19 .. 19.3); vars=(1, 2), lw=0.0, alpha=0.9, c=coli)
end
savefig(fig, "$plot_dir/$(MODELs)_confidence2D")

fig = Plots.plot(; xlab=L"x", ylab=L"y")
Plots.plot!(fig, MODEL.sol(0 .. 10); vars=(1, 3), lw=0.0, alpha=0.9, c=:black, lab=1.0)
for (ci, vi, coli) in zip(MODEL.c, MODEL.confidences, colors)
    Plots.plot!(fig, ci(0 .. 10); vars=(1, 3), lw=0.0, alpha=0.9, c=coli)
end
savefig(fig, "$plot_dir/$(MODELs)_confidence2D_start.pdf")
savefig(fig, "$plot_dir/$(MODELs)_confidence2D_start")

fig = Plots.plot(; xlab=L"x(t)", ylab=L"z(t)", xtickfontsize=14, ytickfontsize=14,
                 xguidefontsize=20, yguidefontsize=20, right_margin=1.5cm)
for (i, s) in enumerate(MODEL.c)
    plot!(fig, s(0 .. 17.8); vars=(1, 3), lw=0, color=colors[i], lc=colors[i], alpha=0.9)
end
clims = (0, 1)
p_all = scatter!([0], [0]; zcolor=[NaN], clims=clims, label="", c=:RdBu, colorbar_title="",
                 label_fontsize=22, background_color_subplot=:transparent, lw=0.0,
                 markerstrokecolor=:transparent, markersize=0, subplot=1)
annotate!(31, 25, text(L"\alpha", 22))
savefig(fig, "$plot_dir/$(MODELs)_confidence2D_cols.pdf")
savefig(fig, "$plot_dir/$(MODELs)_confidence2D_cols")

fig = Plots.plot(; xlab=L"t", ylab=L"z(t)", xtickfontsize=14, ytickfontsize=14, xguidefontsize=20,
                 yguidefontsize=20, right_margin=3cm)
for (i, s) in enumerate(MODEL.c)
    plot!(fig, s(17 .. 20); vars=(0, 3), lw=0, xlims=(17, 20), color=colors[i], lc=colors[i],
          alpha=0.9, right_margin=1cm)
end
clims = (0, 1)
p_all = scatter!([0], [0]; zcolor=[NaN], clims=clims, label="", c=:RdBu, colorbar_title="",
                 label_fontsize=20, background_color_subplot=:transparent, lw=0.0,
                 markerstrokecolor=:transparent, markersize=0, subplot=1)
annotate!(20.7, -170, text(L"\alpha", 22))
savefig(fig, "$plot_dir/$(MODELs)_confidence2D_tz.pdf")
savefig(fig, "$plot_dir/$(MODELs)_confidence2D_tz")

# -----------------------------------------------------------
# p-box plots

PBA.plot(MODEL.pbs[1]; fontsize=30)
PyPlot.xlabel(L"X"; fontsize=30)
PyPlot.savefig("$plot_dir/pbox1.pdf")
PyPlot.savefig("$plot_dir/pbox1")

PBA.plot(MODEL.pbs[2]; fontsize=30)
PyPlot.xlabel(L"Y"; fontsize=30)
PyPlot.savefig("$plot_dir/pbox2.pdf")
PyPlot.savefig("$plot_dir/pbox2")

PBA.plot(MODEL.pbs[3]; fontsize=30)
PyPlot.xticks([-300, -200, -100, 0])
PyPlot.xlabel(L"Z"; fontsize=30)
PyPlot.savefig("$plot_dir/pbox3.pdf")
PyPlot.savefig("$plot_dir/pbox3")
PyPlot.close("all")

m1 = mass(MODEL.pbs[2], interval(-Inf, 0))
m2 = mass(MODEL.pbs[2], interval(100, 200))
m3 = mass(MODEL.pbs[2], interval(300, Inf))

open(name, "w") do io
    println("[-Inf, 0] : $m1")
    println(io, "[-Inf, 0] : $m1")
    println("[100, 200] : $m2")
    println(io, "[100, 200] : $m2")
    println("[300, Inf] : $m3")
    return println(io, "[300, Inf] : $m3")
end
