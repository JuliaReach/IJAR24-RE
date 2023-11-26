using ProbabilisticReachability, ProbabilityBoundsAnalysis, LaTeXStrings
import PyPlot
using PyPlot: using3D, figure, gca

global BENCHMARK_MODE
if !@isdefined(BENCHMARK_MODE)
    BENCHMARK_MODE = true  # run full benchmarks
end

ProbabilityBoundsAnalysis.setSteps(BENCHMARK_MODE ? 200 : 10)

plot_dir = "$(@__DIR__)/plots"
if !isdir(plot_dir)
    mkdir(plot_dir)
end

using3D()

ft_big = 55
ft = 44
ft_ticks = 40

fig_size = (10, 10)

x1 = beta(interval(1, 2), 3)
x2 = gamma(interval(5, 6), 2)
x3 = normal(interval(0, 0.5), interval(2, 3))

fig = figure("x1"; figsize=fig_size)
ProbabilityBoundsAnalysis.plot(x1; name="x1", fontsize=ft)
PyPlot.xlabel(L"$X_1$"; fontsize=ft)
PyPlot.title(L"$F_{X_{1}}$"; fontsize=ft_big)
PyPlot.tight_layout()
PyPlot.savefig("$plot_dir/x1.pdf")
PyPlot.savefig("$plot_dir/x1.png")

fig = figure("x2"; figsize=fig_size)
ProbabilityBoundsAnalysis.plot(x2; name="x2", fontsize=ft)
PyPlot.xlabel(L"$X_2$"; fontsize=ft)
PyPlot.title(L"$F_{X_{2}}$"; fontsize=ft_big)
PyPlot.tight_layout()
PyPlot.savefig("$plot_dir/x2.pdf")
PyPlot.savefig("$plot_dir/x2.png")

fig = figure("x3"; figsize=fig_size)
ProbabilityBoundsAnalysis.plot(x3; name="x3", fontsize=ft)
PyPlot.xticks(range(-10, 10; length=5); fontsize=ft)
PyPlot.xlabel(L"$X_3$"; fontsize=ft)
PyPlot.title(L"$F_{X_{3}}$"; fontsize=ft_big)
PyPlot.tight_layout()
PyPlot.savefig("$plot_dir/x3.pdf")
PyPlot.savefig("$plot_dir/x3.png")
PyPlot.close("all")

cor12 = -0.8;
cor13 = 0.8;
cor23 = -0.8;

C12 = GauCopula(cor12)
C13 = GauCopula(cor13)
C23 = GauCopula(cor23)

F12 = C12(x1, x2)
F13 = C12(x1, x3)
F23 = C12(x2, x3)

name = "biv_x1x2"
fig = figure(name; figsize=fig_size)
ProbabilityBoundsAnalysis.plot(F12; name=name)

PyPlot.zlabel(L"$CDF$"; fontsize=ft, rotation=90)
PyPlot.xlabel(L"$X_1$"; fontsize=ft)
PyPlot.ylabel(L"$X_2$"; fontsize=ft)

ax = gca()
PyPlot.xticks(; fontsize=ft_ticks)
PyPlot.yticks(; fontsize=ft_ticks)
#ax.zticks(fontsize = ft_ticks)
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])

for t in ax.zaxis.get_major_ticks()
    t.label.set_fontsize(ft_ticks)
    #t.label.set_fontsize("")
end
ax.axes.zaxis.set_ticklabels([])

ax.view_init(30, -60)
fig.canvas.draw()
PyPlot.title(L"$F_{X_1,X_2}$"; fontsize=ft_big)
PyPlot.tight_layout()
PyPlot.savefig("$plot_dir/$name.pdf")
PyPlot.savefig("$plot_dir/$name.png")
PyPlot.close("all")

name = "biv_x1x3"
fig = figure(name; figsize=fig_size)
ProbabilityBoundsAnalysis.plot(F13; name=name)

PyPlot.zlabel(L"$CDF$"; fontsize=ft, rotation=90)
PyPlot.xlabel(L"$X_1$"; fontsize=ft)
PyPlot.ylabel(L"$X_3$"; fontsize=ft)

ax = gca()
PyPlot.xticks(; fontsize=ft_ticks)
PyPlot.yticks(; fontsize=ft_ticks)
#ax.zticks(fontsize = ft_ticks)
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])

for t in ax.zaxis.get_major_ticks()
    t.label.set_fontsize(ft_ticks)
    #t.label.set_fontsize("")
end
ax.axes.zaxis.set_ticklabels([])

ax.view_init(25, -60)
fig.canvas.draw()
PyPlot.title(L"$F_{X_1, X_3}$"; fontsize=ft_big)
PyPlot.tight_layout()
PyPlot.savefig("$plot_dir/$name.pdf")
PyPlot.savefig("$plot_dir/$name.png")
PyPlot.close("all")

name = "biv_x2x3"
fig = figure(name; figsize=fig_size)
ProbabilityBoundsAnalysis.plot(F23; name=name)

PyPlot.zlabel(L"$CDF$"; fontsize=ft, rotation=90)
PyPlot.xlabel(L"$X_2$"; fontsize=ft)
PyPlot.ylabel(L"$X_3$"; fontsize=ft)

ax = gca()
PyPlot.xticks(; fontsize=ft_ticks)
PyPlot.yticks(; fontsize=ft_ticks)
#ax.zticks(fontsize = ft_ticks)
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])

for t in ax.zaxis.get_major_ticks()
    t.label.set_fontsize(ft_ticks)
    #t.label.set_fontsize("")
end
ax.axes.zaxis.set_ticklabels([])

ax.view_init(30, -60)
fig.canvas.draw()
PyPlot.title(L"$F_{X_2,X_3}$"; fontsize=ft_big)
PyPlot.tight_layout()
PyPlot.savefig("$plot_dir/$name.pdf")
PyPlot.savefig("$plot_dir/$name.png")
PyPlot.close("all")

######

name = "biv_x1x2_2"
fig = figure(name; figsize=fig_size)
ProbabilityBoundsAnalysis.plot(F12; name=name)

PyPlot.zlabel(L"$CDF$"; fontsize=ft, rotation=90)
PyPlot.xlabel(L"$X_1$"; fontsize=ft)
PyPlot.ylabel(L"$X_2$"; fontsize=ft)

ax = gca()
PyPlot.xticks(; fontsize=ft_ticks)
PyPlot.yticks(; fontsize=ft_ticks)
#ax.zticks(fontsize = ft_ticks)
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])

for t in ax.zaxis.get_major_ticks()
    t.label.set_fontsize(ft_ticks)
    #t.label.set_fontsize("")
end
ax.axes.zaxis.set_ticklabels([])

ax.view_init(30 + 5, 180 - 25)
fig.canvas.draw()
PyPlot.title(L"$F_{X_1, X_2}$"; fontsize=ft_big)
PyPlot.tight_layout()
PyPlot.savefig("$plot_dir/$name.pdf")
PyPlot.savefig("$plot_dir/$name.png")
PyPlot.close("all")

name = "biv_x1x3_2"
fig = figure(name; figsize=fig_size)
ProbabilityBoundsAnalysis.plot(F13; name=name)

PyPlot.zlabel(L"$CDF$"; fontsize=ft, rotation=90)
PyPlot.xlabel(L"$X_1$"; fontsize=ft)
PyPlot.ylabel(L"$X_3$"; fontsize=ft)

ax = gca()
PyPlot.xticks(; fontsize=ft_ticks)
PyPlot.yticks(; fontsize=ft_ticks)
#ax.zticks(fontsize = ft_ticks)
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])

for t in ax.zaxis.get_major_ticks()
    t.label.set_fontsize(ft_ticks)
    #t.label.set_fontsize("")
end
ax.axes.zaxis.set_ticklabels([])

ax.view_init(30 + 5, 180 - 25)
fig.canvas.draw()
PyPlot.title(L"$F_{X_1,X_3}$"; fontsize=ft_big)
PyPlot.tight_layout()
PyPlot.savefig("$plot_dir/$name.pdf")
PyPlot.savefig("$plot_dir/$name.png")
PyPlot.close("all")

name = "biv_x2x3_2"
fig = figure(name; figsize=fig_size)
ProbabilityBoundsAnalysis.plot(F23; name=name)

PyPlot.zlabel(L"$CDF$"; fontsize=ft, rotation=90)
PyPlot.xlabel(L"$X_2$"; fontsize=ft)
PyPlot.ylabel(L"$X_3$"; fontsize=ft)

ax = gca()
PyPlot.xticks(; fontsize=ft_ticks)
PyPlot.yticks(; fontsize=ft_ticks)
#ax.zticks(fontsize = ft_ticks)
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])

for t in ax.zaxis.get_major_ticks()
    t.label.set_fontsize(ft_ticks)
    #t.label.set_fontsize("")
end
ax.axes.zaxis.set_ticklabels([])

ax.view_init(30 + 5, 180 - 25)
fig.canvas.draw()
PyPlot.title(L"$F_{X_2,X_3}$"; fontsize=ft_big)
PyPlot.tight_layout()
PyPlot.savefig("$plot_dir/$name.pdf")
PyPlot.savefig("$plot_dir/$name.png")
PyPlot.close("all")

#=  The hard way

Cov =   [1 0.8 -0.8;
        0.8 1 -0.8;
        -0.8 -0.8 1];

C = GaussianCopula(Cov)

F1 = C([x1, x2, x3])

x1Range = range(x1)
x2Range = range(x2)
x3Range = range(x3)

numel = 100;

x1s = range(x1Range.lo, x1Range.hi, length = numel)
x2s = range(x2Range.lo, x2Range.hi, length = numel)
x3s = range(x3Range.lo, x3Range.hi, length = numel)

bivX1X2 = [cdf(F1,[a,b, Inf]) for a in x1s, b in x2s]
bivX1X3 = [cdf(F1,[a,Inf, b]) for a in x1s, b in x3s]
bivX2X3 = [cdf(F1,[Inf, a, b]) for a in x2s, b in x3s]

=#
