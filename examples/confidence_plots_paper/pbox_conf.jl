using Distributions, PyPlot, LaTeXStrings, IntervalArithmetic, ProbabilityBoundsAnalysis

plot_dir = "$(@__DIR__)/plots"
if !isdir(plot_dir)
    mkdir(plot_dir)
end

ProbabilityBoundsAnalysis.setSteps(1000)

fontsize = 25
tick_font = 20

alpha1 = 0.05
alpha2 = 0.6

A = normal(-0.5 .. 0.5, 1 .. 1.8)

cut1 = ProbabilityBoundsAnalysis.cut(A, alpha1 / 2 .. 1 - alpha1 / 2)
cut2 = ProbabilityBoundsAnalysis.cut(A, alpha2 / 2 .. 1 - alpha2 / 2)

PyPlot.plot(A; fontsize=tick_font)
PyPlot.xlabel(L"x"; fontsize=fontsize)
PyPlot.ylabel("cdf"; fontsize=fontsize)
PyPlot.xticks(; fontsize=tick_font)
PyPlot.yticks(; fontsize=tick_font)

fig = gcf()
fig.set_size_inches(6.4, 4.8)

PyPlot.tight_layout()

xlimm = gca().get_xlim()
ylimm = gca().get_ylim()

xpoint = xlimm[1]
ypoint = ylimm[1]

PyPlot.plot([xpoint; cut1.hi], [1 - alpha1 / 2; 1 - alpha1 / 2], "--"; color="red")
PyPlot.plot([xpoint; cut1.lo], [alpha1 / 2; alpha1 / 2], "--"; color="red")
PyPlot.plot([cut1.lo; cut1.lo], [0; alpha1 / 2], "--"; color="red")
PyPlot.plot([cut1.hi; cut1.hi], [0; 1 - alpha1 / 2], "--"; color="red")

PyPlot.plot([xpoint; cut2.hi], [1 - alpha2 / 2; 1 - alpha2 / 2], "--"; color="blue")
PyPlot.plot([xpoint; cut2.lo], [alpha2 / 2; alpha2 / 2], "--"; color="blue")
PyPlot.plot([cut2.lo; cut2.lo], [0; alpha2 / 2], "--"; color="blue")
PyPlot.plot([cut2.hi; cut2.hi], [0; 1 - alpha2 / 2], "--"; color="blue")

PyPlot.plot([cut2.lo; cut2.hi], [-0.01; -0.01]; color="blue", linewidth=3)
PyPlot.plot([cut1.lo; cut1.hi], [-0.015; -0.015]; color="red", linewidth=3)

PyPlot.plot([xpoint - 0.1; xpoint - 0.1], [alpha2 / 2; 1 - alpha2 / 2]; color="blue", linewidth=3)
PyPlot.plot([xpoint - 0.2; xpoint - 0.2], [alpha1 / 2; 1 - alpha1 / 2]; color="red", linewidth=3)

PyPlot.savefig("$plot_dir/pbox_conf.pdf")
PyPlot.savefig("$plot_dir/pbox_conf.png")
PyPlot.close("all")

function get_alpha_level(x, membs)
    index = findlast(x .∈ membs)
    if isnothing(index)
        index = 0
    end
    return index / length(membs)
end

membs = interval.(A.u[1:500], reverse(A.d[501:end]))
xs = range(-6.1, 6.1; length=500)

alphas = get_alpha_level.(xs, Ref(membs))

PyPlot.PyPlot.plot(xs, alphas; color="red")
fill_between(xs, alphas; color="grey", alpha=0.2)

PyPlot.xlabel(L"x"; fontsize=fontsize)
PyPlot.ylabel(L"\alpha"; fontsize=fontsize)
PyPlot.xticks(; fontsize=tick_font)
PyPlot.yticks(; fontsize=tick_font)

fig = gcf()
fig.set_size_inches(6.4, 4.8)

PyPlot.tight_layout()

PyPlot.savefig("$plot_dir/fuzzy_pbox.pdf")
PyPlot.savefig("$plot_dir/fuzzy_pbox.png")
PyPlot.close("all")
