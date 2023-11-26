using Distributions, PyPlot, LaTeXStrings, IntervalArithmetic

plot_dir = "$(@__DIR__)/plots"
if !isdir(plot_dir)
    mkdir(plot_dir)
end

fontsize = 25
tick_font = 20

A = Normal(0, 1)
xs = range(-3.5, 3.5; length=200)

pdfs = pdf.(A, xs)

alpha1 = 0.05
alpha2 = 0.6

cut1 = quantile(A, alpha1 / 2) .. quantile(A, 1 - alpha1 / 2)
cut2 = quantile(A, alpha2 / 2) .. quantile(A, 1 - alpha2 / 2)

xs1_lo = range(cut1.lo, cut2.lo; length=50)
xs1_hi = range(cut2.hi, cut1.hi; length=50)

xs2 = range(cut2.lo, cut2.hi; length=50)

PyPlot.fill_between(xs1_lo, pdf.(A, xs1_lo); color="red", alpha=0.35)
PyPlot.fill_between(xs1_hi, pdf.(A, xs1_hi); color="red", alpha=0.35)

PyPlot.fill_between(xs2, pdf.(A, xs2); color="blue", alpha=0.5)

PyPlot.plot(xs, pdfs; color="black")
PyPlot.xlabel(L"x"; fontsize=fontsize)
PyPlot.ylabel("pdf"; fontsize=fontsize)
PyPlot.xticks(; fontsize=tick_font)
PyPlot.yticks(range(0, 0.4; length=5); fontsize=tick_font)

xlimm = gca().get_xlim()
ylimm = gca().get_ylim()

PyPlot.tight_layout()

PyPlot.plot([cut2.lo; cut2.hi], [-0.01; -0.01]; color="blue", linewidth=3)
PyPlot.plot([cut1.lo; cut1.hi], [-0.015; -0.015]; color="red", linewidth=3)
PyPlot.xlim(xlimm)
PyPlot.ylim((-0.02, ylimm[2]))

PyPlot.savefig("$plot_dir/density_integration.pdf")
PyPlot.savefig("$plot_dir/density_integration.png")
PyPlot.close("all")

function get_alpha(x, mu, sigma)
    dist = Normal(mu, sigma)
    cdf_val = cdf(dist, x)
    alpha = cdf_val

    if alpha >= 0.5
        alpha = 1 - alpha
    end

    return alpha * 2
end

xs = range(-3.5, 3.5; length=500)

alphas = get_alpha.(xs, 0, 1)

PyPlot.plot(xs, alphas; color="red")
PyPlot.fill_between(xs, alphas; color="grey", alpha=0.2)

xlimm = gca().get_xlim()
ylimm = gca().get_ylim()

xpoint = xlimm[1]
ypoint = ylimm[1]

#=
PyPlot.plot([xpoint; cut1.hi], [0.1; 0.1], "--", color = "black")
PyPlot.plot([cut1.lo; cut1.lo], [0; 0.1], "--", color = "black")
PyPlot.plot([cut1.hi; cut1.hi], [0; 0.1], "--", color = "black")
=#

PyPlot.plot([xpoint; cut2.hi], [alpha2; alpha2], "--"; color="black")
PyPlot.plot([cut2.lo; cut2.lo], [0; alpha2], "--"; color="black")
PyPlot.plot([cut2.hi; cut2.hi], [0; alpha2], "--"; color="black")

PyPlot.plot([cut2.lo; cut2.hi], [-0.01; -0.01]; color="black", linewidth=3)

PyPlot.xlabel(L"x"; fontsize=fontsize)
PyPlot.ylabel(L"\alpha"; fontsize=fontsize)
PyPlot.xticks(; fontsize=tick_font)
PyPlot.yticks(; fontsize=tick_font)
PyPlot.tight_layout()
PyPlot.xlim(xlimm)
PyPlot.ylim((-0.02, ylimm[2]))

PyPlot.savefig("$plot_dir/confidence_precise.pdf")
PyPlot.savefig("$plot_dir/confidence_precise.png")
PyPlot.close("all")
