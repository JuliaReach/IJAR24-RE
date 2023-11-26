using BivariateCopulas, PyPlot

plot_dir = "$(@__DIR__)/plots"
if !isdir(plot_dir)
    mkdir(plot_dir)
end

C1 = Gaussian(0.8)

fontsize = 24
font2 = 30

plotContourCdf(C1)
PyPlot.xlabel("u"; fontsize=font2)
PyPlot.ylabel("v"; fontsize=font2)

PyPlot.xticks(; fontsize=fontsize);
yticks(; fontsize=fontsize);

PyPlot.savefig("$plot_dir/GauCop1.png")
