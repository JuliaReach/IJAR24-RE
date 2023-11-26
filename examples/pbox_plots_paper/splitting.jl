using ProbabilisticReachability, ProbabilityBoundsAnalysis, PyPlot, LaTeXStrings, Distributions,
      IntervalArithmetic

function plot_splits_1D(x; subpl=missing, fontsize=18)
    NumEl = length(x)
    #is = mince(0..1, NumEl)
    is = range(0, 1; length=NumEl + 1)
    cdfs = interval.(is[1:(end - 1)], is[2:end])

    plotBoxes(x, cdfs, subpl; alpha=0.2, linewidth=1.5)

    PyPlot.xticks(; fontsize=fontsize)
    PyPlot.yticks(; fontsize=fontsize)
    PyPlot.xlabel("Dist range"; fontsize=fontsize)
    return PyPlot.ylabel("CDF"; fontsize=fontsize)
end

plot_dir = "$(@__DIR__)/plots"
if !isdir(plot_dir)
    mkdir(plot_dir)
end

PyPlot.using3D()

ft_big = 55
ft = 44
ft_ticks = 40

fig_size = (10, 10)

x1_param_1 = interval(1, 2)
x1_param_2 = 3

x2_param_1 = interval(5, 6)
x2_param_2 = 2

x3_param_1 = interval(0, 0.5)
x3_param_2 = interval(2, 3)

ProbabilityBoundsAnalysis.setSteps(5)
x1 = beta(x1_param_1, x1_param_2)
ProbabilityBoundsAnalysis.setSteps(3)
x2 = gamma(x2_param_1, x2_param_2)
ProbabilityBoundsAnalysis.setSteps(10)
x3 = normal(x3_param_1, x3_param_2)

fe_x1, _ = split(x1)
fe_x2, _ = split(x2)
fe_x3, _ = split(x3)

###
#   Plotting without showing exact p-box
###
fig = figure("x1"; figsize=fig_size)
subpl = fig.add_subplot()
plot_splits_1D(fe_x1; subpl=subpl, fontsize=ft)
PyPlot.xlabel(L"$X_1$"; fontsize=ft)
PyPlot.title(L"$F_{X_{1}}$"; fontsize=ft_big)
PyPlot.tight_layout()

PyPlot.savefig("$plot_dir/x1_split.pdf")
PyPlot.savefig("$plot_dir/x1_split.png")
PyPlot.close("all")

fig = figure("x2"; figsize=fig_size)
subpl = fig.add_subplot()
plot_splits_1D(fe_x2; subpl=subpl, fontsize=ft)
PyPlot.xlabel(L"$X_2$"; fontsize=ft)
PyPlot.title(L"$F_{X_{2}}$"; fontsize=ft_big)
PyPlot.tight_layout()

PyPlot.savefig("$plot_dir/x2_split.pdf")
PyPlot.savefig("$plot_dir/x2_split.png")
PyPlot.close("all")

fig = figure("x3"; figsize=fig_size)
subpl = fig.add_subplot()
plot_splits_1D(fe_x3; subpl=subpl, fontsize=ft)
PyPlot.xticks([-10, -5, 0, 5, 10])
PyPlot.xlabel(L"$X_3$"; fontsize=ft)
PyPlot.title(L"$F_{X_{3}}$"; fontsize=ft_big)
PyPlot.tight_layout()

PyPlot.savefig("$plot_dir/x3_split.pdf")
PyPlot.savefig("$plot_dir/x3_split.png")
PyPlot.close("all")

###
#   Plotting with exact p-box
###
N_points = 1000

# X1
x1_lo = minimum(fe_x1).lo
x1_hi = maximum(fe_x1).hi
x1_range = range(x1_lo, x1_hi; length=N_points)

x1_cdf_1 = cdf.(Beta(inf(x1_param_1), inf(x1_param_2)), x1_range)
x1_cdf_2 = cdf.(Beta(inf(x1_param_1), sup(x1_param_2)), x1_range)
x1_cdf_3 = cdf.(Beta(sup(x1_param_1), inf(x1_param_2)), x1_range)
x1_cdf_4 = cdf.(Beta(sup(x1_param_1), sup(x1_param_2)), x1_range)

x1_cdf_lo = minimum([x1_cdf_1 x1_cdf_2 x1_cdf_3 x1_cdf_4]; dims=2)
x1_cdf_hi = maximum([x1_cdf_1 x1_cdf_2 x1_cdf_3 x1_cdf_4]; dims=2)

fig = figure("x1"; figsize=fig_size)
subpl = fig.add_subplot()
plot_splits_1D(fe_x1; subpl=subpl, fontsize=ft)
PyPlot.plot(x1_range, x1_cdf_lo, "blue")
PyPlot.plot(x1_range, x1_cdf_hi, "blue")

PyPlot.xlabel(L"$X_1$"; fontsize=ft)
PyPlot.title(L"$F_{X_{1}}$"; fontsize=ft_big)
PyPlot.tight_layout()

PyPlot.savefig("$plot_dir/x1_split_approx.pdf")
PyPlot.savefig("$plot_dir/x1_split_approx.png")
PyPlot.close("all")

# X2
x2_lo = minimum(fe_x2).lo
x2_hi = maximum(fe_x2).hi
x2_range = range(x2_lo, x2_hi; length=N_points)

x2_cdf_1 = cdf.(Gamma(inf(x2_param_1), inf(x2_param_2)), x2_range)
x2_cdf_2 = cdf.(Gamma(inf(x2_param_1), sup(x2_param_2)), x2_range)
x2_cdf_3 = cdf.(Gamma(sup(x2_param_1), inf(x2_param_2)), x2_range)
x2_cdf_4 = cdf.(Gamma(sup(x2_param_1), sup(x2_param_2)), x2_range)

x2_cdf_lo = minimum([x2_cdf_1 x2_cdf_2 x2_cdf_3 x2_cdf_4]; dims=2)
x2_cdf_hi = maximum([x2_cdf_1 x2_cdf_2 x2_cdf_3 x2_cdf_4]; dims=2)

fig = figure("x2"; figsize=fig_size)
subpl = fig.add_subplot()
plot_splits_1D(fe_x2; subpl=subpl, fontsize=ft)
PyPlot.plot(x2_range, x2_cdf_lo, "blue")
PyPlot.plot(x2_range, x2_cdf_hi, "blue")

PyPlot.xlabel(L"$X_2$"; fontsize=ft)
PyPlot.title(L"$F_{X_{2}}$"; fontsize=ft_big)
PyPlot.tight_layout()

PyPlot.savefig("$plot_dir/x2_split_approx.pdf")
PyPlot.savefig("$plot_dir/x2_split_approx.png")
PyPlot.close("all")

# X3
x3_lo = minimum(fe_x3).lo
x3_hi = maximum(fe_x3).hi
x3_range = range(x3_lo, x3_hi; length=N_points)

x3_cdf_1 = cdf.(Normal(inf(x3_param_1), inf(x3_param_2)), x3_range)
x3_cdf_2 = cdf.(Normal(inf(x3_param_1), sup(x3_param_2)), x3_range)
x3_cdf_3 = cdf.(Normal(sup(x3_param_1), inf(x3_param_2)), x3_range)
x3_cdf_4 = cdf.(Normal(sup(x3_param_1), sup(x3_param_2)), x3_range)

x3_cdf_lo = minimum([x3_cdf_1 x3_cdf_2 x3_cdf_3 x3_cdf_4]; dims=2)
x3_cdf_hi = maximum([x3_cdf_1 x3_cdf_2 x3_cdf_3 x3_cdf_4]; dims=2)

fig = figure("x3"; figsize=fig_size)
subpl = fig.add_subplot()
plot_splits_1D(fe_x3; subpl=subpl, fontsize=ft)
PyPlot.plot(x3_range, x3_cdf_lo, "blue")
PyPlot.plot(x3_range, x3_cdf_hi, "blue")

PyPlot.xticks([-10, -5, 0, 5, 10])
PyPlot.xlabel(L"$X_3$"; fontsize=ft)
PyPlot.title(L"$F_{X_{3}}$"; fontsize=ft_big)
PyPlot.tight_layout()

PyPlot.savefig("$plot_dir/x3_split_approx.pdf")
PyPlot.savefig("$plot_dir/x3_split_approx.png")
PyPlot.close("all")

#=
C1 = GaussianCopula([1 0.8; 0.8 1])
F12 = C1([x2,x3])

fe, masses = split(F12)
boxes = [Vector(f.v) for f in fe]

bbs = Matrix{Interval{Float64}}(undef,length(boxes),2)
for i = 1:length(boxes)
    bbs[i,1] = boxes[i][1]
    bbs[i,2] = boxes[i][2]
end

plotBoxes(bbs[:,1], bbs[:,2], alpha= 0.1)
=#
