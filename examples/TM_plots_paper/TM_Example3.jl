using ReachabilityAnalysis
using Plots, Plots.Measures, LaTeXStrings
using ReachabilityAnalysis.TaylorModels

plot_dir = "$(@__DIR__)/plots"
if !isdir(plot_dir)
    mkdir(plot_dir)
end

# Construction using a Taylor model with remainder
# =================================================
x1, x2 = set_variables("x"; order=4, numvars=2)
p1 = 2x1 + x2 + 2x1^2 * x2
p2 = 2x2 + 2x1^3 * x2

orderT = 3
t0 = interval(0) # expansion point
tdom = interval(0) # expansion point
T = [TaylorModel1(Taylor1(pol, orderT), interval(-0.05, 0.05), t0, tdom) for pol in [p1, p2]]
T = TaylorModelReachSet(T, tdom);

H0 = overapproximate(T, Hyperrectangle; nsdiv=200)

# Construction using a Taylor model with non-zero time interval and remainder
# ============================================================================

x1, x2 = set_variables("x"; order=4, numvars=2)
p1 = 2x1 + x2 + 2x1^2 * x2
p2 = 2x2 + 2x1^3 * x2

orderT = 3
t0 = interval(0) # expansion point
tdom = interval(0, 0.2) # expansion point

a = Taylor1(p1, orderT)  # p1(x) + t * p2(x)
a.coeffs[2] = p2

b = Taylor1(p2, orderT) # p2(x) + t * p1(x)
b.coeffs[2] = p1

T = [TaylorModel1(p, interval(-0.05, 0.05), t0, tdom) for p in [a, b]]
T = TaylorModelReachSet(T, tdom);

H = overapproximate(T, Hyperrectangle; nsdiv=300)
fig = plot(; xlab           = L"p_1(x)", ylab           = L"p_2(x)", xtickfontsize  = 14,
           ytickfontsize  = 14,
           xguidefontsize = 20, yguidefontsize = 20, right_margin   = 1.2cm,
           bottom_margin  = 0.2cm)
plot!(fig, H; vars=(1, 2), lw=1.0, c=:red, alpha=0.05)
plot!(fig, H0; vars=(1, 2), lw=0.5, c=:blue, alpha=0.03)
fig

savefig(fig, "$plot_dir/Example3.png")
savefig(fig, "$plot_dir/Example3.pdf")
