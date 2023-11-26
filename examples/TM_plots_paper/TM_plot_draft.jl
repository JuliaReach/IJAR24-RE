using ReachabilityAnalysis, Plots
using ReachabilityAnalysis.TaylorModels

R = ReachSet(rand(Zonotope), 0 .. 1)
plot(R; vars=(1, 2))

T = overapproximate(R, TaylorModelReachSet)
plot!(T; vars=(1, 2))

# Construction using a sparse polynomial zonotope
# =================================================

# 2x1 + x2 + 2 x1^2 x2
# 2x2 + 2 x1^3 x2
#
G = [2 1 2 0
     0 2 0 2.0]

E = [1 0 2 3
     0 1 1 1]

c = [0, 0.0]

S = SimpleSparsePolynomialZonotope(c, G, E)

plot(S; nsdiv=200, lw=0.0, color=:blue)

# Construction using a Taylor model
# ====================================
x1, x2 = set_variables("x"; order=4, numvars=2)
p1 = 2x1 + x2 + 2x1^2 * x2
p2 = 2x2 + 2x1^3 * x2

orderT = 3
t0 = interval(0) # expansion point
tdom = interval(0) # expansion point
T = [TaylorModel1(Taylor1(pol, orderT), interval(0), t0, tdom) for pol in [p1, p2]]
T = TaylorModelReachSet(T, tdom);

H = overapproximate(T, Hyperrectangle; nsdiv=50)
plot!(H; vars=(1, 2), lw=0.0, c=:red)
#plot!(T, vars=(1, 2), c=:magenta) # just 1 zonotope

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

H = overapproximate(T, Hyperrectangle; nsdiv=200)
plot(H; vars=(1, 2), lw=0.5, c=:blue, alpha=0.03)
plot!(S; nsdiv=50, lw=0.0, color=:blue)
H0 = copy(H)

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
fig = plot(H; vars=(1, 2), lw=1.0, c=:red, alpha=0.05, xlab=L"x_1", ylab=L"x_2")
plot!(fig, H0; vars=(1, 2), lw=0.5, c=:blue, alpha=0.03)
fig

#plot!(S, nsdiv=50, lw=0.0, color=:blue)
