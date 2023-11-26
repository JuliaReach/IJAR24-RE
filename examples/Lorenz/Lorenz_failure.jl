using ProbabilisticReachability, ProbabilityBoundsAnalysis

import Main.BENCHMARK_MODE
if !@isdefined(BENCHMARK_MODE)
    global BENCHMARK_MODE = true  # run full benchmarks
end

println("======")
println("Lorenz")
println("======")

ProbabilityBoundsAnalysis.setSteps(BENCHMARK_MODE ? 100 : 4)

x1 = U(interval(0.9, 1), 1.1)
x2 = U(-1e-2, interval(0, 1e-2))
x3 = U(-1e-2, 1e-2)

C = IndepCopula(3)

Xpbox = C([x1, x2, x3])

prob = @ivp(x' = lorenz!(x), dim = 3, x(0) ∈ Xpbox)

alg = TMJets21a(; abstol=1e-15, orderT=10, orderQ=2, maxsteps=50_000)

sol = solve(prob; T=20.0, alg=alg)

# p = get_pbox(sol, 18.25)
#
# # To see the p-box:
# # ProbabilityBoundsAnalysis.plot(p[2])
#
# # The range of the p-box is large: [-5.1749, 506.43]
# # However, the p-box still gives good information within this interval
#
# # For example, the proportion of the trajectories which are negative is [0.199, 0.2]
# mass(p[2], interval(-Inf, 0))
#
# # The portion of the trajectories in [100, 200] is in [0.2, 0.4]
# mass(p[2], interval(100, 200))
#
# # The portion of the trajectories that are larger than 300 is in [0.09999, 0.2]
# mass(p[2], interval(300, Inf))

# confidence plot
confidences = BENCHMARK_MODE ? [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.001] : [1.0]
c = confidence(sol, confidences)
