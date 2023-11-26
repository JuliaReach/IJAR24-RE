# # Univariate oscillator

module UnivariateOscillator #jl

# ## Model

# We consider a one-dimensional nonlinear ODE with analytic solution.
# The differential equation is:
#
# ```math
# x'(t) = -x(t) ~ \sin(t)
# ```
# and the solution is:
#
# ```math
# x(t) = x_0 e^{\cos(t) - 1},\qquad t \geq 0.
# ```

using ReachabilityAnalysis,
      ProbabilisticReachability,
      ProbabilityBoundsAnalysis
const PBA = ProbabilityBoundsAnalysis #hide
const IA = PBA.IntervalArithmetic

using ProbabilisticReachability: _failure_probability_analytic

import Main.BENCHMARK_MODE
if !@isdefined(BENCHMARK_MODE)
    global BENCHMARK_MODE = true  # run full benchmarks
end

println("=====================")
println("Univariate oscillator")
println("=====================")

# The corresponding Julia function describing the ODE is:

@taylorize function f!(dx, x, p, t)
    return dx[1] = -x[1] * sin(t)
end

# ## Probabilistic initial-value problem

# We assume that the initial condition belongs to the p-box $x(0) \in U_0([-1, 0], [0, 1])$.
# Such p-box describes a uniform distribution whose end-points are intervals.

## fix the discretization steps
ProbabilityBoundsAnalysis.setSteps(BENCHMARK_MODE ? 100 : 10)

U0 = uniform(-1, 1)

## initial-value problem
prob = @ivp(x' = f!(x), x(0) ∈ U0, dim = 1)

# We can plot the cdf of the initial distribution.
#=
PBA.plot(U0, save=true, name="univariate_oscillator_U0", plotting=false) #hide
PBA.plot(U0, plotting=false)
PBA.gcf() #hide
=#
# ## p-box propagation

# We obtain a flowpipe corresponding to the probabilistic initial-value problem
# defined above from time $t = 0$ to $T = 15$.

alg = TMJets21a(; orderT=10, orderQ=2, abstol=1e-15)
Tmax = 15

sol = solve(prob; alg=alg, T=Tmax)

## number of computed reach-sets
length(sol)

# Data for Example 4 which illustrates the evaluation of a Taylor model over a time domain
# and a spatial domain.
X = sol[1]
println("First time interval: $(tspan(sol[1]))")
println("Taylor model: $(set(X))")

Δt = interval(0, 0.01)
println("Evaluation in the time interval $Δt gives $(evaluate(X, Δt))")
Δd = -1 .. 1
R1 = ReachabilityAnalysis._overapproximate(X, Hyperrectangle; Δt, dom=IntervalBox(Δd))
println("Evaluation in the time interval $Δt and spatial domain $Δd is $(convert(LazySets.Interval, set(R1)))")

Δt = interval(0.04, 0.05)
println("Evaluation in the time interval $Δt gives $(evaluate(X, Δt))")
Δd = -1 .. 1
R1 = ReachabilityAnalysis._overapproximate(X, Hyperrectangle; Δt, dom=IntervalBox(Δd))
println("Evaluation in the time interval $Δt and spatial domain $Δd is $(convert(LazySets.Interval, set(R1)))")

# It is illustrative to plot the computed flowpipe and the known analytic solution
# for a range of initial conditions that cover the range of `U0`.

analytic_sol(x0) = t -> x0 * exp(cos(t) - 1.0)

dt = range(0, Tmax; length=100)
x0vals = range(-1, 1; length=25)

#=
fig = plot() #!jl
plot!(fig, sol, vars=(0, 1), lw=0.) #!jl
[plot!(fig, dt, analytic_sol(x0).(dt), lab="", c=:magenta, lw=2.0, xlab="t", ylab="x(t)") for x0 in x0vals]; #!jl
fig  #!jl
=#

# ## Computation of a failure probability

# Here we approximate the probability that the reach-set at time $t = 4*π ≈ 12.57$ enters into the interval $[-1, -0.5]$.

failure_domain = interval(-1, -0.5)
failure_time_domain = 4 * π

# extended failure domain
failure_domain2 = interval(0, 0.5)
failure_domain_both = UnionSet(LazySets.Interval(failure_domain),
                               LazySets.Interval(failure_domain2))

p_lower_precise = failure_probability(sol, failure_time_domain, failure_domain, U0)

# The analytic probability in this case coincides with the probability in the initial set.

X0 = U0.u[1] .. U0.d[end]
p_lower_analytic_uniform = diam(X0 ∩ failure_domain) / diam(X0)

# It can be verified that such interval is slightly more accurate than
# evaluating the cumulative distribution function of the propagated p-box.

pb = get_pbox(sol, failure_time_domain)

#-

cdf(pb, -0.1)

#-
#=
PBA.plot(pb, save=true, name="univariate_oscillator_t15", plotting=false) #hide
PBA.plot(pb, plotting=false)
PBA.gcf() #hide
=#

# ## Check pbox enclosure
# Check that the Taylor solution is an outer robust approximation of the analytic solution

#=
import Base.⊆

function issubseteq(x:: pbox, y :: pbox)

    xInts = interval.(x.u, x.d)
    yInts = interval.(y.u, y.d)

    return all( xInts .⊆ yInts)

end

⊆(x ::pbox, y :: pbox) = issubseteq(x,y)

Ntimes = 20

#ts = 0.1:Ntimes:14.9
ts = range(0, Tmax, length= Ntimes)

pbs1 = get_pbox.(Ref(sol), ts)

anal_pbs = analytic_sol(U0).(ts)

all_in = all( anal_pbs .⊆ pbs1)
=#

#=
for i =1:Ntimes

    PBA.plot(anal_pbs[i], name ="same_$i", col = "red")
    PBA.plot(pbs1[i], name ="same_$i", col = "blue")
    PBA.gcf() #hide
end
=#

#[PBA.plot(anal_pbs[i], name ="same_$i", col = "red") for i =1:Ntimes]
#[PBA.plot(pbs1[i], name ="same_$i", col = "blue") for i =1:Ntimes]
# ## Gaussian distribution

# Here we consider using a non-uniform distribution over the set of initial states.
# Consider a Gaussian distribution of mean $0$ and and standard deviation

# We create a normal distribution truncated to the range -1 .. 1.

#nn = normal(0, 0.5)
#ii = makepbox(-1 .. 1)
#imp([nn, ii])

# =====

# modify discretization parameter
if BENCHMARK_MODE
    steps_list = range(10; step=20, stop=1000)
else
    steps_list = range(10; step=20, stop=30)
end

p_lower_precise_uniform_tune = Vector{IA.Interval}(undef, length(steps_list))
@inbounds for (i, steps) in enumerate(steps_list)
    ProbabilityBoundsAnalysis.setSteps(steps)
    U0 = uniform(-1, 1)
    p_lower_precise_uniform_tune[i] = failure_probability(sol, failure_time_domain,
                                                          failure_domain_both, U0)
end

# reset steps
ProbabilityBoundsAnalysis.setSteps(BENCHMARK_MODE ? 100 : 10)
U0 = uniform(-1, 1)

# =====

# imprecise distribution
U0_imp = uniform(interval(-1, 0), interval(0, 1))
p_lower_imprecise = failure_probability(sol, failure_time_domain, failure_domain, U0_imp)

pb_imp = get_pbox(sol, failure_time_domain, U0_imp)

# =====

# simplification: skip intersection of failure domain with X0 because it is the identity here
p_lower_analytic_uniform_larger_domain = volume(failure_domain_both) / diam(X0)

p_lower_precise_larger_domain = failure_probability(sol, failure_time_domain, failure_domain_both,
                                                    U0)

# =====

# Define confidence levels
confidence_levels = BENCHMARK_MODE ? [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.001] :
                    [1.0, 0.5]

# confidence flowpipes/percentile sets
conf = confidence(sol, confidence_levels; cons=false)
conf_imp = confidence(sol, confidence_levels; pbinit=U0_imp, cons=false)

# using consonsance
conf_cons = confidence(sol, confidence_levels; cons=true)
conf_imp_cons = confidence(sol, confidence_levels; pbinit=U0_imp, cons=true)

end #jl
