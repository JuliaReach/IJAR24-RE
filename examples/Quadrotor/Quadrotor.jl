# # Quadrotor

module Quadrotor

using ReachabilityAnalysis, ProbabilisticReachability, ProbabilityBoundsAnalysis
const PBA = ProbabilityBoundsAnalysis
const IA = ProbabilityBoundsAnalysis.IntervalArithmetic

using ProbabilisticReachability: normalize
using ProbabilityBoundsAnalysis: setSteps
using ReachabilityAnalysis: is_intersection_empty

# ------------------------------------------------------------------------------
# Algorithm options

import Main.BENCHMARK_MODE
if !@isdefined(BENCHMARK_MODE)
    global BENCHMARK_MODE = true  # run full benchmarks
end

println("=========")
println("Quadrotor")
println("=========")

SAMESTEPS = false
NSTEPS = BENCHMARK_MODE ? 100 : 2

include("Quadrotor_model.jl")
include("Quadrotor_failure.jl")

end
