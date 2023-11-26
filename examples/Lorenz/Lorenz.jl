# # Lorenz

module Lorenz #jl

using ReachabilityAnalysis,
      ProbabilisticReachability,
      ProbabilityBoundsAnalysis
using Test
const PBA = ProbabilityBoundsAnalysis #hide
const IA = ProbabilityBoundsAnalysis.IntervalArithmetic

MODELs = "Lorenz"

include("Lorenz_model.jl")

include("Lorenz_failure.jl")

@time pbs = get_pbox(sol, 18.25)

end
