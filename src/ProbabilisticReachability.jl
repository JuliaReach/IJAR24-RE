module ProbabilisticReachability

using Distributions, IntervalArithmetic, PDMats, Primes, ProbabilityBoundsAnalysis, Random,
      LazySets, LinearAlgebra, ReachabilityAnalysis

using Distributions: length
using IntervalArithmetic: Interval
using LazySets: dim
using ProbabilityBoundsAnalysis: AbstractPbox, setSteps
using ReachabilityAnalysis: AbstractContinuousSystem, ReachSolution, TimeInterval, symBox, zeroI

import Base: rand, split
import Distributions: cdf
import LazySets: sample
import MathematicalSystems: initial_state, system
import ProbabilityBoundsAnalysis: AbstractPbox
import ReachabilityAnalysis: post

const IA = IntervalArithmetic
const PBA = ProbabilityBoundsAnalysis

include("copulas.jl")
include("mvpbox.jl")
include("mvncdf.jl")
include("hvolume.jl")
include("failure.jl")
include("api.jl")
include("confidence.jl")

export FailureProbabilityIVP, GauCopCDF, GaussianCopula, Hvolume, IndepCopula, perfect, opposite,
       MvDistribution, U,
       cdf, confidence, confidence, cons_split, failure_probability, get_initial_pbox, get_pbox,
       get_samples, make_cons, make_step_plots,
       nCopula, nPbox, narrow, plot_samples, rand, randomSet, sample, setSteps, split

end  # module
