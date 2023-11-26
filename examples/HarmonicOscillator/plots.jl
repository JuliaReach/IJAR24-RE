using ProbabilityBoundsAnalysis
const PBA = ProbabilityBoundsAnalysis
using ProbabilisticReachability, LazySets
using Plots, LaTeXStrings
using Plots: plot

MODELs = "HarmonicOscillator"

plot_dir = "$(@__DIR__)/plots"
if !isdir(plot_dir)
    mkdir(plot_dir)
end

include("$(MODELs).jl")
