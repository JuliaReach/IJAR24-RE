global BENCHMARK_MODE
if !@isdefined(BENCHMARK_MODE)
    BENCHMARK_MODE = true  # run full benchmarks
end

if BENCHMARK_MODE
    @info "Running long benchmarks..."
else
    @info "Running fast benchmarks..."
end

PLOT_FILES = ["UnivariateOscillator/plots.jl",
              "HarmonicOscillator/plots.jl",
              "Lorenz/plots.jl",
              "Quadrotor/plots.jl",
              "Duffing/Duffing.jl",
              "TM_plots_paper/TM_Example2.jl",
              "TM_plots_paper/TM_Example3.jl",
              # MV pbox plots
              "pbox_plots_paper/plot_MV_pbox.jl",
              "pbox_plots_paper/split_mv_pbox.jl",
              "pbox_plots_paper/splitting.jl",
              # plots for consonant confidence intervals
              "confidence_plots_paper/density_int.jl",
              "confidence_plots_paper/pbox_conf.jl",
              "confidence_plots_paper/consonant_slice.jl",
              "confidence_plots_paper/consonant_plot.jl"]

foreach(PLOT_FILES) do file
    @info "Running $file..."
    include(file)
    return nothing
end
