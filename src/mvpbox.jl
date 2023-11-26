######
# This file is part of the ProbabilisticReachability package
#
#   Definition of multivariate p-box and related functions
#
######

##
#   A multivariate p-box defined in terms of N pbox marginals and an N copula
##
struct nPbox{N} <: AbstractPbox
    marginals::Vector{pbox}
    copula::nCopula{N}
end

dims(X::nPbox{N}) where {N} = N
const PBOX = Union{pbox,nPbox}

##
#  Constructing a mv-pbox from a copulas using Sklar's Theorem
##
function (obj::nCopula)(marginals::Vector{pbox})
    if dims(obj) != length(marginals)
        throw(ArgumentError("Dimension of copula and number of marginals must coincide"))
    end

    return nPbox{dims(obj)}(marginals, obj)
end

# extend dimension methods
LazySets.dim(X::nPbox{N}) where {N} = N
ReachabilityAnalysis._dim(X::nPbox{N}) where {N} = N

###
#   Multivariate distribution type
###

const dist_type = Distribution{Univariate,Continuous}

struct MvDistribution{N,D<:dist_type}
    marginals::Vector{D}
    copula::nCopula{N}
end

function (obj::nCopula)(marginals::Vector{<:dist_type})
    if dims(obj) != length(marginals)
        throw(ArgumentError("Dimension of copula and number of marginals must coincide"))
    end

    return MvDistribution(marginals, obj)
end

function cdf(X::nPbox, u::Vector{<:Real})
    if dims(X) != length(u)
        throw(ArgumentError("Dimension of p-box and vector must coincide"))
    end

    probs = interval.(cdf.(X.marginals, u))  # Get marginal cdfs

    p_lb = X.copula(getfield.(probs, :lo))  # Lower bound on mv cdf
    p_ub = X.copula(getfield.(probs, :hi))  # Upper bound on mv cdf

    return interval(p_lb, p_ub)
end

###
#   Sampling function for distributions, mv-distributions, p-box and mv-pboxes
###

LazySets._default_sampler(X::pbox) = nothing
LazySets._default_sampler(X::nPbox) = nothing

function sample(X::nPbox, num_samples::Int; kwargs...)

    ## perhaps convert back to a continuous distribution

    distribution = get(kwargs, :sampler, X)
    return rand_me(distribution, num_samples)
end

function sample(X::pbox, num_samples::Int; kwargs...)

    ## perhaps convert back to a continuous distribution

    distribution = get(kwargs, :sampler, X)
    return rand_me(distribution, num_samples)
end

rand_me(X::pbox, num_samples) = [[mid(x)] for x in rand(X, num_samples)]
rand_me(X, num_samples) = [[x] for x in rand(X, num_samples)]

function rand_me(X::nPbox, num_samples)
    fun = (x, y) -> mid(cut(x, y))
    return _rand_me(X, num_samples, fun)
end

function rand_me(X::MvDistribution, num_samples)
    fun = (x, y) -> quantile(x, y)
    return _rand_me(X, num_samples, fun)
end

function _rand_me(X, num_samples, fun)
    C = X.copula
    margs = X.marginals

    CopSamps = ProbabilisticReachability.rand(C, num_samples)

    ss = [CopSamps[i, :] for i in 1:size(CopSamps, 1)]

    return [[fun(margs[i], s[i]) for i in 1:size(CopSamps, 2)] for s in ss]
end

###
#   Get sample ensemble at a specific time
###

function get_samples(sol::ReachabilityAnalysis.ReachSolution, t, dims)
    samps = [reduce(vcat, sim(t)[d] for sim in ensemble(sol)) for d in dims]
    return reduce(hcat, samps)
end

function plot_samples(s::Vector{Float64}, col="red", alpha=1)
    s = sort(s)
    i = range(0, 1; length=length(s))
    return PyPlot.step(s, i; color=col, alpha=alpha)
end

##
#   A multivariate random set. Collection of interval boxes (Focal elements) with probability masses correlated by a copula
##
struct randomSet{N,T}
    fe::Vector{IntervalBox{N,T}}      # Focal elements
    masses::Vector{Float64}
end

###
#  Turns an nPbox to a randomSet
#  Inspired by IntervalArithmetic.mince. Does H-volume (HvolumeFast) calculation in one go.
###
@inline function split(X::nPbox{N}) where {N}
    Cop = X.copula

    if Cop.fam == "π"
        return splitIndependence(X)
    end

    marginalInts = [interval.(x.u, x.d) for x in X.marginals]   # Marginal Focal Elements

    NumEls = length.(marginalInts)

    focalElements = Vector{IntervalBox{N,Float64}}(undef, prod(NumEls))

    for (k, cut_indices) in enumerate(CartesianIndices(Tuple(NumEls)))
        focalElements[k] = IntervalBox([marginalInts[i][cut_indices[i]] for i in 1:N])  # Creates joint Focal Elements
    end

    NumEls = NumEls .+ 1
    mCdfs = [range(0, 1; length=n) for n in NumEls]

    jCdfs = zeros(Tuple(NumEls))
    for (k, cut_indices) in enumerate(CartesianIndices(Tuple(NumEls)))
        jCdfs[k] = Cop([mCdfs[i][cut_indices[i]] for i in eachindex(NumEls)])
    end

    masses = HvolumeFast(jCdfs)

    return focalElements, masses[:]
end

###
#  Turns an nPbox to a randomSet. Special case where random set is independent
###
@inline function splitIndependence(X::nPbox{N}) where {N}
    marginalInts = [interval.(x.u, x.d) for x in X.marginals]   # Marginal Focal Elements

    NumEls = length.(marginalInts)

    marginalMasses = [fill(1 / num, num) for num in NumEls]

    focalElements = Vector{IntervalBox{N,Float64}}(undef, prod(NumEls))
    masses = Vector{Float64}(undef, prod(NumEls))

    for (k, cut_indices) in enumerate(CartesianIndices(Tuple(NumEls)))
        focalElements[k] = IntervalBox([marginalInts[i][cut_indices[i]] for i in 1:N])  # Creates joint Focal Elements
        masses[k] = prod([marginalMasses[i][cut_indices[i]] for i in 1:N])             # Product of marginal masses
    end

    return focalElements, masses
end

###
#   Split, but requires more calls to copula
###
function splitOld(X::nPbox{N}) where {N}
    if X.copula.fam == "π"
        return splitIndependence(X)
    end

    marginalInts = [interval.(x.u, x.d) for x in X.marginals]   # Marginal Focal Elements

    Cop = X.copula

    NumEls = length.(marginalInts)
    marginalCdfs = [mince(interval(0, 1), num) for num in NumEls]    # Marginal cdf partition

    focalElements = Vector{IntervalBox{N,Float64}}(undef, prod(NumEls))
    masses = Vector{Float64}(undef, prod(NumEls))

    for (k, cut_indices) in enumerate(CartesianIndices(Tuple(NumEls)))
        focalElements[k] = IntervalBox([marginalInts[i][cut_indices[i]] for i in 1:N])  # Creates joint Focal Elements
        thisCdfBox = IntervalBox([marginalCdfs[i][cut_indices[i]] for i in 1:N])        # Interval box of cdf partition
        masses[k] = Hvolume(Cop, thisCdfBox)                                           # Masses are the h-volume of the copula in partition
    end

    return focalElements, masses
end

function split(X::pbox)
    FocalElements = interval.(X.u, X.d)
    num = length(FocalElements)
    masses = fill(1 / num, num)

    return FocalElements, masses
end

###
#   Normalise the multivariate pbox to [-1, 1]
###
function normalize(x::nPbox{N}) where {N}
    cop = x.copula
    margs = [pbox() for i in 1:N]

    for i in 1:N
        rangelo = x.marginals[i].u[1]
        rangehi = x.marginals[i].d[end]

        us = (x.marginals[i].u .- rangelo) ./ (rangehi - rangelo) .* 2 .- 1
        ds = (x.marginals[i].d .- rangelo) ./ (rangehi - rangelo) .* 2 .- 1

        margs[i] = pbox(us, ds; shape=x.marginals[i].shape)
    end
    return cop(margs)
end

# one-dimensional case
function normalize(x::pbox)
    rangelo = x.u[1]
    rangehi = x.d[end]

    us = (x.u .- rangelo) ./ (rangehi - rangelo) .* 2 .- 1
    ds = (x.d .- rangelo) ./ (rangehi - rangelo) .* 2 .- 1

    return pbox(us, ds; shape=x.shape)
end

###
#   Range of a multivariate p-box
###

function Base.range(x::nPbox{N}) where {N}
    return IntervalBox(ProbabilityBoundsAnalysis.range.(x.marginals))
end

###
#   For plotting p-box transitions in time
###

subscript(i) = join(Char(0x2080 + d) for d in reverse!(digits(i)))

function make_step_plots(sol, Npboxes; dims=1:length(sol), fontsize=24, alpha=0.4, figure=nothing,
                         poly3d=nothing, xlabel=nothing, ylabel=nothing, zlabel=nothing,
                         xlim=nothing, ylim=nothing)
    t_span = tspan(sol)
    ts = range(t_span.lo, t_span.hi; length=Npboxes)

    pboxes = [get_pbox(sol, t) for t in ts]

    iis = range(0, 1; length=length(pboxes[1][1].u))
    jjs = reverse(iis)

    fig = figure(; figsize=(10, 10))
    for i in dims
        ax = fig.add_subplot(111; projection="3d")

        maxPoint = pboxes[1][i].d[end]
        minPoint = pboxes[1][i].u[1]

        for j in 1:Npboxes
            maxPoint = max(pboxes[j][i].d[end], maxPoint)
            minPoint = min(pboxes[j][i].u[1], minPoint)

            t = ts[j]
            y = pboxes[j][i].u
            z = pboxes[j][i].d
            verts1 = [(y[i], iis[i]) for i in 1:length(iis)]
            z1 = reverse(z)
            verts2 = [(z1[i], jjs[i]) for i in 1:length(jjs)]

            verts = vcat(verts1, verts2)

            ax.add_collection3d(poly3d([verts]; facecolor="blue", alpha=alpha, linestyle="-",
                                       edgecolor="black", linewidths=2); zs=[t], zdir="x")
        end

        xlabel("t"; fontsize=fontsize)
        ylabel("x" * subscript(i) * "(t)"; fontsize=fontsize)
        zlabel("cdf"; fontsize=fontsize)

        xlim([t_span.lo, t_span.hi])
        ylim([minPoint, maxPoint])
    end
    return fig
end
