# ====================================================
# Convenience functions to work with time intervals
# ====================================================

_check_time(t::AbstractFloat, R) = t ∉ tspan(R)
_check_time(t::IA.Interval, R) = !(t ⊆ tspan(R))

_check_time_intersection(t::AbstractFloat, R) = t ∈ tspan(R)
_check_time_intersection(t::IA.Interval, R) = !isdisjoint(t, tspan(R))
_check_time_intersection(ts::AbstractVector, R) = any(_check_time_intersection(t, R) for t in ts)

# ===============================================
# Method overloads for LazySets
# ===============================================

# conversion from IA.Interval to LazySets.Interval
LazySets.isdisjoint(l::IA.Interval, r) = isdisjoint(LazySets.Interval(l), r)
LazySets.issubset(l::IA.Interval, r) = issubset(LazySets.Interval(l), r)

# conversion from IA.IntervalBox to LazySets.Hyperrectangle
LazySets.isdisjoint(l::IntervalBox, r) = isdisjoint(convert(Hyperrectangle, l), r)
LazySets.issubset(l::IntervalBox, r) = issubset(convert(Hyperrectangle, l), r)
# in 1D case convert to Interval
LazySets.isdisjoint(l::IntervalBox, r::LazySets.Interval) = isdisjoint(LazySets.Interval(l.v[1]), r)
LazySets.issubset(l::IntervalBox, r::LazySets.Interval) = issubset(LazySets.Interval(l.v[1]), r)

# only sufficient checks for high-dimensional unions
LazySets.issubset(X::Hyperrectangle, Y::UnionSet) = issubset(X, Y.X) || issubset(X, Y.Y)
function LazySets.issubset(Z::Zonotope,
                           Y::UnionSet{N,<:LazySets.Interval,<:LazySets.Interval}) where {N}
    return issubset(Z, Y.X) || issubset(Z, Y.Y)
end
function LazySets.issubset(Z::Zonotope, Y::UnionSetArray{N,<:LazySets.Interval}) where {N}
    return any(issubset(Z, X) for X in array(Y))
end

function _to_lazyset(x::Union{LazySet,UnionSet,UnionSetArray})
    return x
end

function _to_lazyset(x::IntervalArithmetic.Interval)
    return LazySets.Interval(x)
end

function _to_lazyset(x::IntervalBox)
    return convert(Hyperrectangle, x)
end

function _to_lazyset(x::AbstractVector)
    return [_to_lazyset(xi) for xi in x]
end

lefts(box::IntervalBox) = getfield.(box, :lo)
rights(box::IntervalBox) = getfield.(box, :hi)

# ===============================================
# Interface with ReachabilityAnalysis.jl API
# ===============================================

function get_initial_pbox(F::Flowpipe; normalized=false)
    return normalized ? F.ext[:initial_pbox_normalized] : F.ext[:initial_pbox]
end
function get_initial_pbox(sol::ReachSolution; normalized=false)
    return get_initial_pbox(sol.F; normalized=normalized)
end

for alg in (:TMJets20, :TMJets21a, :TMJets21b)
    eval(quote
             function post(alg::$alg{N}, ivp::IVP{ST,PBOX}, timespan;
                           Δt0::TimeInterval=zeroI,
                           kwargs...) where {N,ST<:AbstractContinuousSystem,PBOX<:AbstractPbox}
                 return _post(alg, ivp, timespan, Δt0; kwargs...)
             end
         end)
end

function _post(alg, ivp, timespan, Δt0; kwargs...)
    X0range = range(ivp.x0)
    ivp_range = IVP(ivp.s, X0range)
    F = post(alg, ivp_range, timespan; Δt0=Δt0, kwargs...)

    # add p-box information on the extension dictionary
    F.ext[:initial_pbox] = ivp.x0
    if haskey(kwargs, :initial_pbox_normalized)
        F.ext[:initial_pbox_normalized] = kwargs[:initial_pbox_normalized]
    else
        F.ext[:initial_pbox_normalized] = normalize(ivp.x0)
    end
    F.ext[:initial_pbox_range] = X0range
    return F
end

# ===============================================
# Utilities for the examples
# ===============================================

macro modelpath(model_path::String, name::String)
    __source__.file === nothing && return nothing
    _dirname = dirname(String(__source__.file))
    dir = isempty(_dirname) ? pwd() : abspath(_dirname)
    return joinpath(dir, name)
end

# it is assumed that: B ⊆ [-1, 1]^n

###
#   Just a work around due to outward directed rounding in tspan()
###
function narrow(x)
    if isscalar(x)
        return x
    end
    return IntervalArithmetic.Interval(nextfloat(x.lo), prevfloat(x.hi))
end

function ReachabilityAnalysis.evaluate(R::TaylorModelReachSet, t, B::IntervalArithmetic.Interval;
                                       kwargs...)
    return evaluate(R, t, IntervalBox(B); kwargs...)
end

function ReachabilityAnalysis.evaluate(R::TaylorModelReachSet, t=tspan(R),
                                       B::IntervalBox{D,N}=symBox(ReachabilityAnalysis.dim(R));
                                       partition=nothing, nsdiv=1, ntdiv=1) where {D,N}
    @assert B ⊆ symBox(D)
    X = set(R)

    # renormalize dt
    tn = t - tstart(R)
    @assert partition == nothing && nsdiv == 1 && ntdiv == 1

    Xdt = evaluate(X, narrow(tn))

    vec = _evaluate(Xdt, B)  # result is a vector of intervals
    res = length(vec) == 1 ? vec[1] : IntervalBox(vec)

    return res
end

_evaluate(X, B::IntervalBox) = evaluate.(X, Ref(B))
