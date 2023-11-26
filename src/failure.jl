######
# This file is part of the ProbabilisticReachability package
#
#   Defines failure problems and solution methods.
#
######

# ============================================================
# Structures to represent probabilistic reachability problems
# ============================================================

# abstract supertype for all verification problem types
abstract type AbstractVerificationProblem end

struct ProbabilisticIVP{T,D,ST<:nPbox{D},PT<:InitialValueProblem{T,ST},US<:LazySet} <:
       AbstractVerificationProblem
    ivp::PT
    unsafe_state::US
end

system(prob::ProbabilisticIVP) = prob.ivp.s
initial_state(prob::ProbabilisticIVP) = prob.ivp.x0
unsafe_state(prob::ProbabilisticIVP) = prob.unsafe_state

# if the unsafe states are not specified => the set complement of the invariant set is used
# (but only an upper bound is obtained)
function ProbabilisticIVP(ivp::PT) where {ST,D,XT<:nPbox{D},PT<:InitialValueProblem{ST,XT}}
    I = stateset(ivp.s)
    Iᶜ = complement(I)
    return ProbabilisticIVP(ivp, Iᶜ)
end

# ===========================================
# Parallel algorithm for failure probability
# ===========================================

# compute each flowpipe in parallel, then analyze the resulting prob distribution
struct ParallelSolver
    #
end

struct FailureProbabilitySolution{IT,FT,MT,ST,AT}
    bounds::IT              # interval bounds on the failure probability
    flowpipes::FT           # flowpipe(s) computed
    masses::MT              # Basic mass asignments of the flowpipes
    simulations::ST         # ensemble simulation (`nothing` if it wasn't computed)
    alg::AT                 # set propagation method
end

function failure_probability(ivp::InitialValueProblem,
                             unsafe_states::LazySet;
                             solver=ParallelSolver(),
                             kwargs...)
    prob = ProbabilisticIVP(ivp, unsafe_states)
    return failure_probability(solver, prob; kwargs...)
end

# Compute each flowpipe in parallel, then find the upper and lower bounds on the
# probability distribution by disjointness and inclusion checks
function failure_probability(::ParallelSolver, prob::ProbabilisticIVP; kwargs...)
    S = system(prob)
    X0p = initial_state(prob)
    U = unsafe_state(prob)

    # transform initial p-box into a random set
    focalElements, masses = split(X0p)

    # parallel flowpipe computation for each focal element
    ivp_FE = IVP(S, focalElements)
    sol = solve(ivp_FE; kwargs...)

    # compute upper and lower bounds
    num_masses = length(masses)
    @assert length(sol) == num_masses

    # inclusion check
    is_in_U(X) = X ⊆ U

    upperPB = 0.0
    lowerPB = 0.0
    @inbounds for i in 1:num_masses
        # does the flowpipe intersect the unsafe states?
        hit = !isdisjoint(sol[i], U)

        if hit
            upperPB += masses[i]

            # is there a reach-set completely inside the unsafe states?
            if any(is_in_U, sol[i])
                lowerPB += masses[i]
            end
        end
    end
    bounds = interval(lowerPB, upperPB)

    return FailureProbabilitySolution(bounds, flowpipe(sol), masses, ensemble(sol), sol.alg)
end

# ========================================================================
# Failure probability for Taylor model reach-sets: one-dimensional case
# ========================================================================

# failure probability at a time point, one-dimensional case
function failure_probability(R::TaylorModelReachSet, t, U, pb::pbox; normalize=true, cons=false)
    U = _to_lazyset(U)

    # check sufficient conditions
    if _check_time(t, R) || isdisjoint(box_approximation(R), U)
        return interval(0)
    end

    # transform initial p-box into a random set
    pbn = normalize ? ProbabilisticReachability.normalize(pb) : pb

    if cons
        focal_elems, masses = cons_split(pbn)
    else
        focal_elems, masses = split(pbn)
    end

    return _failure_probability(R, t, U, focal_elems, masses)
end

# ===============================================================
# Failure probability for Taylor model reach-sets: common code
# ===============================================================

function _failure_probability(R::TaylorModelReachSet, t, U, focal_elems, masses)
    upperPB = 0.0
    lowerPB = 0.0

    for i in 1:length(focal_elems)
        Bi = focal_elems[i]
        x = evaluate(R, t, Bi)
        x = _to_lazyset(x)

        # does the reach-set intersect the unsafe states?
        if !isdisjoint(x, U)
            upperPB += masses[i]

            # is the interval included in the unsafe states?
            if x ⊆ U
                lowerPB += masses[i]
            end
        end
    end
    bounds = interval(lowerPB, upperPB)

    return bounds
end

# ======================================================================
# Failure probability for Taylor model solutions
# ======================================================================

# no p-box -> use the distribution in the initial set
function failure_probability(sol::ReachSolution, t::N, U; kwargs...) where {N}
    pb0 = get_initial_pbox(sol)
    return failure_probability(sol, t, _to_lazyset(U), pb0; kwargs...)
end

function failure_probability(sol::ReachSolution, t, U, pb::PBOX; kwargs...)
    return failure_probability(array(sol), t, _to_lazyset(U), pb; kwargs...)
end

# one-dimensional case
function failure_probability(F::AbstractVector{TaylorModelReachSet{N}}, t, U, pb::pbox;
                             normalize=true, maxitive=false, cons=false) where {N}
    # transform initial p-box into a random set
    pbn = normalize ? ProbabilisticReachability.normalize(pb) : pb

    if cons
        focal_elems, masses = cons_split(pbn)
    else
        focal_elems, masses = split(pbn)
    end

    # remove those reach-sets that don't intersect time t
    idx = findall(R -> _check_time_intersection(t, R), F)

    if isempty(idx)
        return interval(0)

    elseif length(idx) == 1
        return _failure_probability(F[idx[1]], t, U, focal_elems, masses)

    elseif any(issubset.(set.(box_approximation.(F[idx])), Ref(U)))
        return interval(1)

    elseif maxitive
        failure_probs = [_failure_probability(F[ids], narrow(t ∩ tspan(F[ids])), U, focal_elems,
                                              masses) for ids in idx]

        return maximum(failure_probs)
    else
        return _failure_probability_common(view(F, idx), t, U, focal_elems, masses)
    end
end

# n-dimensional case
function failure_probability(F::AbstractVector{TaylorModelReachSet{N}}, t, U, pb::nPbox;
                             normalize=true, maxitive=false, cons=false) where {N}
    # transform initial p-box into a random set
    pbn = normalize ? ProbabilisticReachability.normalize(pb) : pb

    if cons
        focal_elems, masses = cons_split(pbn)
    else
        focal_elems, masses = split(pbn)
    end

    # remove those reach-sets that don't intersect time t
    idx = findall(R -> _check_time_intersection(t, R), F)

    if isempty(idx)
        return interval(0)

    elseif length(idx) == 1
        return _failure_probability(F[idx[1]], t, U, focal_elems, masses)

    elseif !(U isa AbstractVector) && any(issubset.(box_approximation.(F[idx]), Ref(U)))
        # algorithm already computes the box approx, better to use it
        return interval(1)
    elseif maxitive
        failure_probs = [_failure_probability(F[ids], narrow(t ∩ tspan(F[ids])), U, focal_elems,
                                              masses) for ids in idx]

        return maximum(failure_probs)
    else
        return _failure_probability_common(view(F, idx), t, U, focal_elems, masses)
    end
end

function _failure_probability_common(F::AbstractVector{<:TaylorModelReachSet}, t, U, focal_elems,
                                     masses; ST=Zonotope)

    # precompute interesting R's with non-empty intersection in space
    idx = Int[]
    tint = IntervalArithmetic.Interval[]
    for (j, R) in enumerate(F)
        tj = t ∩ tspan(R)
        x = overapproximate(R, ST; Δt=tj)
        if !isdisjoint(x, U)
            push!(idx, j)
            push!(tint, tj)
        end
    end
    F2 = view(F, idx)

    upperPB = 0.0
    lowerPB = 0.0

    for i in 1:length(focal_elems)
        fe_i = focal_elems[i]

        intersects = false

        for (j, R) in enumerate(F2)
            tj = tint[j]
            x = overapproximate(R, ST; Δt=tj, dom=fe_i)

            # does the reach-set intersect the unsafe states?
            if !intersects && !isdisjoint(x, U)
                upperPB += masses[i]
                intersects = true
            end

            # is the interval included in the unsafe states?
            if intersects && x ⊆ U
                lowerPB += masses[i]
                break
            end
        end
    end
    bounds = interval(lowerPB, upperPB)

    return bounds
end

# more general version with multiple properties (interpreted as the conjunction)
function _failure_probability_common(F::AbstractVector{<:TaylorModelReachSet},
                                     ts::AbstractVector, Us::AbstractVector, focal_elems, masses;
                                     ST=Zonotope)

    # precompute interesting R's with non-empty intersection in time and space
    m = length(ts)
    @assert m == length(Us) "incompatible lengths of $ts and $Us"
    Fs = Vector{typeof(view(F, [1]))}(undef, m)
    tints = Vector{Vector{IntervalArithmetic.Interval}}(undef, m)
    @inbounds for i in 1:m
        idx = Int[]
        tint = IntervalArithmetic.Interval[]
        U = Us[i]
        t = ts[i]
        for (j, R) in enumerate(F)
            tj = t ∩ tspan(R)
            if isempty(tj)
                continue
            end
            x = overapproximate(R, ST; Δt=tj)
            if !isdisjoint(x, U)
                push!(idx, j)
                push!(tint, tj)
            end
        end
        Fs[i] = view(F, idx)
        tints[i] = tint
    end

    upperPB = 0.0
    lowerPB = 0.0

    for i in 1:length(focal_elems)
        fe_i = focal_elems[i]

        intersects = false

        @inbounds for i in 1:m
            U = Us[i]
            tint = tints[i]
            F2 = Fs[i]
            included = false

            for (j, R) in enumerate(F2)
                tj = tint[j]
                x = overapproximate(R, ST; Δt=tj, dom=fe_i)

                # does the reach set intersect the failure domain?
                if !intersects && !isdisjoint(x, U)
                    upperPB += masses[i]
                    intersects = true
                end

                # is the reach set included in the failure domain?
                if intersects && x ⊆ U
                    lowerPB += masses[i]
                    included = true
                    break
                end
            end

            if included
                break
            end
        end
    end
    bounds = interval(lowerPB, upperPB)

    return bounds
end

# =========================================
# Evaluate function for analytic solution
# =========================================

# assumes f = f(t, x)
function _failure_probability_analytic(f, t, U, focal_elems, masses)
    upperPB = 0.0
    lowerPB = 0.0

    @inbounds for i in 1:length(focal_elems)
        x = f(t, focal_elems[i])

        # does the reach-set intersect the unsafe states?
        if !isdisjoint(x, U)
            upperPB += masses[i]

            # is the interval included in the unsafe states?
            if x ⊆ U
                lowerPB += masses[i]
            end
        end
    end
    bounds = interval(lowerPB, upperPB)

    return bounds
end
